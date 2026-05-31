/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/external/hdf5.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/molecular_weights.hpp>
#include <gauxc/runtime_environment.hpp>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

#ifdef GAUXC_BENCHMARK_HAS_CUEST
#include <cuest.h>
#include <cuda_runtime_api.h>
#include <integratorxx/quadratures/radial/becke.hpp>
#include <integratorxx/quadratures/radial/mhl.hpp>
#include <integratorxx/quadratures/radial/muraknowles.hpp>
#include <integratorxx/quadratures/radial/treutlerahlrichs.hpp>
#endif

#include <highfive/H5File.hpp>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_NO_CUDA
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using matrix_type = Eigen::MatrixXd;

struct Options {
  std::string implementation = "native-gauxc";
  std::string category = "XC";
  std::string quantity = "VXC";
  std::string file;
  int repeats = 3;
  int warmup = 1;
  std::size_t batch_size = 512;
  double basis_tol = 1e-10;
  std::string grid = "GM3";
  std::string pruning = "UNPRUNED";
  std::string rad_quad = "MURAKNOWLES";
  std::string exec_space = "Device";
  std::string functional = "BLYP";
  std::string nlc_math_mode = "FP64";
  std::string gradient_mode = "HF";
  std::string integrator_kernel = "Default";
  std::string lwd_kernel = "Default";
  std::string reduction_kernel = "Default";
  bool force_rks = true;
};

struct Result {
  std::string status = "ok";
  std::string status_detail;
  std::string mode = "RKS";
  int natoms = 0;
  int nbf = 0;
  std::vector<double> seconds;
  double value = 0.0;
  double norm = 0.0;
};

void uppercase(std::string& value) {
  std::transform(value.begin(), value.end(), value.begin(), ::toupper);
}

std::string csv_escape(std::string value) {
  const bool needs_quotes = value.find_first_of(",\n\"") != std::string::npos;
  if(!needs_quotes) return value;
  std::string out = "\"";
  for(char c : value) {
    if(c == '"') out += "\"\"";
    else out += c;
  }
  out += "\"";
  return out;
}

std::string basename(std::string path) {
  const auto pos = path.find_last_of("/");
  return pos == std::string::npos ? path : path.substr(pos + 1);
}

Options parse_options(int argc, char** argv) {
  Options options;
  for(int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto require_value = [&](const char* name) -> std::string {
      if(i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + name);
      return argv[++i];
    };

    if(arg == "--implementation") options.implementation = require_value("--implementation");
    else if(arg == "--category") options.category = require_value("--category");
    else if(arg == "--quantity") options.quantity = require_value("--quantity");
    else if(arg == "--file") options.file = require_value("--file");
    else if(arg == "--repeats") options.repeats = std::stoi(require_value("--repeats"));
    else if(arg == "--warmup") options.warmup = std::stoi(require_value("--warmup"));
    else if(arg == "--batch-size") options.batch_size = static_cast<std::size_t>(std::stoull(require_value("--batch-size")));
    else if(arg == "--basis-tol") options.basis_tol = std::stod(require_value("--basis-tol"));
    else if(arg == "--grid") options.grid = require_value("--grid");
    else if(arg == "--pruning") options.pruning = require_value("--pruning");
    else if(arg == "--rad-quad") options.rad_quad = require_value("--rad-quad");
    else if(arg == "--exec-space") options.exec_space = require_value("--exec-space");
    else if(arg == "--functional") options.functional = require_value("--functional");
    else if(arg == "--nlc-math-mode") options.nlc_math_mode = require_value("--nlc-math-mode");
    else if(arg == "--gradient-mode") options.gradient_mode = require_value("--gradient-mode");
    else if(arg == "--integrator-kernel") options.integrator_kernel = require_value("--integrator-kernel");
    else if(arg == "--lwd-kernel") options.lwd_kernel = require_value("--lwd-kernel");
    else if(arg == "--reduction-kernel") options.reduction_kernel = require_value("--reduction-kernel");
    else if(arg == "--mode") {
      const auto mode = require_value("--mode");
      options.force_rks = mode != "UKS" and mode != "uks";
    } else {
      throw std::runtime_error("Unknown option: " + arg);
    }
  }

  uppercase(options.category);
  uppercase(options.quantity);
  uppercase(options.grid);
  uppercase(options.pruning);
  uppercase(options.rad_quad);
  uppercase(options.exec_space);
  uppercase(options.functional);
  uppercase(options.nlc_math_mode);
  uppercase(options.gradient_mode);

  if(options.file.empty()) throw std::runtime_error("--file is required");
  if(options.repeats <= 0) throw std::runtime_error("--repeats must be positive");
  if(options.warmup < 0) throw std::runtime_error("--warmup must be non-negative");
  return options;
}

matrix_type read_matrix(HighFive::File& file, const std::string& name) {
  auto dset = file.getDataSet(name);
  auto dims = dset.getDimensions();
  if(dims.size() != 2 or dims[0] != dims[1]) {
    throw std::runtime_error("Expected square matrix dataset: " + name);
  }
  matrix_type matrix(dims[0], dims[1]);
  dset.read(matrix.data());
  return matrix;
}

std::vector<double> read_row_major_dataset(HighFive::File& file, const std::string& name, std::size_t& rows, std::size_t& cols) {
  auto dset = file.getDataSet(name);
  auto dims = dset.getDimensions();
  if(dims.size() != 2) throw std::runtime_error("Expected rank-2 matrix dataset: " + name);
  rows = dims[0];
  cols = dims[1];
  std::vector<std::vector<double>> nested;
  dset.read(nested);
  if(nested.size() != rows) throw std::runtime_error("Unexpected row count while reading " + name);
  std::vector<double> row_major(rows * cols);
  for(std::size_t i = 0; i < rows; ++i) {
    if(nested[i].size() != cols) throw std::runtime_error("Unexpected column count while reading " + name);
    for(std::size_t j = 0; j < cols; ++j) {
      row_major[i * cols + j] = nested[i][j];
    }
  }
  return row_major;
}

matrix_type density_from_row_major_coefficients(const std::vector<double>& coeff, std::size_t nocc, std::size_t nbf, double occupation_factor) {
  matrix_type density = matrix_type::Zero(nbf, nbf);
  for(std::size_t occ = 0; occ < nocc; ++occ) {
    const auto* row = coeff.data() + occ * nbf;
    for(std::size_t mu = 0; mu < nbf; ++mu) {
      for(std::size_t nu = 0; nu < nbf; ++nu) {
        density(mu, nu) += occupation_factor * row[mu] * row[nu];
      }
    }
  }
  return density;
}

double mpi_max_seconds(double local_seconds) {
#ifdef GAUXC_HAS_MPI
  double global_seconds = 0.0;
  MPI_Reduce(&local_seconds, &global_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  return global_seconds;
#else
  return local_seconds;
#endif
}

template <typename F>
double time_call(F&& func) {
#ifdef GAUXC_HAS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto start = std::chrono::high_resolution_clock::now();
  func();
#ifdef GAUXC_HAS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto stop = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(stop - start).count();
}

void summarize_seconds(const std::vector<double>& seconds, double& min_s, double& median_s, double& mean_s) {
  if(seconds.empty()) {
    min_s = median_s = mean_s = 0.0;
    return;
  }
  auto sorted = seconds;
  std::sort(sorted.begin(), sorted.end());
  min_s = sorted.front();
  median_s = sorted[sorted.size() / 2];
  mean_s = std::accumulate(sorted.begin(), sorted.end(), 0.0) / static_cast<double>(sorted.size());
}

void emit_csv(const Options& options, const Result& result) {
  double min_s = 0.0, median_s = 0.0, mean_s = 0.0;
  summarize_seconds(result.seconds, min_s, median_s, mean_s);
  std::cout << "system,file,implementation,category,quantity,status,status_detail,mode,natoms,nbf,repeats,warmup,min_s,median_s,mean_s,value,norm\n";
  std::cout << csv_escape(basename(options.file)) << ','
            << csv_escape(options.file) << ','
            << csv_escape(options.implementation) << ','
            << csv_escape(options.category) << ','
            << csv_escape(options.quantity) << ','
            << csv_escape(result.status) << ','
            << csv_escape(result.status_detail) << ','
            << csv_escape(result.mode) << ','
            << result.natoms << ','
            << result.nbf << ','
            << options.repeats << ','
            << options.warmup << ','
            << std::scientific << std::setprecision(12)
            << min_s << ','
            << median_s << ','
            << mean_s << ','
            << result.value << ','
            << result.norm << '\n';
}

ExchCXX::Functional functional_from_string(const std::string& functional) {
  if(functional == "BLYP") return ExchCXX::Functional::BLYP;
  if(functional == "PBE") return ExchCXX::Functional::PBE;
  if(functional == "PBE0") return ExchCXX::Functional::PBE0;
  if(functional == "SVWN5") return ExchCXX::Functional::SVWN5;
  throw std::runtime_error("Unsupported benchmark functional: " + functional);
}

Result metadata_only_result(const Options& options, const std::string& status, const std::string& detail);


#ifdef GAUXC_BENCHMARK_HAS_CUEST

void check_cuda(cudaError_t status, const char* context) {
  if(status != cudaSuccess) {
    throw std::runtime_error(std::string(context) + ": " + cudaGetErrorString(status));
  }
}

void check_cuest(cuestStatus_t status, const char* context) {
  if(status != CUEST_STATUS_SUCCESS) {
    throw std::runtime_error(std::string(context) + " failed with cuEST status " + std::to_string(static_cast<int>(status)));
  }
}

struct DeviceBuffer {
  void* ptr = nullptr;
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t bytes) { reset(bytes); }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  DeviceBuffer(DeviceBuffer&& other) noexcept : ptr(other.ptr) { other.ptr = nullptr; }
  DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
    if(this != &other) {
      if(ptr) cudaFree(ptr);
      ptr = other.ptr;
      other.ptr = nullptr;
    }
    return *this;
  }
  ~DeviceBuffer() { if(ptr) cudaFree(ptr); }
  void reset(std::size_t bytes) {
    if(ptr) cudaFree(ptr);
    ptr = nullptr;
    if(bytes) check_cuda(cudaMalloc(&ptr, bytes), "cudaMalloc");
  }
  template <typename T> T* as() { return static_cast<T*>(ptr); }
};

struct WorkspaceHolder {
  std::vector<unsigned char> host;
  DeviceBuffer device;
  cuestWorkspace_t workspace{};

  explicit WorkspaceHolder(const cuestWorkspaceDescriptor_t& descriptor) :
    host(descriptor.hostBufferSizeInBytes), device(descriptor.deviceBufferSizeInBytes) {
    workspace.hostBuffer = reinterpret_cast<uintptr_t>(host.data());
    workspace.hostBufferSizeInBytes = host.size();
    workspace.deviceBuffer = reinterpret_cast<uintptr_t>(device.ptr);
    workspace.deviceBufferSizeInBytes = descriptor.deviceBufferSizeInBytes;
  }
};

struct CuestHandleDeleter { void operator()(void* handle) const { if(handle) cuestDestroy(handle); } };
struct CuestParametersDeleter {
  cuestParametersType_t type;
  void operator()(void* handle) const { if(handle) cuestParametersDestroy(type, handle); }
};
struct CuestShellDeleter { void operator()(void* handle) const { if(handle) cuestAOShellDestroy(handle); } };
struct CuestBasisDeleter { void operator()(void* handle) const { if(handle) cuestAOBasisDestroy(handle); } };
struct CuestAtomGridDeleter { void operator()(void* handle) const { if(handle) cuestAtomGridDestroy(handle); } };
struct CuestMolecularGridDeleter { void operator()(void* handle) const { if(handle) cuestMolecularGridDestroy(handle); } };
struct CuestXCPlanDeleter { void operator()(void* handle) const { if(handle) cuestXCIntPlanDestroy(handle); } };

template <typename Deleter>
using unique_cuest = std::unique_ptr<void, Deleter>;

unique_cuest<CuestParametersDeleter> make_parameters(cuestParametersType_t type) {
  void* raw = nullptr;
  check_cuest(cuestParametersCreate(type, &raw), "cuestParametersCreate");
  return unique_cuest<CuestParametersDeleter>(raw, CuestParametersDeleter{type});
}

std::vector<double> factor_density_to_row_major_coefficients(const matrix_type& density, bool rks_total_density) {
  const auto scaled = rks_total_density ? 0.5 * density : density;
  Eigen::SelfAdjointEigenSolver<matrix_type> solver(scaled);
  if(solver.info() != Eigen::Success) throw std::runtime_error("density eigensolve failed");
  const auto& values = solver.eigenvalues();
  const auto& vectors = solver.eigenvectors();
  std::vector<int> active;
  for(int i = 0; i < values.size(); ++i) {
    if(values(i) > 1e-10) active.push_back(i);
  }
  if(active.empty()) throw std::runtime_error("density factorization produced no occupied vectors");
  std::vector<double> coeff(active.size() * density.rows());
  for(std::size_t row = 0; row < active.size(); ++row) {
    const int idx = active[row];
    const double scale = std::sqrt(values(idx));
    for(int col = 0; col < density.rows(); ++col) {
      coeff[row * density.rows() + col] = scale * vectors(col, idx);
    }
  }
  return coeff;
}

std::vector<double> read_row_major_matrix(HighFive::File& file, const std::string& name, std::size_t& rows, std::size_t& cols) {
  return read_row_major_dataset(file, name, rows, cols);
}

template <typename RadialType>
void fill_radial(std::int64_t nrad, double scale, std::vector<double>& nodes, std::vector<double>& weights) {
  RadialType radial(nrad, scale);
  nodes.resize(radial.npts());
  weights.resize(radial.npts());
  for(std::size_t i = 0; i < radial.npts(); ++i) {
    nodes[i] = radial.points(i);
    weights[i] = radial.weights(i);
  }
}

void fill_radial(GauXC::RadialQuad radial_quad, std::int64_t nrad, double scale,
                 std::vector<double>& nodes, std::vector<double>& weights) {
  switch(radial_quad) {
    case GauXC::RadialQuad::Becke:
      fill_radial<IntegratorXX::Becke<double, double>>(nrad, scale, nodes, weights);
      break;
    case GauXC::RadialQuad::MuraKnowles:
      fill_radial<IntegratorXX::MuraKnowles<double, double>>(nrad, scale, nodes, weights);
      break;
    case GauXC::RadialQuad::MurrayHandyLaming:
      fill_radial<IntegratorXX::MurrayHandyLaming<double, double>>(nrad, scale, nodes, weights);
      break;
    case GauXC::RadialQuad::TreutlerAhlrichs:
      fill_radial<IntegratorXX::TreutlerAhlrichs<double, double>>(nrad, scale, nodes, weights);
      break;
    default:
      throw std::runtime_error("unsupported radial quadrature for cuEST benchmark");
  }
}

cuestXCIntPlanParametersFunctional_t cuest_functional_from_string(const std::string& functional) {
  if(functional == "BLYP") return CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_BLYP;
  if(functional == "PBE") return CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE;
  if(functional == "PBE0") return CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE0;
  if(functional == "SVWN5") return CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_SVWN5;
  throw std::runtime_error("unsupported cuEST benchmark functional: " + functional);
}

std::vector<std::uint64_t> shells_per_atom_from_basis(const GauXC::Molecule& mol, const GauXC::BasisSet<double>& basis) {
  std::vector<std::uint64_t> shells_per_atom(mol.natoms(), 0);
  std::size_t atom_index = 0;
  for(const auto& shell : basis) {
    const double* origin = shell.O_data();
    while(atom_index < mol.natoms()) {
      const auto& atom = mol[atom_index];
      const double dx = origin[0] - atom.x;
      const double dy = origin[1] - atom.y;
      const double dz = origin[2] - atom.z;
      if(std::sqrt(dx * dx + dy * dy + dz * dz) < 1e-8) break;
      ++atom_index;
    }
    if(atom_index >= mol.natoms()) throw std::runtime_error("could not assign basis shell to atom");
    shells_per_atom[atom_index] += 1;
  }
  return shells_per_atom;
}

template <typename T>
DeviceBuffer copy_to_device(const std::vector<T>& host) {
  DeviceBuffer device(host.size() * sizeof(T));
  check_cuda(cudaMemcpy(device.ptr, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy H2D");
  return device;
}

std::vector<double> copy_from_device(DeviceBuffer& device, std::size_t count) {
  std::vector<double> host(count);
  check_cuda(cudaMemcpy(host.data(), device.ptr, count * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
  return host;
}

template <typename RuntimeIgnored>
Result run_cuest_native(const Options& options) {
  using namespace GauXC;

  if(options.category == "EXCHANGE") {
    return metadata_only_result(options, "missing_input", "cuEST DF exchange requires an auxiliary basis and AO pair list not present in current fixtures");
  }

  Result result;

  std::map<std::string, AtomicGridSizeDefault> grid_map = {
    {"FINE", AtomicGridSizeDefault::FineGrid},
    {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
    {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid},
    {"GM3", AtomicGridSizeDefault::GM3},
    {"GM5", AtomicGridSizeDefault::GM5}
  };
  std::map<std::string, RadialQuad> rad_quad_map = {
    {"BECKE", RadialQuad::Becke},
    {"MURAKNOWLES", RadialQuad::MuraKnowles},
    {"MK", RadialQuad::MuraKnowles},
    {"TREUTLERAHLRICHS", RadialQuad::TreutlerAhlrichs},
    {"TA", RadialQuad::TreutlerAhlrichs},
    {"MURRAYHANDYLAMING", RadialQuad::MurrayHandyLaming},
    {"MHL", RadialQuad::MurrayHandyLaming}
  };
  const auto radial_quad = rad_quad_map.at(options.rad_quad);
  const auto grid_size = grid_map.at(options.grid);

  Molecule mol;
  BasisSet<double> basis;
  read_hdf5_record(mol, options.file, "/MOLECULE");
  read_hdf5_record(basis, options.file, "/BASIS");
  result.natoms = static_cast<int>(mol.natoms());
  result.nbf = static_cast<int>(basis.nbf());
  result.mode = options.force_rks ? "RKS" : "UKS";

  HighFive::File file(options.file, HighFive::File::ReadOnly);
  const bool has_spin_density = file.exist("/DENSITY_Z");
  auto P = read_matrix(file, has_spin_density ? "/DENSITY_SCALAR" : "/DENSITY");
  matrix_type Pz;
  if(!options.force_rks && has_spin_density) Pz = read_matrix(file, "/DENSITY_Z");

  auto handle_params = make_parameters(CUEST_HANDLE_PARAMETERS);
  cuestHandle_t raw_handle = nullptr;
  check_cuest(cuestCreate(handle_params.get(), &raw_handle), "cuestCreate");
  unique_cuest<CuestHandleDeleter> handle(raw_handle);

  auto shell_params = make_parameters(CUEST_AOSHELL_PARAMETERS);
  std::vector<unique_cuest<CuestShellDeleter>> shell_storage;
  std::vector<cuestAOShell_t> shell_handles;
  shell_storage.reserve(basis.size());
  shell_handles.reserve(basis.size());
  for(const auto& shell : basis) {
    cuestAOShell_t raw_shell = nullptr;
    const int32_t is_pure = shell.pure() || shell.l() < 2;
    check_cuest(cuestAOShellCreate(handle.get(), is_pure, shell.l(), shell.nprim(),
      shell.alpha_data(), shell.coeff_data(), shell_params.get(), &raw_shell), "cuestAOShellCreate");
    shell_storage.emplace_back(raw_shell);
    shell_handles.push_back(raw_shell);
  }

  const auto shells_per_atom = shells_per_atom_from_basis(mol, basis);
  auto basis_params = make_parameters(CUEST_AOBASIS_PARAMETERS);
  cuestWorkspaceDescriptor_t basis_persistent_desc{}, basis_temp_desc{};
  check_cuest(cuestAOBasisCreateWorkspaceQuery(handle.get(), mol.natoms(), shells_per_atom.data(),
    shell_handles.data(), basis_params.get(), &basis_persistent_desc, &basis_temp_desc, nullptr),
    "cuestAOBasisCreateWorkspaceQuery");
  WorkspaceHolder basis_persistent(basis_persistent_desc);
  WorkspaceHolder basis_temp(basis_temp_desc);
  cuestAOBasis_t raw_basis = nullptr;
  check_cuest(cuestAOBasisCreate(handle.get(), mol.natoms(), shells_per_atom.data(), shell_handles.data(),
    basis_params.get(), &basis_persistent.workspace, &basis_temp.workspace, &raw_basis), "cuestAOBasisCreate");
  unique_cuest<CuestBasisDeleter> cuest_basis(raw_basis);

  auto atom_grid_params = make_parameters(CUEST_ATOMGRID_PARAMETERS);
  std::vector<unique_cuest<CuestAtomGridDeleter>> atom_grid_storage;
  std::vector<cuestAtomGrid_t> atom_grid_handles;
  std::vector<std::vector<double>> radial_nodes_storage;
  std::vector<std::vector<double>> radial_weights_storage;
  std::vector<std::vector<std::uint64_t>> angular_storage;
  atom_grid_storage.reserve(mol.natoms());
  atom_grid_handles.reserve(mol.natoms());
  radial_nodes_storage.reserve(mol.natoms());
  radial_weights_storage.reserve(mol.natoms());
  angular_storage.reserve(mol.natoms());
  for(const auto& atom : mol) {
    auto [radial_size, angular_size] = default_grid_size(atom.Z, radial_quad, grid_size);
    radial_nodes_storage.emplace_back();
    radial_weights_storage.emplace_back();
    fill_radial(radial_quad, radial_size.get(), default_radial_scaling_factor(radial_quad, atom.Z).get(),
      radial_nodes_storage.back(), radial_weights_storage.back());
    angular_storage.emplace_back(radial_nodes_storage.back().size(), static_cast<std::uint64_t>(angular_size.get()));
    cuestAtomGrid_t raw_grid = nullptr;
    check_cuest(cuestAtomGridCreate(handle.get(), radial_nodes_storage.back().size(),
      radial_nodes_storage.back().data(), radial_weights_storage.back().data(), angular_storage.back().data(),
      atom_grid_params.get(), &raw_grid), "cuestAtomGridCreate");
    atom_grid_storage.emplace_back(raw_grid);
    atom_grid_handles.push_back(raw_grid);
  }

  std::vector<double> xyz;
  xyz.reserve(3 * mol.natoms());
  for(const auto& atom : mol) {
    xyz.push_back(atom.x);
    xyz.push_back(atom.y);
    xyz.push_back(atom.z);
  }

  auto molecular_grid_params = make_parameters(CUEST_MOLECULARGRID_PARAMETERS);
  cuestWorkspaceDescriptor_t grid_persistent_desc{}, grid_temp_desc{};
  check_cuest(cuestMolecularGridCreateWorkspaceQuery(handle.get(), mol.natoms(), atom_grid_handles.data(), xyz.data(),
    molecular_grid_params.get(), &grid_persistent_desc, &grid_temp_desc, nullptr),
    "cuestMolecularGridCreateWorkspaceQuery");
  WorkspaceHolder grid_persistent(grid_persistent_desc);
  WorkspaceHolder grid_temp(grid_temp_desc);
  cuestMolecularGrid_t raw_molecular_grid = nullptr;
  check_cuest(cuestMolecularGridCreate(handle.get(), mol.natoms(), atom_grid_handles.data(), xyz.data(),
    molecular_grid_params.get(), &grid_persistent.workspace, &grid_temp.workspace, &raw_molecular_grid),
    "cuestMolecularGridCreate");
  unique_cuest<CuestMolecularGridDeleter> molecular_grid(raw_molecular_grid);

  auto xc_plan_params = make_parameters(CUEST_XCINTPLAN_PARAMETERS);
  cuestWorkspaceDescriptor_t plan_persistent_desc{}, plan_temp_desc{};
  check_cuest(cuestXCIntPlanCreateWorkspaceQuery(handle.get(), cuest_basis.get(), molecular_grid.get(),
    cuest_functional_from_string(options.functional), xc_plan_params.get(), &plan_persistent_desc, &plan_temp_desc, nullptr),
    "cuestXCIntPlanCreateWorkspaceQuery");
  WorkspaceHolder plan_persistent(plan_persistent_desc);
  WorkspaceHolder plan_temp(plan_temp_desc);
  cuestXCIntPlan_t raw_plan = nullptr;
  check_cuest(cuestXCIntPlanCreate(handle.get(), cuest_basis.get(), molecular_grid.get(),
    cuest_functional_from_string(options.functional), xc_plan_params.get(), &plan_persistent.workspace, &plan_temp.workspace, &raw_plan),
    "cuestXCIntPlanCreate");
  unique_cuest<CuestXCPlanDeleter> xc_plan(raw_plan);

  std::vector<double> coeff_alpha;
  std::vector<double> coeff_beta;
  std::uint64_t nocc_alpha = 0;
  std::uint64_t nocc_beta = 0;
  if(options.force_rks && file.exist("/CUEST_MO_COEFF_OCC")) {
    std::size_t rows = 0, cols = 0;
    coeff_alpha = read_row_major_matrix(file, "/CUEST_MO_COEFF_OCC", rows, cols);
    if(cols != basis.nbf()) throw std::runtime_error("CUEST_MO_COEFF_OCC column count does not match basis size");
    nocc_alpha = rows;
  } else if(!options.force_rks && file.exist("/CUEST_MO_COEFF_OCC_A") && file.exist("/CUEST_MO_COEFF_OCC_B")) {
    std::size_t rows_a = 0, cols_a = 0, rows_b = 0, cols_b = 0;
    coeff_alpha = read_row_major_matrix(file, "/CUEST_MO_COEFF_OCC_A", rows_a, cols_a);
    coeff_beta = read_row_major_matrix(file, "/CUEST_MO_COEFF_OCC_B", rows_b, cols_b);
    if(cols_a != basis.nbf() || cols_b != basis.nbf()) throw std::runtime_error("CUEST_MO_COEFF_OCC spin column count does not match basis size");
    nocc_alpha = rows_a;
    nocc_beta = rows_b;
  } else if(options.force_rks || !has_spin_density) {
    coeff_alpha = factor_density_to_row_major_coefficients(P, true);
    nocc_alpha = coeff_alpha.size() / basis.nbf();
  } else {
    coeff_alpha = factor_density_to_row_major_coefficients(0.5 * (P + Pz), false);
    coeff_beta = factor_density_to_row_major_coefficients(0.5 * (P - Pz), false);
    nocc_alpha = coeff_alpha.size() / basis.nbf();
    nocc_beta = coeff_beta.size() / basis.nbf();
  }
  DeviceBuffer coeff_alpha_device = copy_to_device(coeff_alpha);
  DeviceBuffer coeff_beta_device;
  if(!coeff_beta.empty()) coeff_beta_device = copy_to_device(coeff_beta);

  cuestWorkspaceDescriptor_t variable_buffer{};
  variable_buffer.deviceBufferSizeInBytes = 2ull * 1024ull * 1024ull * 1024ull;

  auto run_with_workspace = [&](auto query, auto compute, std::size_t output_count, double* energy) {
    cuestWorkspaceDescriptor_t temp_desc{};
    check_cuest(query(temp_desc), "cuEST compute workspace query");
    WorkspaceHolder temp(temp_desc);
    DeviceBuffer output_device(output_count * sizeof(double));
    auto call = [&]() {
      check_cuest(compute(temp.workspace, output_device), "cuEST compute");
      check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
      auto output = copy_from_device(output_device, output_count);
      result.norm = std::sqrt(std::inner_product(output.begin(), output.end(), output.begin(), 0.0));
      if(energy) result.value = *energy;
      else result.value = 0.0;
    };
    for(int i = 0; i < options.warmup; ++i) call();
    result.seconds.reserve(options.repeats);
    for(int i = 0; i < options.repeats; ++i) {
      result.seconds.push_back(time_call(call));
    }
  };

  if(options.category == "XC") {
    if(options.quantity == "EXC" || options.quantity == "VXC") {
      auto compute_params = make_parameters(CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS);
      double energy = 0.0;
      run_with_workspace(
        [&](cuestWorkspaceDescriptor_t& temp_desc) {
          return cuestXCPotentialRKSComputeWorkspaceQuery(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp_desc, nocc_alpha, coeff_alpha_device.as<double>(), &energy, nullptr);
        },
        [&](cuestWorkspace_t& temp, DeviceBuffer& output) {
          return cuestXCPotentialRKSCompute(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp, nocc_alpha, coeff_alpha_device.as<double>(), &energy, output.as<double>());
        },
        basis.nbf() * basis.nbf(), &energy);
      if(options.quantity == "EXC") result.norm = 0.0;
    } else if(options.quantity == "EXC_GRAD") {
      auto compute_params = make_parameters(CUEST_XCDERIVATIVERKSCOMPUTE_PARAMETERS);
      run_with_workspace(
        [&](cuestWorkspaceDescriptor_t& temp_desc) {
          return cuestXCDerivativeRKSComputeWorkspaceQuery(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp_desc, nocc_alpha, coeff_alpha_device.as<double>(), nullptr);
        },
        [&](cuestWorkspace_t& temp, DeviceBuffer& output) {
          return cuestXCDerivativeRKSCompute(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp, nocc_alpha, coeff_alpha_device.as<double>(), output.as<double>());
        },
        3 * mol.natoms(), nullptr);
    } else {
      throw std::runtime_error("unsupported cuEST XC quantity: " + options.quantity);
    }
  } else if(options.category == "NLC") {
    if(options.quantity == "NLC" || options.quantity == "VNLC") {
      auto compute_params = make_parameters(CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS);
      const double scale = 1.0;
      const double vv10_c = 0.0093;
      const double vv10_b = 6.3;
      check_cuest(cuestParametersConfigure(CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS, compute_params.get(),
        CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_SCALE, &scale, sizeof(double)), "configure VV10 scale");
      check_cuest(cuestParametersConfigure(CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS, compute_params.get(),
        CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_C, &vv10_c, sizeof(double)), "configure VV10 c");
      check_cuest(cuestParametersConfigure(CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS, compute_params.get(),
        CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B, &vv10_b, sizeof(double)), "configure VV10 b");
      double energy = 0.0;
      run_with_workspace(
        [&](cuestWorkspaceDescriptor_t& temp_desc) {
          return cuestNonlocalXCPotentialRKSComputeWorkspaceQuery(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp_desc, nocc_alpha, coeff_alpha_device.as<double>(), &energy, nullptr);
        },
        [&](cuestWorkspace_t& temp, DeviceBuffer& output) {
          return cuestNonlocalXCPotentialRKSCompute(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp, nocc_alpha, coeff_alpha_device.as<double>(), &energy, output.as<double>());
        },
        basis.nbf() * basis.nbf(), &energy);
      if(options.quantity == "NLC") result.norm = 0.0;
    } else if(options.quantity == "NLC_GRAD") {
      auto compute_params = make_parameters(CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS);
      const double scale = 1.0;
      const double vv10_c = 0.0093;
      const double vv10_b = 6.3;
      check_cuest(cuestParametersConfigure(CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS, compute_params.get(),
        CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_VV10_SCALE, &scale, sizeof(double)), "configure VV10 grad scale");
      check_cuest(cuestParametersConfigure(CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS, compute_params.get(),
        CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_VV10_C, &vv10_c, sizeof(double)), "configure VV10 grad c");
      check_cuest(cuestParametersConfigure(CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS, compute_params.get(),
        CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_VV10_B, &vv10_b, sizeof(double)), "configure VV10 grad b");
      run_with_workspace(
        [&](cuestWorkspaceDescriptor_t& temp_desc) {
          return cuestNonlocalXCDerivativeRKSComputeWorkspaceQuery(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp_desc, nocc_alpha, coeff_alpha_device.as<double>(), nullptr);
        },
        [&](cuestWorkspace_t& temp, DeviceBuffer& output) {
          return cuestNonlocalXCDerivativeRKSCompute(handle.get(), xc_plan.get(), compute_params.get(), &variable_buffer,
            &temp, nocc_alpha, coeff_alpha_device.as<double>(), output.as<double>());
        },
        3 * mol.natoms(), nullptr);
    } else {
      throw std::runtime_error("unsupported cuEST NLC quantity: " + options.quantity);
    }
  }

  return result;
}

#endif

template <typename RuntimeType>
Result run_native(const Options& options, RuntimeType& runtime) {
  using namespace GauXC;

  Result result;

  std::map<std::string, AtomicGridSizeDefault> grid_map = {
    {"FINE", AtomicGridSizeDefault::FineGrid},
    {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
    {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid},
    {"GM3", AtomicGridSizeDefault::GM3},
    {"GM5", AtomicGridSizeDefault::GM5}
  };
  std::map<std::string, PruningScheme> pruning_map = {
    {"UNPRUNED", PruningScheme::Unpruned},
    {"ROBUST", PruningScheme::Robust},
    {"TREUTLER", PruningScheme::Treutler}
  };
  std::map<std::string, RadialQuad> rad_quad_map = {
    {"BECKE", RadialQuad::Becke},
    {"MURAKNOWLES", RadialQuad::MuraKnowles},
    {"MK", RadialQuad::MuraKnowles},
    {"TREUTLERAHLRICHS", RadialQuad::TreutlerAhlrichs},
    {"TA", RadialQuad::TreutlerAhlrichs},
    {"MURRAYHANDYLAMING", RadialQuad::MurrayHandyLaming},
    {"MHL", RadialQuad::MurrayHandyLaming}
  };

  Molecule mol;
  BasisSet<double> basis;
  read_hdf5_record(mol, options.file, "/MOLECULE");
  read_hdf5_record(basis, options.file, "/BASIS");
  for(auto& shell : basis) shell.set_shell_tolerance(options.basis_tol);

  HighFive::File file(options.file, HighFive::File::ReadOnly);
  const bool has_spin_density = file.exist("/DENSITY_Z");
  const bool uks = has_spin_density and not options.force_rks;
  result.mode = uks ? "UKS" : "RKS";
  result.natoms = static_cast<int>(mol.natoms());
  result.nbf = static_cast<int>(basis.nbf());

  const auto density_name = has_spin_density ? std::string("/DENSITY_SCALAR") : std::string("/DENSITY");
  auto P = read_matrix(file, density_name);
  matrix_type Pz;
  if(uks) Pz = read_matrix(file, "/DENSITY_Z");

  const auto molgrid = MolGridFactory::create_default_molgrid(
    mol, pruning_map.at(options.pruning), BatchSize(options.batch_size),
    rad_quad_map.at(options.rad_quad), grid_map.at(options.grid));

  const auto lb_exec_space = options.exec_space == "DEVICE" ? ExecutionSpace::Device : ExecutionSpace::Host;
  LoadBalancerFactory lb_factory(lb_exec_space, "Default");
  auto load_balancer = lb_factory.get_shared_instance(runtime, mol, molgrid, basis);

  MolecularWeightsFactory mw_factory(lb_exec_space, "Default", MolecularWeightsSettings{});
  auto molecular_weights = mw_factory.get_instance();
  molecular_weights.modify_weights(*load_balancer);

  auto functional = functional_type(ExchCXX::Backend::builtin, functional_from_string(options.functional),
    uks ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized);
  XCIntegratorFactory<matrix_type> integrator_factory(lb_exec_space, "Replicated",
    options.integrator_kernel, options.lwd_kernel, options.reduction_kernel);
  auto integrator = integrator_factory.get_instance(functional, load_balancer);

  IntegratorSettingsNLC nlc_settings;
  nlc_settings.math_mode = options.nlc_math_mode == "FLOATPAIR" ? NLCMathMode::FloatPair : NLCMathMode::NativeFP64;
  nlc_settings.include_weight_derivatives = options.gradient_mode == "FULL";

  auto evaluate = [&]() {
    if(options.category == "XC") {
      if(options.quantity == "EXC") {
        result.value = uks ? integrator.eval_exc(P, Pz) : integrator.eval_exc(P);
        result.norm = 0.0;
      } else if(options.quantity == "VXC") {
        if(uks) {
          auto out = integrator.eval_exc_vxc(P, Pz);
          result.value = std::get<0>(out);
          result.norm = std::get<1>(out).norm() + std::get<2>(out).norm();
        } else {
          auto out = integrator.eval_exc_vxc(P);
          result.value = std::get<0>(out);
          result.norm = std::get<1>(out).norm();
        }
      } else if(options.quantity == "EXC_GRAD") {
        auto grad = uks ? integrator.eval_exc_grad(P, Pz) : integrator.eval_exc_grad(P);
        result.value = 0.0;
        result.norm = std::sqrt(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0.0));
      } else {
        throw std::runtime_error("Unknown XC quantity: " + options.quantity);
      }
    } else if(options.category == "NLC") {
      if(options.quantity == "NLC") {
        result.value = uks ? integrator.eval_nlc(P, Pz, nlc_settings) : integrator.eval_nlc(P, nlc_settings);
        result.norm = 0.0;
      } else if(options.quantity == "VNLC") {
        if(uks) {
          auto out = integrator.eval_nlc_vnlc(P, Pz, nlc_settings);
          result.value = std::get<0>(out);
          result.norm = std::get<1>(out).norm() + std::get<2>(out).norm();
        } else {
          auto out = integrator.eval_nlc_vnlc(P, nlc_settings);
          result.value = std::get<0>(out);
          result.norm = std::get<1>(out).norm();
        }
      } else if(options.quantity == "NLC_GRAD") {
        auto grad = uks ? integrator.eval_nlc_grad(P, Pz, nlc_settings) : integrator.eval_nlc_grad(P, nlc_settings);
        result.value = 0.0;
        result.norm = std::sqrt(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0.0));
      } else {
        throw std::runtime_error("Unknown NLC quantity: " + options.quantity);
      }
    } else if(options.category == "EXCHANGE") {
      if(options.quantity != "K") throw std::runtime_error("Unknown exchange quantity: " + options.quantity);
      auto K = integrator.eval_exx(P, IntegratorSettingsEXX{});
      result.value = 0.0;
      result.norm = K.norm();
    } else {
      throw std::runtime_error("Unknown category: " + options.category);
    }
  };

  for(int i = 0; i < options.warmup; ++i) evaluate();
  result.seconds.reserve(options.repeats);
  for(int i = 0; i < options.repeats; ++i) {
    const auto local_seconds = time_call(evaluate);
    result.seconds.push_back(mpi_max_seconds(local_seconds));
  }

  return result;
}

template <typename RuntimeType>
Result run_gauxc_cuest_api(const Options& options, RuntimeType& runtime) {
  using namespace GauXC;

  Result result;

  std::map<std::string, AtomicGridSizeDefault> grid_map = {
    {"FINE", AtomicGridSizeDefault::FineGrid},
    {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
    {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid},
    {"GM3", AtomicGridSizeDefault::GM3},
    {"GM5", AtomicGridSizeDefault::GM5}
  };
  std::map<std::string, PruningScheme> pruning_map = {
    {"UNPRUNED", PruningScheme::Unpruned},
    {"ROBUST", PruningScheme::Robust},
    {"TREUTLER", PruningScheme::Treutler}
  };
  std::map<std::string, RadialQuad> rad_quad_map = {
    {"BECKE", RadialQuad::Becke},
    {"MURAKNOWLES", RadialQuad::MuraKnowles},
    {"MK", RadialQuad::MuraKnowles},
    {"TREUTLERAHLRICHS", RadialQuad::TreutlerAhlrichs},
    {"TA", RadialQuad::TreutlerAhlrichs},
    {"MURRAYHANDYLAMING", RadialQuad::MurrayHandyLaming},
    {"MHL", RadialQuad::MurrayHandyLaming}
  };

  Molecule mol;
  BasisSet<double> basis;
  read_hdf5_record(mol, options.file, "/MOLECULE");
  read_hdf5_record(basis, options.file, "/BASIS");
  for(auto& shell : basis) shell.set_shell_tolerance(options.basis_tol);
  result.mode = "RKS";
  result.natoms = static_cast<int>(mol.natoms());
  result.nbf = static_cast<int>(basis.nbf());

  HighFive::File file(options.file, HighFive::File::ReadOnly);
  if(!file.exist("/CUEST_MO_COEFF_OCC")) {
    return metadata_only_result(options, "missing_input", "CUEST_MO_COEFF_OCC is required for the GauXC cuEST API benchmark path");
  }
  std::size_t nocc = 0, nbf = 0;
  const auto coeff = read_row_major_dataset(file, "/CUEST_MO_COEFF_OCC", nocc, nbf);
  if(nbf != basis.nbf()) throw std::runtime_error("CUEST_MO_COEFF_OCC column count does not match basis size");
  const bool density_factorized_input = file.exist("/EXPORT_PROVENANCE");
  const double occupation_factor = density_factorized_input ? 1.0 : 2.0;
  const auto P = density_from_row_major_coefficients(coeff, nocc, nbf, occupation_factor);

  const auto molgrid = MolGridFactory::create_default_molgrid(
    mol, pruning_map.at(options.pruning), BatchSize(options.batch_size),
    rad_quad_map.at(options.rad_quad), grid_map.at(options.grid));

  const auto lb_exec_space = options.exec_space == "DEVICE" ? ExecutionSpace::Device : ExecutionSpace::Host;
  LoadBalancerFactory lb_factory(lb_exec_space, "Default");
  auto load_balancer = lb_factory.get_shared_instance(runtime, mol, molgrid, basis);

  MolecularWeightsFactory mw_factory(lb_exec_space, "Default", MolecularWeightsSettings{});
  auto molecular_weights = mw_factory.get_instance();
  molecular_weights.modify_weights(*load_balancer);

  auto functional = functional_type(ExchCXX::Backend::builtin, functional_from_string(options.functional), ExchCXX::Spin::Unpolarized);
  XCIntegratorFactory<matrix_type> integrator_factory(lb_exec_space, "Replicated",
    options.integrator_kernel, options.lwd_kernel, options.reduction_kernel);
  auto integrator = integrator_factory.get_instance(functional, load_balancer);

  IntegratorSettingsNLC nlc_settings;
  nlc_settings.math_mode = options.nlc_math_mode == "FLOATPAIR" ? NLCMathMode::FloatPair : NLCMathMode::NativeFP64;
  nlc_settings.include_weight_derivatives = options.gradient_mode == "FULL";

  auto evaluate = [&]() {
    if(options.category == "XC") {
      if(options.quantity == "EXC") {
        result.value = integrator.eval_exc(P);
        result.norm = 0.0;
      } else if(options.quantity == "VXC") {
        auto out = integrator.eval_exc_vxc(P);
        result.value = std::get<0>(out);
        result.norm = std::get<1>(out).norm();
      } else if(options.quantity == "EXC_GRAD") {
        auto grad = integrator.eval_exc_grad(P);
        result.value = 0.0;
        result.norm = std::sqrt(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0.0));
      } else {
        throw std::runtime_error("Unknown XC quantity: " + options.quantity);
      }
    } else if(options.category == "NLC") {
      if(options.quantity == "NLC") {
        result.value = integrator.eval_nlc(P, nlc_settings);
        result.norm = 0.0;
      } else if(options.quantity == "VNLC") {
        auto out = integrator.eval_nlc_vnlc(P, nlc_settings);
        result.value = std::get<0>(out);
        result.norm = std::get<1>(out).norm();
      } else if(options.quantity == "NLC_GRAD") {
        auto grad = integrator.eval_nlc_grad(P, nlc_settings);
        result.value = 0.0;
        result.norm = std::sqrt(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0.0));
      } else {
        throw std::runtime_error("Unknown NLC quantity: " + options.quantity);
      }
    } else if(options.category == "EXCHANGE") {
      if(options.quantity != "K") throw std::runtime_error("Unknown exchange quantity: " + options.quantity);
      auto K = integrator.eval_exx(P, IntegratorSettingsEXX{});
      result.value = 0.0;
      result.norm = K.norm();
    } else {
      throw std::runtime_error("Unknown category: " + options.category);
    }
  };

  for(int i = 0; i < options.warmup; ++i) evaluate();
  result.seconds.reserve(options.repeats);
  for(int i = 0; i < options.repeats; ++i) {
    const auto local_seconds = time_call(evaluate);
    result.seconds.push_back(mpi_max_seconds(local_seconds));
  }

  return result;
}

Result metadata_only_result(const Options& options, const std::string& status, const std::string& detail) {
  Result result;
  result.status = status;
  result.status_detail = detail;
  try {
    GauXC::Molecule mol;
    GauXC::BasisSet<double> basis;
    read_hdf5_record(mol, options.file, "/MOLECULE");
    read_hdf5_record(basis, options.file, "/BASIS");
    result.natoms = static_cast<int>(mol.natoms());
    result.nbf = static_cast<int>(basis.nbf());
    HighFive::File file(options.file, HighFive::File::ReadOnly);
    result.mode = file.exist("/DENSITY_Z") && !options.force_rks ? "UKS" : "RKS";
  } catch(const std::exception& ex) {
    result.status_detail += std::string("; metadata read failed: ") + ex.what();
  }
  return result;
}

Result run_case(const Options& options) {
  if(options.implementation == "gauxc-cuest-api") {
#ifdef GAUXC_HAS_DEVICE
    if(options.exec_space == "DEVICE") {
      auto runtime = GauXC::DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9);
      return run_gauxc_cuest_api(options, runtime);
    }
#endif
    auto runtime = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
    return run_gauxc_cuest_api(options, runtime);
  }
  if(options.implementation == "cuest") {
#ifdef GAUXC_BENCHMARK_HAS_CUEST
    return run_cuest_native<void>(options);
#else
    return metadata_only_result(options, "unavailable", "libcuEST was not found when building cuest_comparison_case");
#endif
  }

#ifdef GAUXC_HAS_DEVICE
  if(options.exec_space == "DEVICE") {
    auto runtime = GauXC::DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9);
    return run_native(options, runtime);
  }
#endif
  auto runtime = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
  return run_native(options, runtime);
}

} // namespace

int main(int argc, char** argv) {
#ifdef GAUXC_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  int status = 0;
  try {
    const auto options = parse_options(argc, argv);
    const auto result = run_case(options);
#ifdef GAUXC_HAS_MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) emit_csv(options, result);
#else
    emit_csv(options, result);
#endif
  } catch(const std::exception& ex) {
    std::cerr << "cuest_comparison_case error: " << ex.what() << std::endl;
    status = 1;
  }

#ifdef GAUXC_HAS_MPI
  MPI_Finalize();
#endif
  return status;
}
