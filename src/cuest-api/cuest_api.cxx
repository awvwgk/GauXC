#define GAUXC_CUEST_NO_ALIASES
#include <gauxc/c/cuest.h>

#include <cstring>
#include <exception>
#include <vector>

namespace {

enum class Kind : uint32_t { Context = 1, Parameters, Shell, Basis, AtomGrid, MolecularGrid, XCPlan };

struct Header { Kind kind; };
struct Context { Header header{Kind::Context}; gauxcMathMode_t math_mode = GAUXC_CUEST_DEFAULT_MATH_MODE; uint32_t cc_major = 0; uint32_t cc_minor = 0; };
struct Parameters {
  Header header{Kind::Parameters};
  gauxcParametersType_t type;
  double vv10_scale = 0.0;
  double vv10_c = 0.0;
  double vv10_b = 0.0;
  uint64_t int8_slice_count = 5;
  uint64_t int8_modulus_count = 8;
};
struct Shell { Header header{Kind::Shell}; int32_t is_pure; uint64_t l; std::vector<double> exponents; std::vector<double> coefficients; };
struct Basis { Header header{Kind::Basis}; uint64_t num_atoms; std::vector<uint64_t> shells_per_atom; std::vector<Shell> shells; uint64_t nao = 0; };
struct AtomGrid { Header header{Kind::AtomGrid}; std::vector<double> radial_nodes; std::vector<double> radial_weights; std::vector<uint64_t> angular_points; };
struct MolecularGrid { Header header{Kind::MolecularGrid}; uint64_t num_atoms; std::vector<double> xyz; };
struct XCPlan { Header header{Kind::XCPlan}; gauxcXCIntPlanParametersFunctional_t functional; uint64_t nao; };

template <typename T>
T* checked(void* ptr, Kind kind) {
  if(ptr == nullptr) return nullptr;
  auto* header = static_cast<Header*>(ptr);
  return header->kind == kind ? static_cast<T*>(ptr) : nullptr;
}

uint64_t shell_nao(const Shell& shell) {
  if(shell.is_pure && shell.l > 1) return 2 * shell.l + 1;
  return (shell.l + 1) * (shell.l + 2) / 2;
}

void zero_workspace(gauxcWorkspaceDescriptor_t* descriptor) {
  descriptor->hostBufferSizeInBytes = 0;
  descriptor->deviceBufferSizeInBytes = 0;
}

gauxcStatus_t validate_compute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(checked<XCPlan>(plan, Kind::XCPlan) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(variableBufferSize == nullptr || temporaryWorkspaceDescriptor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numOccupied == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  zero_workspace(temporaryWorkspaceDescriptor);
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t validate_uks_compute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(checked<XCPlan>(plan, Kind::XCPlan) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(variableBufferSize == nullptr || temporaryWorkspaceDescriptor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numOccupiedAlpha == 0 || numOccupiedBeta == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  zero_workspace(temporaryWorkspaceDescriptor);
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t validate_exchange(gauxcHandle_t handle, const gauxcDFIntPlan_t plan, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(checked<XCPlan>(plan, Kind::XCPlan) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(variableBufferSize == nullptr || temporaryWorkspaceDescriptor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numOccupied == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  zero_workspace(temporaryWorkspaceDescriptor);
  return GAUXC_CUEST_STATUS_SUCCESS;
}

} // namespace

extern "C" {

gauxcStatus_t gauxcCreate(gauxcHandleParameters_t params, gauxcHandle_t* handle) {
  if(params == nullptr || handle == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  try { *handle = new Context{}; return GAUXC_CUEST_STATUS_SUCCESS; }
  catch(const std::exception&) { return GAUXC_CUEST_STATUS_EXCEPTION; }
  catch(...) { return GAUXC_CUEST_STATUS_UNKNOWN_ERROR; }
}

gauxcStatus_t gauxcDestroy(gauxcHandle_t handle) {
  auto* context = checked<Context>(handle, Kind::Context);
  if(context == nullptr) return handle == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete context;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcSetMathMode(gauxcHandle_t handle, gauxcMathMode_t mode) {
  auto* context = checked<Context>(handle, Kind::Context);
  if(context == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(mode != GAUXC_CUEST_DEFAULT_MATH_MODE && mode != GAUXC_CUEST_NATIVE_FP64_MATH_MODE) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  context->math_mode = mode;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcGetMathMode(gauxcHandle_t handle, gauxcMathMode_t* mode) {
  auto* context = checked<Context>(handle, Kind::Context);
  if(context == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(mode == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  *mode = context->math_mode;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcSetComputeCapabilityTarget(gauxcHandle_t handle, uint32_t major, uint32_t minor) {
  auto* context = checked<Context>(handle, Kind::Context);
  if(context == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  context->cc_major = major;
  context->cc_minor = minor;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcGetComputeCapabilityTarget(gauxcHandle_t handle, uint32_t* major, uint32_t* minor) {
  auto* context = checked<Context>(handle, Kind::Context);
  if(context == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(major == nullptr || minor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  *major = context->cc_major;
  *minor = context->cc_minor;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcGetMajorVersion(gauxcHandle_t handle, uint32_t* major) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(major == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  *major = GAUXC_CUEST_VER_MAJOR;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcGetMinorVersion(gauxcHandle_t handle, uint32_t* minor) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(minor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  *minor = GAUXC_CUEST_VER_MINOR;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcParametersCreate(gauxcParametersType_t type, void** outParameters) {
  if(outParameters == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  *outParameters = new Parameters{{Kind::Parameters}, type};
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcParametersDestroy(gauxcParametersType_t, void* parameters) {
  auto* p = checked<Parameters>(parameters, Kind::Parameters);
  if(p == nullptr) return parameters == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete p;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcParametersQuery(gauxcParametersType_t type, const void* parameters, int attribute, void* value, size_t valueSize) {
  auto* p = checked<Parameters>(const_cast<void*>(parameters), Kind::Parameters);
  if(p == nullptr || value == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(p->type != type) return GAUXC_CUEST_STATUS_INVALID_PARAMETER;
  if(type == GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS || type == GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS) {
    if(valueSize != sizeof(double)) return GAUXC_CUEST_STATUS_INVALID_SIZE;
    double result = 0.0;
    if(attribute == 0) result = p->vv10_scale;
    else if(attribute == 1) result = p->vv10_c;
    else if(attribute == 2) result = p->vv10_b;
    else return GAUXC_CUEST_STATUS_INVALID_ATTRIBUTE;
    std::memcpy(value, &result, sizeof(double));
  } else if(type == GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS) {
    if(valueSize != sizeof(uint64_t)) return GAUXC_CUEST_STATUS_INVALID_SIZE;
    uint64_t result = 0;
    if(attribute == 0) result = p->int8_slice_count;
    else if(attribute == 1) result = p->int8_modulus_count;
    else return GAUXC_CUEST_STATUS_INVALID_ATTRIBUTE;
    std::memcpy(value, &result, sizeof(uint64_t));
  }
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcParametersConfigure(gauxcParametersType_t type, void* parameters, int attribute, const void* value, size_t valueSize) {
  auto* p = checked<Parameters>(parameters, Kind::Parameters);
  if(p == nullptr || value == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(p->type != type) return GAUXC_CUEST_STATUS_INVALID_PARAMETER;
  if(type == GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS || type == GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS) {
    if(valueSize != sizeof(double)) return GAUXC_CUEST_STATUS_INVALID_SIZE;
    double input = 0.0;
    std::memcpy(&input, value, sizeof(double));
    if(attribute == 0) p->vv10_scale = input;
    else if(attribute == 1) p->vv10_c = input;
    else if(attribute == 2) p->vv10_b = input;
    else return GAUXC_CUEST_STATUS_INVALID_ATTRIBUTE;
  } else if(type == GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS) {
    if(valueSize != sizeof(uint64_t)) return GAUXC_CUEST_STATUS_INVALID_SIZE;
    uint64_t input = 0;
    std::memcpy(&input, value, sizeof(uint64_t));
    if(attribute == 0) p->int8_slice_count = input;
    else if(attribute == 1) p->int8_modulus_count = input;
    else return GAUXC_CUEST_STATUS_INVALID_ATTRIBUTE;
  }
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAOShellCreate(gauxcHandle_t handle, int32_t isPure, uint64_t L, uint64_t numPrimitive, const double* exponents, const double* coefficients, const gauxcAOShellParameters_t, gauxcAOShell_t* outShell) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(exponents == nullptr || coefficients == nullptr || outShell == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numPrimitive == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  auto* shell = new Shell{{Kind::Shell}, isPure, L, {exponents, exponents + numPrimitive}, {coefficients, coefficients + numPrimitive}};
  *outShell = shell;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAOShellDestroy(gauxcAOShell_t shell) {
  auto* ptr = checked<Shell>(shell, Kind::Shell);
  if(ptr == nullptr) return shell == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete ptr;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAOBasisCreateWorkspaceQuery(gauxcHandle_t handle, uint64_t numAtoms, const uint64_t* numShellsPerAtom, const gauxcAOShell_t* shells, const gauxcAOBasisParameters_t, gauxcWorkspaceDescriptor_t* persistentWorkspaceDescriptor, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, gauxcAOBasis_t*) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(numShellsPerAtom == nullptr || shells == nullptr || persistentWorkspaceDescriptor == nullptr || temporaryWorkspaceDescriptor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numAtoms == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  zero_workspace(persistentWorkspaceDescriptor);
  zero_workspace(temporaryWorkspaceDescriptor);
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAOBasisCreate(gauxcHandle_t handle, uint64_t numAtoms, const uint64_t* numShellsPerAtom, const gauxcAOShell_t* shells, const gauxcAOBasisParameters_t parameters, gauxcWorkspace_t*, gauxcWorkspace_t*, gauxcAOBasis_t* outBasis) {
  gauxcWorkspaceDescriptor_t persistent{}, temporary{};
  auto status = gauxcAOBasisCreateWorkspaceQuery(handle, numAtoms, numShellsPerAtom, shells, parameters, &persistent, &temporary, nullptr);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(outBasis == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  auto* basis = new Basis{{Kind::Basis}, numAtoms, {}, {}, 0};
  basis->shells_per_atom.assign(numShellsPerAtom, numShellsPerAtom + numAtoms);
  uint64_t shell_count = 0;
  for(uint64_t i = 0; i < numAtoms; ++i) shell_count += numShellsPerAtom[i];
  for(uint64_t i = 0; i < shell_count; ++i) {
    auto* shell = checked<Shell>(shells[i], Kind::Shell);
    if(shell == nullptr) { delete basis; return GAUXC_CUEST_STATUS_INVALID_TYPE; }
    basis->shells.push_back(*shell);
    basis->nao += shell_nao(*shell);
  }
  *outBasis = basis;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAOBasisDestroy(gauxcAOBasis_t basis) {
  auto* ptr = checked<Basis>(basis, Kind::Basis);
  if(ptr == nullptr) return basis == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete ptr;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAtomGridCreate(gauxcHandle_t handle, uint64_t numRadialPoints, const double* radialNodes, const double* radialWeights, const uint64_t* numAngularPoints, const gauxcAtomGridParameters_t, gauxcAtomGrid_t* outAtomGrid) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(radialNodes == nullptr || radialWeights == nullptr || numAngularPoints == nullptr || outAtomGrid == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numRadialPoints == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  auto* grid = new AtomGrid{{Kind::AtomGrid}, {radialNodes, radialNodes + numRadialPoints}, {radialWeights, radialWeights + numRadialPoints}, {numAngularPoints, numAngularPoints + numRadialPoints}};
  *outAtomGrid = grid;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcAtomGridDestroy(gauxcAtomGrid_t atomGrid) {
  auto* ptr = checked<AtomGrid>(atomGrid, Kind::AtomGrid);
  if(ptr == nullptr) return atomGrid == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete ptr;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcMolecularGridCreateWorkspaceQuery(gauxcHandle_t handle, uint64_t numAtoms, const gauxcAtomGrid_t* atomGrid, const double* xyz, const gauxcMolecularGridParameters_t, gauxcWorkspaceDescriptor_t* persistentWorkspaceDescriptor, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, gauxcMolecularGrid_t*) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(atomGrid == nullptr || xyz == nullptr || persistentWorkspaceDescriptor == nullptr || temporaryWorkspaceDescriptor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  if(numAtoms == 0) return GAUXC_CUEST_STATUS_INVALID_ARGUMENT;
  zero_workspace(persistentWorkspaceDescriptor);
  zero_workspace(temporaryWorkspaceDescriptor);
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcMolecularGridCreate(gauxcHandle_t handle, uint64_t numAtoms, const gauxcAtomGrid_t* atomGrid, const double* xyz, const gauxcMolecularGridParameters_t parameters, gauxcWorkspace_t*, gauxcWorkspace_t*, gauxcMolecularGrid_t* outGrid) {
  gauxcWorkspaceDescriptor_t persistent{}, temporary{};
  auto status = gauxcMolecularGridCreateWorkspaceQuery(handle, numAtoms, atomGrid, xyz, parameters, &persistent, &temporary, nullptr);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(outGrid == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  for(uint64_t i = 0; i < numAtoms; ++i) if(checked<AtomGrid>(atomGrid[i], Kind::AtomGrid) == nullptr) return GAUXC_CUEST_STATUS_INVALID_TYPE;
  auto* grid = new MolecularGrid{{Kind::MolecularGrid}, numAtoms, {xyz, xyz + 3 * numAtoms}};
  *outGrid = grid;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcMolecularGridDestroy(gauxcMolecularGrid_t grid) {
  auto* ptr = checked<MolecularGrid>(grid, Kind::MolecularGrid);
  if(ptr == nullptr) return grid == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete ptr;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcXCIntPlanCreateWorkspaceQuery(gauxcHandle_t handle, const gauxcAOBasis_t basis, const gauxcMolecularGrid_t grid, gauxcXCIntPlanParametersFunctional_t, const gauxcXCIntPlanParameters_t, gauxcWorkspaceDescriptor_t* persistentWorkspaceDescriptor, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, gauxcXCIntPlan_t*) {
  if(checked<Context>(handle, Kind::Context) == nullptr) return GAUXC_CUEST_STATUS_INVALID_HANDLE;
  if(checked<Basis>(basis, Kind::Basis) == nullptr || checked<MolecularGrid>(grid, Kind::MolecularGrid) == nullptr) return GAUXC_CUEST_STATUS_INVALID_TYPE;
  if(persistentWorkspaceDescriptor == nullptr || temporaryWorkspaceDescriptor == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  zero_workspace(persistentWorkspaceDescriptor);
  zero_workspace(temporaryWorkspaceDescriptor);
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcXCIntPlanCreate(gauxcHandle_t handle, const gauxcAOBasis_t basis, const gauxcMolecularGrid_t grid, gauxcXCIntPlanParametersFunctional_t functional, const gauxcXCIntPlanParameters_t parameters, gauxcWorkspace_t*, gauxcWorkspace_t*, gauxcXCIntPlan_t* outPlan) {
  gauxcWorkspaceDescriptor_t persistent{}, temporary{};
  auto status = gauxcXCIntPlanCreateWorkspaceQuery(handle, basis, grid, functional, parameters, &persistent, &temporary, nullptr);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(outPlan == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  auto* basis_ptr = checked<Basis>(basis, Kind::Basis);
  *outPlan = new XCPlan{{Kind::XCPlan}, functional, basis_ptr->nao};
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcXCIntPlanDestroy(gauxcXCIntPlan_t plan) {
  auto* ptr = checked<XCPlan>(plan, Kind::XCPlan);
  if(ptr == nullptr) return plan == nullptr ? GAUXC_CUEST_STATUS_NULL_POINTER : GAUXC_CUEST_STATUS_INVALID_TYPE;
  delete ptr;
  return GAUXC_CUEST_STATUS_SUCCESS;
}

gauxcStatus_t gauxcXCPotentialRKSComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcXCPotentialRKSComputeParameters_t, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double*, double*, double*) {
  return validate_compute(handle, plan, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied);
}

gauxcStatus_t gauxcXCPotentialRKSCompute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcXCPotentialRKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t*, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix) {
  gauxcWorkspaceDescriptor_t temporary{};
  auto status = gauxcXCPotentialRKSComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, &temporary, numOccupied, coefficientMatrix, outXCEnergy, outXCPotentialMatrix);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(coefficientMatrix == nullptr || outXCEnergy == nullptr || outXCPotentialMatrix == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  return GAUXC_CUEST_STATUS_UNSUPPORTED_ARGUMENT;
}

gauxcStatus_t gauxcNonlocalXCPotentialRKSComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialRKSComputeParameters_t, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double*, double*, double*) {
  return validate_compute(handle, plan, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied);
}

gauxcStatus_t gauxcNonlocalXCPotentialRKSCompute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialRKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t*, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix) {
  gauxcWorkspaceDescriptor_t temporary{};
  auto status = gauxcNonlocalXCPotentialRKSComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, &temporary, numOccupied, coefficientMatrix, outXCEnergy, outXCPotentialMatrix);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(coefficientMatrix == nullptr || outXCEnergy == nullptr || outXCPotentialMatrix == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  return GAUXC_CUEST_STATUS_UNSUPPORTED_ARGUMENT;
}

gauxcStatus_t gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialUKSComputeParameters_t, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double*, const double*, double*, double*) {
  return validate_uks_compute(handle, plan, variableBufferSize, temporaryWorkspaceDescriptor, numOccupiedAlpha, numOccupiedBeta);
}

gauxcStatus_t gauxcNonlocalXCPotentialUKSCompute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialUKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t*, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix) {
  gauxcWorkspaceDescriptor_t temporary{};
  auto status = gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, &temporary, numOccupiedAlpha, numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrix);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(coefficientMatrixAlpha == nullptr || coefficientMatrixBeta == nullptr || outXCEnergy == nullptr || outXCPotentialMatrix == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  return GAUXC_CUEST_STATUS_UNSUPPORTED_ARGUMENT;
}

gauxcStatus_t gauxcDFSymmetricExchangeComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcDFIntPlan_t plan, const gauxcDFSymmetricExchangeComputeParameters_t, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double*, double*) {
  return validate_exchange(handle, plan, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied);
}

gauxcStatus_t gauxcDFSymmetricExchangeCompute(gauxcHandle_t handle, const gauxcDFIntPlan_t plan, const gauxcDFSymmetricExchangeComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t*, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix) {
  gauxcWorkspaceDescriptor_t temporary{};
  auto status = gauxcDFSymmetricExchangeComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, &temporary, numOccupied, coefficientMatrix, outExchangeMatrix);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return status;
  if(coefficientMatrix == nullptr || outExchangeMatrix == nullptr) return GAUXC_CUEST_STATUS_NULL_POINTER;
  return GAUXC_CUEST_STATUS_UNSUPPORTED_ARGUMENT;
}

#ifdef GAUXC_HAS_CUEST_ABI
cuestStatus_t cuestCreate(cuestHandleParameters_t params, cuestHandle_t* handle) { return gauxcCreate(params, handle); }
cuestStatus_t cuestDestroy(cuestHandle_t handle) { return gauxcDestroy(handle); }
cuestStatus_t cuestGetMajorVersion(cuestHandle_t handle, uint32_t* major) { return gauxcGetMajorVersion(handle, major); }
cuestStatus_t cuestGetMinorVersion(cuestHandle_t handle, uint32_t* minor) { return gauxcGetMinorVersion(handle, minor); }
cuestStatus_t cuestNonlocalXCPotentialRKSCompute(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialRKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix) {
  return gauxcNonlocalXCPotentialRKSCompute(handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, outXCEnergy, outXCPotentialMatrix);
}
cuestStatus_t cuestNonlocalXCPotentialRKSComputeWorkspaceQuery(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialRKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix) {
  return gauxcNonlocalXCPotentialRKSComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, coefficientMatrix, outXCEnergy, outXCPotentialMatrix);
}
cuestStatus_t cuestNonlocalXCPotentialUKSCompute(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialUKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspace_t* temporaryWorkspace, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix) {
  return gauxcNonlocalXCPotentialUKSCompute(handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupiedAlpha, numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrix);
}
cuestStatus_t cuestNonlocalXCPotentialUKSComputeWorkspaceQuery(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialUKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix) {
  return gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupiedAlpha, numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrix);
}
cuestStatus_t cuestDFSymmetricExchangeCompute(cuestHandle_t handle, const cuestDFIntPlan_t plan, const cuestDFSymmetricExchangeComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix) {
  return gauxcDFSymmetricExchangeCompute(handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, outExchangeMatrix);
}
cuestStatus_t cuestDFSymmetricExchangeComputeWorkspaceQuery(cuestHandle_t handle, const cuestDFIntPlan_t plan, const cuestDFSymmetricExchangeComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix) {
  return gauxcDFSymmetricExchangeComputeWorkspaceQuery(handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, coefficientMatrix, outExchangeMatrix);
}
#endif

} // extern "C"
