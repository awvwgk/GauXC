#include <catch2/catch.hpp>

#include <gauxc/c/gauxc_config.h>

#ifdef GAUXC_HAS_CUEST_API

#include <gauxc/c/cuest.h>

TEST_CASE("cuEST-compatible API context", "[cuest-api]") {
  void* params = nullptr;
  REQUIRE(gauxcParametersCreate(GAUXC_CUEST_HANDLE_PARAMETERS, &params) == GAUXC_CUEST_STATUS_SUCCESS);

  gauxcHandle_t handle = nullptr;
  REQUIRE(gauxcCreate(params, &handle) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(handle != nullptr);

  uint32_t major = 99;
  uint32_t minor = 99;
  REQUIRE(gauxcGetMajorVersion(handle, &major) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcGetMinorVersion(handle, &minor) == GAUXC_CUEST_STATUS_SUCCESS);
  CHECK(major == GAUXC_CUEST_VER_MAJOR);
  CHECK(minor == GAUXC_CUEST_VER_MINOR);

  REQUIRE(gauxcSetMathMode(handle, GAUXC_CUEST_NATIVE_FP64_MATH_MODE) == GAUXC_CUEST_STATUS_SUCCESS);
  gauxcMathMode_t mode = GAUXC_CUEST_DEFAULT_MATH_MODE;
  REQUIRE(gauxcGetMathMode(handle, &mode) == GAUXC_CUEST_STATUS_SUCCESS);
  CHECK(mode == GAUXC_CUEST_NATIVE_FP64_MATH_MODE);

  REQUIRE(gauxcDestroy(handle) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcParametersDestroy(GAUXC_CUEST_HANDLE_PARAMETERS, params) == GAUXC_CUEST_STATUS_SUCCESS);
}

TEST_CASE("cuEST-compatible API stores NLC parameters", "[cuest-api]") {
  void* params = nullptr;
  REQUIRE(gauxcParametersCreate(GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS, &params) == GAUXC_CUEST_STATUS_SUCCESS);

  const double b = 6.3;
  REQUIRE(gauxcParametersConfigure(
    GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS,
    params,
    GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B,
    &b,
    sizeof(double)) == GAUXC_CUEST_STATUS_SUCCESS);

  double read_b = 0.0;
  REQUIRE(gauxcParametersQuery(
    GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS,
    params,
    GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B,
    &read_b,
    sizeof(double)) == GAUXC_CUEST_STATUS_SUCCESS);
  CHECK(read_b == b);

  REQUIRE(gauxcParametersDestroy(GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS, params) == GAUXC_CUEST_STATUS_SUCCESS);
}

TEST_CASE("cuEST-compatible API stores EXX parameters", "[cuest-api]") {
  void* params = nullptr;
  REQUIRE(gauxcParametersCreate(GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, &params) == GAUXC_CUEST_STATUS_SUCCESS);

  const uint64_t slices = 6;
  REQUIRE(gauxcParametersConfigure(
    GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS,
    params,
    GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT,
    &slices,
    sizeof(uint64_t)) == GAUXC_CUEST_STATUS_SUCCESS);

  uint64_t read_slices = 0;
  REQUIRE(gauxcParametersQuery(
    GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS,
    params,
    GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT,
    &read_slices,
    sizeof(uint64_t)) == GAUXC_CUEST_STATUS_SUCCESS);
  CHECK(read_slices == slices);

  REQUIRE(gauxcParametersDestroy(GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, params) == GAUXC_CUEST_STATUS_SUCCESS);
}

TEST_CASE("cuEST-compatible API owns shell basis grid plan handles", "[cuest-api]") {
  void* params = nullptr;
  REQUIRE(gauxcParametersCreate(GAUXC_CUEST_HANDLE_PARAMETERS, &params) == GAUXC_CUEST_STATUS_SUCCESS);
  gauxcHandle_t handle = nullptr;
  REQUIRE(gauxcCreate(params, &handle) == GAUXC_CUEST_STATUS_SUCCESS);

  const double exponents[] = {1.0};
  const double coefficients[] = {1.0};
  gauxcAOShell_t shell = nullptr;
  REQUIRE(gauxcAOShellCreate(handle, 0, 0, 1, exponents, coefficients, nullptr, &shell) == GAUXC_CUEST_STATUS_SUCCESS);

  const uint64_t shells_per_atom[] = {1};
  gauxcWorkspace_t workspace{};
  gauxcAOBasis_t basis = nullptr;
  REQUIRE(gauxcAOBasisCreate(handle, 1, shells_per_atom, &shell, nullptr, &workspace, &workspace, &basis) == GAUXC_CUEST_STATUS_SUCCESS);

  const double radial_nodes[] = {0.5};
  const double radial_weights[] = {1.0};
  const uint64_t angular_points[] = {6};
  gauxcAtomGrid_t atom_grid = nullptr;
  REQUIRE(gauxcAtomGridCreate(handle, 1, radial_nodes, radial_weights, angular_points, nullptr, &atom_grid) == GAUXC_CUEST_STATUS_SUCCESS);

  const double xyz[] = {0.0, 0.0, 0.0};
  gauxcMolecularGrid_t molecular_grid = nullptr;
  REQUIRE(gauxcMolecularGridCreate(handle, 1, &atom_grid, xyz, nullptr, &workspace, &workspace, &molecular_grid) == GAUXC_CUEST_STATUS_SUCCESS);

  gauxcXCIntPlan_t plan = nullptr;
  REQUIRE(gauxcXCIntPlanCreate(handle, basis, molecular_grid, GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE, nullptr, &workspace, &workspace, &plan) == GAUXC_CUEST_STATUS_SUCCESS);

  gauxcWorkspaceDescriptor_t variable_buffer{};
  gauxcWorkspaceDescriptor_t temporary{};
  REQUIRE(gauxcXCPotentialRKSComputeWorkspaceQuery(handle, plan, nullptr, &variable_buffer, &temporary, 1, nullptr, nullptr, nullptr) == GAUXC_CUEST_STATUS_SUCCESS);
  CHECK(temporary.hostBufferSizeInBytes == 0);
  CHECK(temporary.deviceBufferSizeInBytes == 0);

  REQUIRE(gauxcNonlocalXCPotentialRKSComputeWorkspaceQuery(handle, plan, nullptr, &variable_buffer, &temporary, 1, nullptr, nullptr, nullptr) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery(handle, plan, nullptr, &variable_buffer, &temporary, 1, 1, nullptr, nullptr, nullptr, nullptr) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcDFSymmetricExchangeComputeWorkspaceQuery(handle, plan, nullptr, &variable_buffer, &temporary, 1, nullptr, nullptr) == GAUXC_CUEST_STATUS_SUCCESS);
  CHECK(gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery(handle, plan, nullptr, &variable_buffer, &temporary, 0, 1, nullptr, nullptr, nullptr, nullptr) == GAUXC_CUEST_STATUS_INVALID_ARGUMENT);
  CHECK(gauxcDFSymmetricExchangeComputeWorkspaceQuery(handle, plan, nullptr, &variable_buffer, &temporary, 0, nullptr, nullptr) == GAUXC_CUEST_STATUS_INVALID_ARGUMENT);

  REQUIRE(gauxcXCIntPlanDestroy(plan) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcMolecularGridDestroy(molecular_grid) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcAtomGridDestroy(atom_grid) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcAOBasisDestroy(basis) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcAOShellDestroy(shell) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcDestroy(handle) == GAUXC_CUEST_STATUS_SUCCESS);
  REQUIRE(gauxcParametersDestroy(GAUXC_CUEST_HANDLE_PARAMETERS, params) == GAUXC_CUEST_STATUS_SUCCESS);
}

#endif
