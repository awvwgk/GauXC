#pragma once

#include <gauxc/c/gauxc_config.h>

#ifndef GAUXC_HAS_C
#error "gauxc/c/cuest.h requires GAUXC_HAS_C"
#endif

#ifndef GAUXC_HAS_CUDA
#error "gauxc/c/cuest.h requires GAUXC_HAS_CUDA"
#endif

#ifndef GAUXC_HAS_CUEST_API
#error "gauxc/c/cuest.h requires GAUXC_HAS_CUEST_API"
#endif

#include <stddef.h>
#include <stdint.h>

#ifndef GAUXC_CUEST_API
#define GAUXC_CUEST_API
#endif

#define GAUXC_CUEST_VER_MAJOR 0
#define GAUXC_CUEST_VER_MINOR 1
#define GAUXC_CUEST_VER_PATCH 1
#define GAUXC_CUEST_VERSION (GAUXC_CUEST_VER_MAJOR * 10000 + GAUXC_CUEST_VER_MINOR * 100 + GAUXC_CUEST_VER_PATCH)

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  GAUXC_CUEST_STATUS_SUCCESS = 0,
  GAUXC_CUEST_STATUS_EXCEPTION = 1,
  GAUXC_CUEST_STATUS_NULL_POINTER = 2,
  GAUXC_CUEST_STATUS_INVALID_ARGUMENT = 3,
  GAUXC_CUEST_STATUS_INVALID_SIZE = 4,
  GAUXC_CUEST_STATUS_INVALID_TYPE = 5,
  GAUXC_CUEST_STATUS_INVALID_PARAMETER = 6,
  GAUXC_CUEST_STATUS_INVALID_ATTRIBUTE = 7,
  GAUXC_CUEST_STATUS_INVALID_HANDLE = 8,
  GAUXC_CUEST_STATUS_UNKNOWN_ERROR = 9,
  GAUXC_CUEST_STATUS_UNSUPPORTED_ARGUMENT = 10,
  GAUXC_CUEST_STATUS_UNSUPPORTED_ARCHITECTURE = 11
} gauxcStatus_t;

typedef void* gauxcHandle_t;
typedef void* gauxcAOShell_t;
typedef void* gauxcAOBasis_t;
typedef void* gauxcAtomGrid_t;
typedef void* gauxcMolecularGrid_t;
typedef void* gauxcXCIntPlan_t;
typedef void* gauxcDFIntPlan_t;

typedef void* gauxcHandleParameters_t;
typedef void* gauxcAOShellParameters_t;
typedef void* gauxcAOBasisParameters_t;
typedef void* gauxcAtomGridParameters_t;
typedef void* gauxcMolecularGridParameters_t;
typedef void* gauxcXCIntPlanParameters_t;
typedef void* gauxcXCPotentialRKSComputeParameters_t;
typedef void* gauxcXCPotentialUKSComputeParameters_t;
typedef void* gauxcNonlocalXCPotentialRKSComputeParameters_t;
typedef void* gauxcNonlocalXCPotentialUKSComputeParameters_t;
typedef void* gauxcDFSymmetricExchangeComputeParameters_t;

typedef struct {
  uintptr_t hostBuffer;
  size_t hostBufferSizeInBytes;
  uintptr_t deviceBuffer;
  size_t deviceBufferSizeInBytes;
} gauxcWorkspace_t;

typedef struct {
  size_t hostBufferSizeInBytes;
  size_t deviceBufferSizeInBytes;
} gauxcWorkspaceDescriptor_t;

typedef enum {
  GAUXC_CUEST_HANDLE_PARAMETERS = 0,
  GAUXC_CUEST_AOSHELL_PARAMETERS = 1,
  GAUXC_CUEST_AOBASIS_PARAMETERS = 2,
  GAUXC_CUEST_ATOMGRID_PARAMETERS = 6,
  GAUXC_CUEST_MOLECULARGRID_PARAMETERS = 7,
  GAUXC_CUEST_DFINTPLAN_PARAMETERS = 16,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS = 17,
  GAUXC_CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS = 18,
  GAUXC_CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS = 19,
  GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS = 21,
  GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS = 22,
  GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS = 23
} gauxcParametersType_t;

typedef enum {
  GAUXC_CUEST_DEFAULT_MATH_MODE = 0,
  GAUXC_CUEST_NATIVE_FP64_MATH_MODE = 1
} gauxcMathMode_t;

typedef enum {
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HF = 0,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP1 = 1,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP5 = 2,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97 = 3,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_BLYP = 4,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06L = 5,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE = 6,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE0 = 7,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_R2SCAN = 8,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_SVWN5 = 9,
  GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97MV = 10
} gauxcXCIntPlanParametersFunctional_t;

typedef enum {
  GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_SCALE = 0,
  GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_C = 1,
  GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B = 2
} gauxcNonlocalXCPotentialRKSComputeParametersAttributes_t;

typedef enum {
  GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_SCALE = 0,
  GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_C = 1,
  GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_B = 2
} gauxcNonlocalXCPotentialUKSComputeParametersAttributes_t;

typedef enum {
  GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT = 0,
  GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_MODULUS_COUNT = 1
} gauxcDFSymmetricExchangeComputeParametersAttributes_t;

GAUXC_CUEST_API gauxcStatus_t gauxcCreate(gauxcHandleParameters_t params, gauxcHandle_t* handle);
GAUXC_CUEST_API gauxcStatus_t gauxcDestroy(gauxcHandle_t handle);
GAUXC_CUEST_API gauxcStatus_t gauxcSetMathMode(gauxcHandle_t handle, gauxcMathMode_t mode);
GAUXC_CUEST_API gauxcStatus_t gauxcGetMathMode(gauxcHandle_t handle, gauxcMathMode_t* mode);
GAUXC_CUEST_API gauxcStatus_t gauxcSetComputeCapabilityTarget(gauxcHandle_t handle, uint32_t major, uint32_t minor);
GAUXC_CUEST_API gauxcStatus_t gauxcGetComputeCapabilityTarget(gauxcHandle_t handle, uint32_t* major, uint32_t* minor);
GAUXC_CUEST_API gauxcStatus_t gauxcGetMajorVersion(gauxcHandle_t handle, uint32_t* major);
GAUXC_CUEST_API gauxcStatus_t gauxcGetMinorVersion(gauxcHandle_t handle, uint32_t* minor);

GAUXC_CUEST_API gauxcStatus_t gauxcParametersCreate(gauxcParametersType_t type, void** outParameters);
GAUXC_CUEST_API gauxcStatus_t gauxcParametersDestroy(gauxcParametersType_t type, void* parameters);
GAUXC_CUEST_API gauxcStatus_t gauxcParametersQuery(gauxcParametersType_t type, const void* parameters, int attribute, void* value, size_t valueSize);
GAUXC_CUEST_API gauxcStatus_t gauxcParametersConfigure(gauxcParametersType_t type, void* parameters, int attribute, const void* value, size_t valueSize);

GAUXC_CUEST_API gauxcStatus_t gauxcAOShellCreate(gauxcHandle_t handle, int32_t isPure, uint64_t L, uint64_t numPrimitive, const double* exponents, const double* coefficients, const gauxcAOShellParameters_t parameters, gauxcAOShell_t* outShell);
GAUXC_CUEST_API gauxcStatus_t gauxcAOShellDestroy(gauxcAOShell_t shell);
GAUXC_CUEST_API gauxcStatus_t gauxcAOBasisCreateWorkspaceQuery(gauxcHandle_t handle, uint64_t numAtoms, const uint64_t* numShellsPerAtom, const gauxcAOShell_t* shells, const gauxcAOBasisParameters_t parameters, gauxcWorkspaceDescriptor_t* persistentWorkspaceDescriptor, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, gauxcAOBasis_t* outBasis);
GAUXC_CUEST_API gauxcStatus_t gauxcAOBasisCreate(gauxcHandle_t handle, uint64_t numAtoms, const uint64_t* numShellsPerAtom, const gauxcAOShell_t* shells, const gauxcAOBasisParameters_t parameters, gauxcWorkspace_t* persistentWorkspace, gauxcWorkspace_t* temporaryWorkspace, gauxcAOBasis_t* outBasis);
GAUXC_CUEST_API gauxcStatus_t gauxcAOBasisDestroy(gauxcAOBasis_t basis);

GAUXC_CUEST_API gauxcStatus_t gauxcAtomGridCreate(gauxcHandle_t handle, uint64_t numRadialPoints, const double* radialNodes, const double* radialWeights, const uint64_t* numAngularPoints, const gauxcAtomGridParameters_t parameters, gauxcAtomGrid_t* outAtomGrid);
GAUXC_CUEST_API gauxcStatus_t gauxcAtomGridDestroy(gauxcAtomGrid_t atomGrid);
GAUXC_CUEST_API gauxcStatus_t gauxcMolecularGridCreateWorkspaceQuery(gauxcHandle_t handle, uint64_t numAtoms, const gauxcAtomGrid_t* atomGrid, const double* xyz, const gauxcMolecularGridParameters_t parameters, gauxcWorkspaceDescriptor_t* persistentWorkspaceDescriptor, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, gauxcMolecularGrid_t* outGrid);
GAUXC_CUEST_API gauxcStatus_t gauxcMolecularGridCreate(gauxcHandle_t handle, uint64_t numAtoms, const gauxcAtomGrid_t* atomGrid, const double* xyz, const gauxcMolecularGridParameters_t parameters, gauxcWorkspace_t* persistentWorkspace, gauxcWorkspace_t* temporaryWorkspace, gauxcMolecularGrid_t* outGrid);
GAUXC_CUEST_API gauxcStatus_t gauxcMolecularGridDestroy(gauxcMolecularGrid_t grid);

GAUXC_CUEST_API gauxcStatus_t gauxcXCIntPlanCreateWorkspaceQuery(gauxcHandle_t handle, const gauxcAOBasis_t basis, const gauxcMolecularGrid_t grid, gauxcXCIntPlanParametersFunctional_t functional, const gauxcXCIntPlanParameters_t parameters, gauxcWorkspaceDescriptor_t* persistentWorkspaceDescriptor, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, gauxcXCIntPlan_t* outPlan);
GAUXC_CUEST_API gauxcStatus_t gauxcXCIntPlanCreate(gauxcHandle_t handle, const gauxcAOBasis_t basis, const gauxcMolecularGrid_t grid, gauxcXCIntPlanParametersFunctional_t functional, const gauxcXCIntPlanParameters_t parameters, gauxcWorkspace_t* persistentWorkspace, gauxcWorkspace_t* temporaryWorkspace, gauxcXCIntPlan_t* outPlan);
GAUXC_CUEST_API gauxcStatus_t gauxcXCIntPlanDestroy(gauxcXCIntPlan_t plan);

GAUXC_CUEST_API gauxcStatus_t gauxcXCPotentialRKSCompute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcXCPotentialRKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API gauxcStatus_t gauxcXCPotentialRKSComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcXCPotentialRKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API gauxcStatus_t gauxcNonlocalXCPotentialRKSCompute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialRKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API gauxcStatus_t gauxcNonlocalXCPotentialRKSComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialRKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API gauxcStatus_t gauxcNonlocalXCPotentialUKSCompute(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialUKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t* temporaryWorkspace, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API gauxcStatus_t gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcXCIntPlan_t plan, const gauxcNonlocalXCPotentialUKSComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix);

GAUXC_CUEST_API gauxcStatus_t gauxcDFSymmetricExchangeCompute(gauxcHandle_t handle, const gauxcDFIntPlan_t plan, const gauxcDFSymmetricExchangeComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix);
GAUXC_CUEST_API gauxcStatus_t gauxcDFSymmetricExchangeComputeWorkspaceQuery(gauxcHandle_t handle, const gauxcDFIntPlan_t plan, const gauxcDFSymmetricExchangeComputeParameters_t parameters, const gauxcWorkspaceDescriptor_t* variableBufferSize, gauxcWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix);

#ifdef GAUXC_HAS_CUEST_ABI
typedef gauxcStatus_t cuestStatus_t;
typedef gauxcHandle_t cuestHandle_t;
typedef gauxcHandleParameters_t cuestHandleParameters_t;
typedef gauxcMathMode_t cuestMathMode_t;
typedef gauxcXCIntPlan_t cuestXCIntPlan_t;
typedef gauxcDFIntPlan_t cuestDFIntPlan_t;
typedef gauxcWorkspace_t cuestWorkspace_t;
typedef gauxcWorkspaceDescriptor_t cuestWorkspaceDescriptor_t;
typedef gauxcDFSymmetricExchangeComputeParameters_t cuestDFSymmetricExchangeComputeParameters_t;
typedef gauxcNonlocalXCPotentialRKSComputeParameters_t cuestNonlocalXCPotentialRKSComputeParameters_t;
typedef gauxcNonlocalXCPotentialUKSComputeParameters_t cuestNonlocalXCPotentialUKSComputeParameters_t;
GAUXC_CUEST_API cuestStatus_t cuestCreate(cuestHandleParameters_t params, cuestHandle_t* handle);
GAUXC_CUEST_API cuestStatus_t cuestDestroy(cuestHandle_t handle);
GAUXC_CUEST_API cuestStatus_t cuestGetMajorVersion(cuestHandle_t handle, uint32_t* major);
GAUXC_CUEST_API cuestStatus_t cuestGetMinorVersion(cuestHandle_t handle, uint32_t* minor);
GAUXC_CUEST_API cuestStatus_t cuestNonlocalXCPotentialRKSCompute(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialRKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API cuestStatus_t cuestNonlocalXCPotentialRKSComputeWorkspaceQuery(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialRKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API cuestStatus_t cuestNonlocalXCPotentialUKSCompute(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialUKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspace_t* temporaryWorkspace, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API cuestStatus_t cuestNonlocalXCPotentialUKSComputeWorkspaceQuery(cuestHandle_t handle, const cuestXCIntPlan_t plan, const cuestNonlocalXCPotentialUKSComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupiedAlpha, uint64_t numOccupiedBeta, const double* coefficientMatrixAlpha, const double* coefficientMatrixBeta, double* outXCEnergy, double* outXCPotentialMatrix);
GAUXC_CUEST_API cuestStatus_t cuestDFSymmetricExchangeCompute(cuestHandle_t handle, const cuestDFIntPlan_t plan, const cuestDFSymmetricExchangeComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspace_t* temporaryWorkspace, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix);
GAUXC_CUEST_API cuestStatus_t cuestDFSymmetricExchangeComputeWorkspaceQuery(cuestHandle_t handle, const cuestDFIntPlan_t plan, const cuestDFSymmetricExchangeComputeParameters_t parameters, const cuestWorkspaceDescriptor_t* variableBufferSize, cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor, uint64_t numOccupied, const double* coefficientMatrix, double* outExchangeMatrix);
#endif

#ifndef GAUXC_CUEST_NO_ALIASES
typedef gauxcStatus_t cuestStatus_t;
typedef gauxcHandle_t cuestHandle_t;
typedef gauxcAOShell_t cuestAOShell_t;
typedef gauxcAOBasis_t cuestAOBasis_t;
typedef gauxcAtomGrid_t cuestAtomGrid_t;
typedef gauxcMolecularGrid_t cuestMolecularGrid_t;
typedef gauxcXCIntPlan_t cuestXCIntPlan_t;
typedef gauxcDFIntPlan_t cuestDFIntPlan_t;
typedef gauxcWorkspace_t cuestWorkspace_t;
typedef gauxcWorkspaceDescriptor_t cuestWorkspaceDescriptor_t;
typedef gauxcParametersType_t cuestParametersType_t;
typedef gauxcMathMode_t cuestMathMode_t;
typedef gauxcXCIntPlanParametersFunctional_t cuestXCIntPlanParametersFunctional_t;
typedef gauxcDFSymmetricExchangeComputeParameters_t cuestDFSymmetricExchangeComputeParameters_t;
typedef gauxcNonlocalXCPotentialRKSComputeParameters_t cuestNonlocalXCPotentialRKSComputeParameters_t;
typedef gauxcNonlocalXCPotentialUKSComputeParameters_t cuestNonlocalXCPotentialUKSComputeParameters_t;
#define CUEST_STATUS_SUCCESS GAUXC_CUEST_STATUS_SUCCESS
#define CUEST_DEFAULT_MATH_MODE GAUXC_CUEST_DEFAULT_MATH_MODE
#define CUEST_NATIVE_FP64_MATH_MODE GAUXC_CUEST_NATIVE_FP64_MATH_MODE
#define CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE GAUXC_CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE
#define CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_SCALE GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_SCALE
#define CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_C GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_C
#define CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B GAUXC_CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B
#define CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_SCALE GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_SCALE
#define CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_C GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_C
#define CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_B GAUXC_CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_B
#define CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT
#define CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_MODULUS_COUNT GAUXC_CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_MODULUS_COUNT
#define cuestCreate gauxcCreate
#define cuestDestroy gauxcDestroy
#define cuestSetMathMode gauxcSetMathMode
#define cuestGetMathMode gauxcGetMathMode
#define cuestGetMajorVersion gauxcGetMajorVersion
#define cuestGetMinorVersion gauxcGetMinorVersion
#define cuestParametersCreate gauxcParametersCreate
#define cuestParametersDestroy gauxcParametersDestroy
#define cuestAOShellCreate gauxcAOShellCreate
#define cuestAOShellDestroy gauxcAOShellDestroy
#define cuestAOBasisCreate gauxcAOBasisCreate
#define cuestAOBasisDestroy gauxcAOBasisDestroy
#define cuestXCIntPlanCreate gauxcXCIntPlanCreate
#define cuestXCIntPlanDestroy gauxcXCIntPlanDestroy
#define cuestNonlocalXCPotentialRKSCompute gauxcNonlocalXCPotentialRKSCompute
#define cuestNonlocalXCPotentialRKSComputeWorkspaceQuery gauxcNonlocalXCPotentialRKSComputeWorkspaceQuery
#define cuestNonlocalXCPotentialUKSCompute gauxcNonlocalXCPotentialUKSCompute
#define cuestNonlocalXCPotentialUKSComputeWorkspaceQuery gauxcNonlocalXCPotentialUKSComputeWorkspaceQuery
#define cuestDFSymmetricExchangeCompute gauxcDFSymmetricExchangeCompute
#define cuestDFSymmetricExchangeComputeWorkspaceQuery gauxcDFSymmetricExchangeComputeWorkspaceQuery
#endif

#ifdef __cplusplus
}
#endif
