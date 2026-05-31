#define GAUXC_CUEST_NO_ALIASES
#include <gauxc/c/cuest.h>

#ifndef GAUXC_HAS_CUEST_ABI
#error "GAUXC_HAS_CUEST_ABI is required for this compile test"
#endif

int main() {
  void* params = reinterpret_cast<void*>(1);
  cuestHandle_t handle = nullptr;
  auto status = cuestCreate(params, &handle);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return 1;
  uint32_t major = 0;
  status = cuestGetMajorVersion(handle, &major);
  if(status != GAUXC_CUEST_STATUS_SUCCESS) return 2;
  cuestWorkspaceDescriptor_t variable_buffer{};
  cuestWorkspaceDescriptor_t temporary{};
  status = cuestDFSymmetricExchangeComputeWorkspaceQuery(handle, nullptr, nullptr, &variable_buffer, &temporary, 1, nullptr, nullptr);
  if(status != GAUXC_CUEST_STATUS_INVALID_HANDLE) return 3;
  status = cuestNonlocalXCPotentialUKSComputeWorkspaceQuery(handle, nullptr, nullptr, &variable_buffer, &temporary, 1, 1, nullptr, nullptr, nullptr, nullptr);
  if(status != GAUXC_CUEST_STATUS_INVALID_HANDLE) return 4;
  return cuestDestroy(handle) == GAUXC_CUEST_STATUS_SUCCESS ? 0 : 5;
}
