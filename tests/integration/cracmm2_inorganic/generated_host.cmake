# Supply a CRACMM2 bundle emitted by mkpp.cli compile --host-api. Build its wrapper and its kernel
# sources together so the ABI and species ordering match.
set(MKPP_CRACMM2_HOST_API_DIR
    ""
    CACHE PATH "Generated CRACMM2 --host-api output directory")
if(MKPP_CRACMM2_HOST_API_DIR)
  foreach(required IN ITEMS mkpp.h mkpp.hpp mkpp_c_api.cpp cracmm2.hpp
                            cracmm2/host_declarations.hpp)
    if(NOT EXISTS "${MKPP_CRACMM2_HOST_API_DIR}/${required}")
      message(FATAL_ERROR "Missing ${required} in MKPP_CRACMM2_HOST_API_DIR")
    endif()
  endforeach()
  include("${CMAKE_CURRENT_LIST_DIR}/../../../src/mkpp_host/cmake/MKPPCompiledMechanism.cmake")
  mkpp_add_compiled_mechanism(cracmm2_host_abi_kernels "${MKPP_CRACMM2_HOST_API_DIR}" cracmm2)
  add_executable(test_cracmm2_generated_host "${CMAKE_CURRENT_LIST_DIR}/test_generated_host.cpp"
                                             "${MKPP_CRACMM2_HOST_API_DIR}/mkpp_c_api.cpp")
  target_include_directories(test_cracmm2_generated_host BEFORE
                             PRIVATE "${MKPP_CRACMM2_HOST_API_DIR}")
  target_link_libraries(test_cracmm2_generated_host PRIVATE MKPP::cracmm2_inorganic MKPP::mkpp_host
                                                            cracmm2_host_abi_kernels GTest::gtest)
  add_dependencies(cracmm2_inorganic_checks test_cracmm2_generated_host)
  add_test(NAME cracmm2_generated_host_thermodynamics COMMAND test_cracmm2_generated_host)
  set_tests_properties(cracmm2_generated_host_thermodynamics PROPERTIES LABELS cracmm2_inorganic)
endif()
