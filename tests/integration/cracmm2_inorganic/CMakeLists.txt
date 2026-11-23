add_executable(test_cracmm2_inorganic test_cracmm2_inorganic.cpp test_state_binding.cpp)
target_link_libraries(test_cracmm2_inorganic PRIVATE MKPP::cracmm2_inorganic GTest::gtest_main)
add_test(NAME CRACMM2Inorganic_host COMMAND test_cracmm2_inorganic)
set_tests_properties(CRACMM2Inorganic_host PROPERTIES LABELS cracmm2_inorganic)
add_custom_target(cracmm2_inorganic_checks DEPENDS test_cracmm2_inorganic)

add_executable(test_cracmm2_fine_modes test_fine_modes.cpp)
target_link_libraries(test_cracmm2_fine_modes PRIVATE MKPP::cracmm2_inorganic GTest::gtest_main)
add_dependencies(cracmm2_inorganic_checks test_cracmm2_fine_modes)
add_test(NAME CRACMM2Inorganic_fine_modes COMMAND test_cracmm2_fine_modes)
set_tests_properties(CRACMM2Inorganic_fine_modes PROPERTIES LABELS cracmm2_inorganic)

option(MKPP_TEST_LITE_FORTRAN "Compare the CRACMM2 host adapter with the original Lite Fortran" OFF)
if(MKPP_TEST_LITE_FORTRAN)
  enable_language(Fortran)
  set(reference_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../src/isorropialite/isolite1_0_src")
  add_executable(
    test_cracmm2_lite_reference test_reference.cpp lite_reference.f90 "${reference_dir}/isocom.f"
                                "${reference_dir}/isofwd.f" "${reference_dir}/isorev.f")
  target_include_directories(test_cracmm2_lite_reference PRIVATE "${reference_dir}")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_cracmm2_lite_reference
                           PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>)
  endif()
  target_link_libraries(test_cracmm2_lite_reference PRIVATE MKPP::cracmm2_inorganic)
  add_dependencies(cracmm2_inorganic_checks test_cracmm2_lite_reference)
  foreach(case_name ammonium-nitrate mixed-ions cold-inorganic sulfate-rich)
    add_test(NAME CRACMM2Inorganic_reference_${case_name} COMMAND test_cracmm2_lite_reference
                                                                  ${case_name})
    set_tests_properties(CRACMM2Inorganic_reference_${case_name} PROPERTIES LABELS
                                                                            cracmm2_inorganic)
  endforeach()
endif()

include("${CMAKE_CURRENT_LIST_DIR}/generated_host.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmaq_ii.cmake")
