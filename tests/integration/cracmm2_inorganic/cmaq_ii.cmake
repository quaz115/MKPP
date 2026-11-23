# Optional external reference from USEPA/CMAQ 5.5.0.4, revision
# 1b70e0a395fc42cc34a4355807d39df1d3dbae84. No CMAQ source is redistributed.
set(MKPP_CMAQ_ISORROPIA_II_DIR
    ""
    CACHE PATH "Pinned CMAQ 5.5.0.4 CCTM/src/aero/aero6 directory")
if(MKPP_CMAQ_ISORROPIA_II_DIR)
  enable_language(Fortran)

  # Accept a sparse checkout or an exported source directory, with LF or CRLF.
  function(mkpp_check_cmaq_ii_source name expected)
    file(READ "${MKPP_CMAQ_ISORROPIA_II_DIR}/${name}" content)
    string(REPLACE "\r\n" "\n" content "${content}")
    string(SHA256 digest "${content}")
    if(NOT "${digest}" STREQUAL "${expected}")
      message(FATAL_ERROR "CMAQ ISORROPIA II reference ${name} differs from revision 1b70e0a3")
    endif()
  endfunction()
  mkpp_check_cmaq_ii_source(isocom.f
                            4a7cabcac3e02ba14e426758d31dc7f54c20520f28c001d89d65d5c5fe70b054)
  mkpp_check_cmaq_ii_source(isofwd.f
                            67b49a9bba1a4740b3f95e6cd0b2c77fb6da8c7f3ffa83201c1427f7911fdb0f)
  mkpp_check_cmaq_ii_source(isorev.f
                            b0c966abbff7c3f5003c22281008490542c57d3932fa43cd42f3023878cc24ee)
  mkpp_check_cmaq_ii_source(isrpia.inc
                            940017c28fd58074366e480b48a63c60e0fe87836ba359bd8456ffb42d17ce5d)

  # Keep ISOROPIA and its COMMON blocks out of the separate original-Lite test.
  add_executable(
    test_cracmm2_cmaq_ii
    "${CMAKE_CURRENT_LIST_DIR}/test_cmaq_ii.cpp" "${CMAKE_CURRENT_LIST_DIR}/cmaq_ii_reference.f90"
    "${MKPP_CMAQ_ISORROPIA_II_DIR}/isocom.f" "${MKPP_CMAQ_ISORROPIA_II_DIR}/isofwd.f"
    "${MKPP_CMAQ_ISORROPIA_II_DIR}/isorev.f")
  target_include_directories(test_cracmm2_cmaq_ii PRIVATE "${MKPP_CMAQ_ISORROPIA_II_DIR}")
  target_compile_options(test_cracmm2_cmaq_ii
                         PRIVATE $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=legacy>)
  target_link_libraries(test_cracmm2_cmaq_ii PRIVATE MKPP::cracmm2_inorganic)
  add_dependencies(cracmm2_inorganic_checks test_cracmm2_cmaq_ii)
  foreach(climate IN ITEMS cold warm)
    add_test(NAME cracmm2_cmaq_ii_${climate} COMMAND test_cracmm2_cmaq_ii ${climate})
    set_tests_properties(cracmm2_cmaq_ii_${climate} PROPERTIES LABELS cracmm2_inorganic)
  endforeach()
endif()
