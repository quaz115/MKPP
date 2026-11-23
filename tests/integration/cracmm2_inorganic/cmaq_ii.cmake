# Optional external references. Each executable links one ISOROPIA and its COMMON blocks.
set(MKPP_CMAQ_ISORROPIA_II_DIR
    ""
    CACHE PATH "Pinned CMAQ 5.5.0.4 CCTM/src/aero/aero6 directory")
set(MKPP_EPFL_ISORROPIA_II_DIR
    ""
    CACHE PATH "Pinned EPFL-LAPI ISORROPIA II 2.3 iso_2_3Src directory")
if(MKPP_CMAQ_ISORROPIA_II_DIR OR MKPP_EPFL_ISORROPIA_II_DIR)
  enable_language(Fortran)
endif()

# Accept a sparse checkout or an exported source directory, with LF or CRLF.
function(mkpp_check_ii_source directory name expected)
  set(source "${directory}/${name}")
  set_property(
    DIRECTORY
    APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS "${source}")
  file(READ "${source}" content)
  string(REPLACE "\r\n" "\n" content "${content}")
  string(SHA256 digest "${content}")
  if(NOT "${digest}" STREQUAL "${expected}")
    message(FATAL_ERROR "ISORROPIA II reference ${source} differs from the pinned source")
  endif()
endfunction()

function(
  mkpp_add_ii_reference
  tag
  label
  directory
  common
  forward
  reverse)
  set(target "test_cracmm2_${tag}_ii")
  add_executable(
    ${target}
    "${CMAKE_CURRENT_LIST_DIR}/test_cmaq_ii.cpp" "${CMAKE_CURRENT_LIST_DIR}/cmaq_ii_reference.f90"
    "${directory}/${common}" "${directory}/${forward}" "${directory}/${reverse}")
  target_include_directories(${target} PRIVATE "${directory}")
  target_compile_definitions(
    ${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:MKPP_II_REFERENCE_NAME=\"${label}\">")
  target_compile_options(${target} PRIVATE $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=legacy>)
  target_link_libraries(${target} PRIVATE MKPP::cracmm2_inorganic)
  add_dependencies(cracmm2_inorganic_checks ${target})
  foreach(climate IN ITEMS cold warm)
    add_test(NAME cracmm2_${tag}_ii_${climate} COMMAND ${target} ${climate})
    set_tests_properties(cracmm2_${tag}_ii_${climate} PROPERTIES LABELS cracmm2_inorganic)
  endforeach()
endfunction()

if(MKPP_CMAQ_ISORROPIA_II_DIR)
  # USEPA/CMAQ 5.5.0.4, 1b70e0a395fc42cc34a4355807d39df1d3dbae84
  set(cmaq_dir "${MKPP_CMAQ_ISORROPIA_II_DIR}")
  mkpp_check_ii_source("${cmaq_dir}" isocom.f
                       4a7cabcac3e02ba14e426758d31dc7f54c20520f28c001d89d65d5c5fe70b054)
  mkpp_check_ii_source("${cmaq_dir}" isofwd.f
                       67b49a9bba1a4740b3f95e6cd0b2c77fb6da8c7f3ffa83201c1427f7911fdb0f)
  mkpp_check_ii_source("${cmaq_dir}" isorev.f
                       b0c966abbff7c3f5003c22281008490542c57d3932fa43cd42f3023878cc24ee)
  mkpp_check_ii_source("${cmaq_dir}" isrpia.inc
                       940017c28fd58074366e480b48a63c60e0fe87836ba359bd8456ffb42d17ce5d)
  mkpp_add_ii_reference(cmaq "CMAQ 5.5.0.4 II" "${cmaq_dir}" isocom.f isofwd.f isorev.f)
endif()

if(MKPP_EPFL_ISORROPIA_II_DIR)
  # EPFL-LAPI/ISORROPIA-II 2.3, 17102b249d9279d4565da3c652ea6ef57c3a3b34
  set(epfl_dir "${MKPP_EPFL_ISORROPIA_II_DIR}")
  mkpp_check_ii_source("${epfl_dir}" isocom.for
                       f2642903648e311ba7839ee27eee13ba17c6826cf93d1ea529b87ba3aeb33bde)
  mkpp_check_ii_source("${epfl_dir}" isofwd.FOR
                       7e5dad7ab5b94ec22afe56439fc0486f2e8ff84b5b5765e07de911218670b1f4)
  mkpp_check_ii_source("${epfl_dir}" isorev.for
                       f33d21c83790767a5dc2774d629500e685a7055485d7d8baec2aebfbdf2f2ff0)
  mkpp_check_ii_source("${epfl_dir}" isrpia.inc
                       eb8231851d5cc12333abef6d98f528d0560abdcff6c7f7af5d7c5a35391cdfed)
  mkpp_add_ii_reference(epfl "EPFL II 2.3" "${epfl_dir}" isocom.for isofwd.FOR isorev.for)
endif()
