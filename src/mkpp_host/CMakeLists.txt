cmake_minimum_required(VERSION 3.24)
project(mkpp_host LANGUAGES CXX C)

# See the note in the MKPP root CMakeLists: a parent build that added MKPP via add_subdirectory
# already provides the Kokkos targets, and re-running find_package would resolve against a different
# installation.
if(NOT TARGET Kokkos::kokkos)
  find_package(Kokkos 5.1.1 REQUIRED)
endif()

add_library(mkpp_host INTERFACE)
add_library(MKPP::mkpp_host ALIAS mkpp_host)
target_compile_features(mkpp_host INTERFACE cxx_std_23)

target_include_directories(
  mkpp_host
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../../mkpp-generated>
            $<INSTALL_INTERFACE:include>)
target_link_libraries(mkpp_host INTERFACE Kokkos::kokkos)

if(MKPP_ENABLE_ADJOINT)
  target_compile_definitions(mkpp_host INTERFACE MKPP_ENABLE_ADJOINT)
endif()

if(MKPP_ENABLE_KOKKOS_KERNELS)
  if(TARGET KokkosKernels::kokkoskernels)
    target_link_libraries(mkpp_host INTERFACE KokkosKernels::kokkoskernels)
  elseif(TARGET kokkoskernels)
    target_link_libraries(mkpp_host INTERFACE kokkoskernels)
  endif()
  target_compile_definitions(mkpp_host INTERFACE MKPP_HAS_KOKKOS_KERNELS)
endif()

option(MKPP_ENABLE_CRACMM2_ISORROPIA_LITE
       "Build the host-side bulk CRACMM2 inorganic equilibrium adapter" OFF)
if(MKPP_ENABLE_CRACMM2_ISORROPIA_LITE)
  set(lite_dir "${CMAKE_CURRENT_SOURCE_DIR}/../isorropialite")
  if(NOT EXISTS "${lite_dir}/src/Solver.cpp")
    message(
      FATAL_ERROR "Initialize the pinned submodule: git submodule update --init src/isorropialite")
  endif()

  # Build the C ABI from the pinned solver sources. Its standalone project also configures a CLI, a
  # Fortran module and tests; none is needed by a host caller.
  add_library(
    mkpp_cracmm2_isorropia_lite STATIC
    "${lite_dir}/src/Solver.cpp" "${lite_dir}/src/Thermodynamics.cpp"
    "${lite_dir}/src/WaterActivities.cpp" "${lite_dir}/src/ForwardSolvers.cpp"
    "${lite_dir}/src/ActivityCoefficients_KM.cpp" "${lite_dir}/src/ReverseSolvers.cpp")
  add_library(MKPP::cracmm2_inorganic ALIAS mkpp_cracmm2_isorropia_lite)
  target_include_directories(mkpp_cracmm2_isorropia_lite PUBLIC "${lite_dir}/include")
  target_link_libraries(mkpp_cracmm2_isorropia_lite PUBLIC mkpp_host)
  set_target_properties(mkpp_cracmm2_isorropia_lite PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
