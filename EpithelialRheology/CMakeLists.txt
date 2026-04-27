cmake_minimum_required(VERSION 3.16)
project(EpithelialRheology CXX)

# ── System library search paths ────────────────────────────────────
link_directories(/usr/lib/x86_64-linux-gnu)

# ── Find Boost (needed by Chaste targets) ─────────────────────────
find_package(Boost REQUIRED COMPONENTS filesystem system serialization program_options)

# ── Find VTK (needed by Chaste targets) ───────────────────────────
find_package(VTK REQUIRED HINTS /opt/conda/lib/cmake/vtk-7.1)

# ── Find pre-built Chaste ──────────────────────────────────────────
set(Chaste_DIR "/opt/conda/build" CACHE PATH "Path to ChasteConfig.cmake")
find_package(Chaste COMPONENTS cell_based REQUIRED)

# ── CxxTest setup ──────────────────────────────────────────────────
set(CXXTEST_DIR "${Chaste_DIR}/cxxtest")
set(CXXTESTGEN "${CXXTEST_DIR}/cxxtestgen.py")

# ── Include paths ──────────────────────────────────────────────────
include_directories(SYSTEM ${Chaste_INCLUDE_DIRS})
include_directories(SYSTEM ${Chaste_THIRD_PARTY_INCLUDE_DIRS})
include_directories(SYSTEM ${CXXTEST_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# ── Project library (custom components) ────────────────────────────
file(GLOB PROJECT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
if(PROJECT_SOURCES)
    add_library(chaste_project_EpithelialRheology SHARED ${PROJECT_SOURCES})
    target_link_libraries(chaste_project_EpithelialRheology PUBLIC ${Chaste_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES})
    target_link_options(chaste_project_EpithelialRheology PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")
endif()

# ── Generate test runner from .hpp via cxxtestgen ──────────────────
set(TEST_HPP "${CMAKE_CURRENT_SOURCE_DIR}/test/TestEpithelialFillSimulation.hpp")
set(TEST_RUNNER "${CMAKE_CURRENT_BINARY_DIR}/TestEpithelialFillSimulation_runner.cpp")

add_custom_command(
    OUTPUT ${TEST_RUNNER}
    COMMAND ${Python3_EXECUTABLE} ${CXXTESTGEN} --error-printer -o ${TEST_RUNNER} ${TEST_HPP}
    DEPENDS ${TEST_HPP}
    COMMENT "Generating CxxTest runner"
)

# ── Build the test executable ──────────────────────────────────────
add_executable(TestEpithelialFillSimulation ${TEST_RUNNER})
target_include_directories(TestEpithelialFillSimulation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test
)
target_link_libraries(TestEpithelialFillSimulation
    chaste_project_EpithelialRheology
    ${Chaste_LIBRARIES}
    ${Chaste_THIRD_PARTY_LIBRARIES}
)
# Allow unresolved symbols in shared libs (VTK rendering/GL deps not needed at runtime)
target_link_options(TestEpithelialFillSimulation PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")

# ── Register with CTest ───────────────────────────────────────────
enable_testing()
add_test(NAME TestEpithelialFillSimulation COMMAND TestEpithelialFillSimulation)

# === Control: 15x4 Fill Simulation Test ===
set(TEST_HPP_CTRL "${CMAKE_CURRENT_SOURCE_DIR}/test/TestControl_15x4.hpp")
set(TEST_RUNNER_CTRL "${CMAKE_CURRENT_BINARY_DIR}/TestControl_15x4_runner.cpp")

add_custom_command(
    OUTPUT ${TEST_RUNNER_CTRL}
    COMMAND ${Python3_EXECUTABLE} ${CXXTESTGEN} --error-printer -o ${TEST_RUNNER_CTRL} ${TEST_HPP_CTRL}
    DEPENDS ${TEST_HPP_CTRL}
    COMMENT "Generating CxxTest runner for TestControl_15x4"
)

add_executable(TestControl_15x4 ${TEST_RUNNER_CTRL})
target_include_directories(TestControl_15x4 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test
)
target_link_libraries(TestControl_15x4
    chaste_project_EpithelialRheology
    ${Chaste_LIBRARIES}
    ${Chaste_THIRD_PARTY_LIBRARIES}
)
target_link_options(TestControl_15x4 PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")

add_test(NAME TestControl_15x4 COMMAND TestControl_15x4)

# === Proliferation parameter sweep ===
set(TEST_HPP_PROLIF "${CMAKE_CURRENT_SOURCE_DIR}/test/TestProliferation.hpp")
set(TEST_RUNNER_PROLIF "${CMAKE_CURRENT_BINARY_DIR}/TestProliferation_runner.cpp")

add_custom_command(
    OUTPUT ${TEST_RUNNER_PROLIF}
    COMMAND ${Python3_EXECUTABLE} ${CXXTESTGEN} --error-printer -o ${TEST_RUNNER_PROLIF} ${TEST_HPP_PROLIF}
    DEPENDS ${TEST_HPP_PROLIF}
    COMMENT "Generating CxxTest runner for TestProliferation"
)

add_executable(TestProliferation ${TEST_RUNNER_PROLIF})
target_include_directories(TestProliferation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test
)
target_link_libraries(TestProliferation
    chaste_project_EpithelialRheology
    ${Chaste_LIBRARIES}
    ${Chaste_THIRD_PARTY_LIBRARIES}
)
target_link_options(TestProliferation PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")

add_test(NAME TestProliferation COMMAND TestProliferation)

# === Deformation energy parameter sweep ===
set(TEST_HPP_DEFORM "${CMAKE_CURRENT_SOURCE_DIR}/test/TestDeformationEnergy.hpp")
set(TEST_RUNNER_DEFORM "${CMAKE_CURRENT_BINARY_DIR}/TestDeformationEnergy_runner.cpp")

add_custom_command(
    OUTPUT ${TEST_RUNNER_DEFORM}
    COMMAND ${Python3_EXECUTABLE} ${CXXTESTGEN} --error-printer -o ${TEST_RUNNER_DEFORM} ${TEST_HPP_DEFORM}
    DEPENDS ${TEST_HPP_DEFORM}
    COMMENT "Generating CxxTest runner for TestDeformationEnergy"
)

add_executable(TestDeformationEnergy ${TEST_RUNNER_DEFORM})
target_include_directories(TestDeformationEnergy PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test
)
target_link_libraries(TestDeformationEnergy
    chaste_project_EpithelialRheology
    ${Chaste_LIBRARIES}
    ${Chaste_THIRD_PARTY_LIBRARIES}
)
target_link_options(TestDeformationEnergy PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")

add_test(NAME TestDeformationEnergy COMMAND TestDeformationEnergy)

# === Adhesion (line tension) parameter sweep ===
set(TEST_HPP_ADHESION "${CMAKE_CURRENT_SOURCE_DIR}/test/TestAdhesion.hpp")
set(TEST_RUNNER_ADHESION "${CMAKE_CURRENT_BINARY_DIR}/TestAdhesion_runner.cpp")

add_custom_command(
    OUTPUT ${TEST_RUNNER_ADHESION}
    COMMAND ${Python3_EXECUTABLE} ${CXXTESTGEN} --error-printer -o ${TEST_RUNNER_ADHESION} ${TEST_HPP_ADHESION}
    DEPENDS ${TEST_HPP_ADHESION}
    COMMENT "Generating CxxTest runner for TestAdhesion"
)

add_executable(TestAdhesion ${TEST_RUNNER_ADHESION})
target_include_directories(TestAdhesion PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test
)
target_link_libraries(TestAdhesion
    chaste_project_EpithelialRheology
    ${Chaste_LIBRARIES}
    ${Chaste_THIRD_PARTY_LIBRARIES}
)
target_link_options(TestAdhesion PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")

add_test(NAME TestAdhesion COMMAND TestAdhesion)

# === BIOE242: 30x8 Fill Simulation Test ===
set(TEST_HPP_30x8 "${CMAKE_SOURCE_DIR}/../BIOE242/TestEpithelialFillSimulation_30x8.hpp")
set(TEST_RUNNER_30x8 "${CMAKE_BINARY_DIR}/TestEpithelialFillSimulation_30x8_runner.cpp")

add_custom_command(
    OUTPUT ${TEST_RUNNER_30x8}
    COMMAND ${Python3_EXECUTABLE} ${CXXTESTGEN} --error-printer -o ${TEST_RUNNER_30x8} ${TEST_HPP_30x8}
    DEPENDS ${TEST_HPP_30x8}
    COMMENT "Generating CxxTest runner for 30x8 test"
)

add_executable(TestEpithelialFillSimulation_30x8 ${TEST_RUNNER_30x8})
target_include_directories(TestEpithelialFillSimulation_30x8 PRIVATE
    ${CMAKE_SOURCE_DIR}/../BIOE242
)
target_link_libraries(TestEpithelialFillSimulation_30x8
    chaste_project_EpithelialRheology
    ${Chaste_LIBRARIES}
    ${Chaste_THIRD_PARTY_LIBRARIES}
)
target_link_options(TestEpithelialFillSimulation_30x8 PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")

add_test(NAME TestEpithelialFillSimulation_30x8 COMMAND TestEpithelialFillSimulation_30x8)
