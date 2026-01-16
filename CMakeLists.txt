cmake_minimum_required(VERSION 3.16)
project(StandaloneEEC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

find_package(Eigen3 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS MathMore GenVector Minuit2)
include(${ROOT_USE_FILE})

set(COMMON_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -pipe
        -Wpedantic
        -Wall
        -Wextra
        -Wundef
        -Wpointer-arith
        -Wcast-align
        -Wsign-compare
        -Wmissing-noreturn
        -Wmissing-format-attribute
        -Wunreachable-code
        -Wdisabled-optimization
        -Werror
        -Wno-suggest-attribute=format
    )
    # Apply release-only tuning without generator expressions (avoids literal `$<0:...` issues)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -flto")
    add_compile_definitions(NDEBUG)
    add_compile_options(-g3)
endif()

file(GLOB EEC_SOURCES SRothman/EECs/src/*.cc)
add_library(EEC SHARED ${EEC_SOURCES})
target_include_directories(EEC PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(EEC PUBLIC Eigen3::Eigen ${ROOT_LIBRARIES})

add_library(EEC_byhand SHARED ${EEC_SOURCES})
target_include_directories(EEC_byhand PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(EEC_byhand PUBLIC Eigen3::Eigen ${ROOT_LIBRARIES})
target_compile_definitions(EEC_byhand PUBLIC CHECK_BY_HAND)

file(GLOB SIMONTOOLS_SOURCES SRothman/SimonTools/src/*.cc)
add_library(SimonTools SHARED ${SIMONTOOLS_SOURCES})
target_include_directories(SimonTools PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(SimonTools PUBLIC Eigen3::Eigen ${ROOT_LIBRARIES})

file(GLOB MATCHING_SOURCES SRothman/Matching/src/*.cc)
add_library(Matching SHARED ${MATCHING_SOURCES})
target_include_directories(Matching PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(Matching PUBLIC Eigen3::Eigen ${ROOT_LIBRARIES})

set(ALL_LIBS EEC SimonTools Matching)
set(ALL_LIBS_BYHAND EEC_byhand SimonTools Matching)

function(add_eec_exe target_name source_file use_byhand)
    add_executable(${target_name} ${source_file})
    target_include_directories(${target_name} PRIVATE ${COMMON_INCLUDE_DIRS})
    if(${use_byhand})
        target_link_libraries(${target_name} PRIVATE ${ALL_LIBS_BYHAND})
        target_compile_definitions(${target_name} PRIVATE CHECK_BY_HAND)
    else()
        target_link_libraries(${target_name} PRIVATE ${ALL_LIBS})
    endif()
endfunction()

add_eec_exe(res4_benchmark drivers/res4_benchmark.cc OFF)
add_eec_exe(res4_unittest drivers/res4_unittest.cc OFF)
add_eec_exe(res4_checkbyhand drivers/res4_checkbyhand.cc ON)
add_eec_exe(CAres4_benchmark drivers/CAres4_benchmark.cc OFF)
add_eec_exe(CAres4_unittest drivers/CAres4_unittest.cc OFF)
add_eec_exe(CAres4_checkbyhand drivers/CAres4_checkbyhand.cc ON)
add_eec_exe(res3_benchmark drivers/res3_benchmark.cc OFF)
add_eec_exe(res3_unittest drivers/res3_unittest.cc OFF)
add_eec_exe(res3_checkbyhand drivers/res3_checkbyhand.cc ON)
add_eec_exe(CAres3_benchmark drivers/CAres3_benchmark.cc OFF)
add_eec_exe(CAres3_unittest drivers/CAres3_unittest.cc OFF)
add_eec_exe(CAres3_checkbyhand drivers/CAres3_checkbyhand.cc ON)
add_eec_exe(proj_unittest drivers/proj_unittest.cc OFF)
add_eec_exe(test_matching_v2 drivers/test_matching_v2.cc OFF)
add_eec_exe(testshower drivers/testshower.cc OFF)
