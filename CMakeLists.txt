#-----------------
# cmake setup
#-----------------
cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR pandora)
project(${PROJECT_NAME_STR} C CXX)

find_package(Threads REQUIRED)

# make sure that we have support for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage -std=c++11")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# add definitions to g++
if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-Wall -ansi -Wno-deprecated -pthread --std=c++11)
endif()

#no idea what this does
#if(MSVC)
#    #vc 2012 fix for vararg templates
#    set(MSVC_COMPILER_DEFS "-D_VARIADIC_MAX=10")
#endif() 

#--------------------
# start project
#--------------------
#set(PROJECT_NAME_STR pandora)
#project(${PROJECT_NAME_STR} C CXX)

#-------------------
# set common include folder for module
#-------------------
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
# also add gtest subdir
add_subdirectory(${EXT_PROJECTS_DIR}/gtest)

#-------------------
# Module source
#-------------------
#set(PROJECT_LIB_NAME ${PROJECT_NAME_STR}_lib)

#include_directories(${COMMON_INCLUDES})
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB HDR_FILES ${COMMON_INCLUDES}/*.h)
#add_library(${PROJECT_LIB_NAME} ${SRC_FILES} ${HDR_FILES})


#-------------------
# Test
#-------------------
#add_definitions(${MSVC_COMPILER_DEFS})
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)

message(STATUS "test section COMMON_INCLUDES: ${COMMON_INCLUDES}")
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} ${SOURCE_DIR})

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES} ${SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)

# link to gtest library
if(NOT WIN32 OR MINGW)
    target_link_libraries(${PROJECT_TEST_NAME}
        ${GTEST_LIBS_DIR}/libgtest.a
        ${GTEST_LIBS_DIR}/libgtest_main.a
        )
else()
    target_link_libraries(${PROJECT_TEST_NAME}
        debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
        optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
        )
    target_link_libraries(${PROJECT_TEST_NAME}
        debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
        optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
        )
endif()

# link to thread thing
target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT})
#target_link_libraries(${PROJECT_TEST_NAME} ${PROJECT_LIB_NAME})

# run the test
add_test(test1 ${PROJECT_TEST_NAME})
