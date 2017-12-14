cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR pandora)
project(${PROJECT_NAME_STR} C CXX)

#------------------
# set up
#------------------
find_package(Threads REQUIRED)
#set( CMAKE_VERBOSE_MAKEFILE on )
find_package (ZLIB)
if (NOT ZLIB_FOUND)
	message (FATAL_ERROR "zlib must be installed before configuration &
building can proceed")
endif()

# make sure that we have support for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -Wextra -O2 -g")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lz")

#if(CMAKE_COMPILER_IS_GNUCXX)
#    add_definitions(-Wall -pthread -Wextra -O0 -g)
#endif()

#if(MSVC)
#    #vc 2012 fix for vararg templates
#    set(MSVC_COMPILER_DEFS "-D_VARIADIC_MAX=10")
#endif() 

#-------------------
# set common include folder for module
#-------------------
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(COMMON_SRC ${PROJECT_SOURCE_DIR}/src)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(COMMON_TEST ${PROJECT_SOURCE_DIR}/test)

#-------------------
# Module source
#-------------------
include_directories(${COMMON_INCLUDES} ${COMMON_SRC})# ${EXT_PROJECTS_DIR}/FQFeeder/include)
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
#${EXT_PROJECTS_DIR}/FQFeeder/src/FastxParser.cpp) 
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARY})

#set(PROJECT_LIB_NAME ${PROJECT_NAME_STR}_impl)

#include_directories(${COMMON_INCLUDES})
#file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
#add_library(${PROJECT_LIB_NAME} ${SRC_FILES})

#-------------------
# Test
#-------------------
add_subdirectory(${EXT_PROJECTS_DIR}/gtest)

add_definitions(${MSVC_COMPILER_DEFS})

enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)

include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} ${COMMON_SRC}
${COMMON_TEST})# ${EXT_PROJECTS_DIR}/FQFeeder/include)

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp ${COMMON_SRC}/*.cpp)
#${EXT_PROJECTS_DIR}/FQFeeder/src/FastxParser.cpp)
list(REMOVE_ITEM TEST_SRC_FILES ${COMMON_SRC}/main.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)


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

target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARY})

add_test(testAll ${PROJECT_TEST_NAME})
