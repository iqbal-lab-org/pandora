cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR pandora)
project(${PROJECT_NAME_STR} C CXX)

#------------------
# set up
#------------------
find_package(Threads REQUIRED)
#set( CMAKE_VERBOSE_MAKEFILE on )
find_package(ZLIB)
if (NOT ZLIB_FOUND)
    message(FATAL_ERROR "zlib must be installed before configuration &
building can proceed")
endif ()

# make sure that we have support for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else ()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -Wextra -DBOOST_SYSTEM_NO_DEPRECATED -O2 -g")

#-------------------
# Module source
#-------------------
include(${PROJECT_SOURCE_DIR}/ext/boost.cmake)

set(EXTERN_LIBS boost_filesystem boost_system)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src)
file(GLOB SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})
add_dependencies(${PROJECT_NAME} boost)
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC
        ${CMAKE_THREAD_LIBS_INIT}
        ${ZLIB_LIBRARY}
        ${EXTERN_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include)

enable_testing()
add_subdirectory(test)
