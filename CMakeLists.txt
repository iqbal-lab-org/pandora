cmake_minimum_required(VERSION 3.12)  # required by hunter ZLIB installation

# include hunter
option(HUNTER_STATUS_DEBUG "Hunter debug" OFF)  # comment if does not want hunter debug on
set(HUNTER_ROOT ${CMAKE_BINARY_DIR}/hunter)
include("cmake/HunterGate.cmake")
HunterGate(
        URL "https://github.com/cpp-pm/hunter/archive/v0.23.289.tar.gz"
        SHA1 "7d7323704780200a1575fc089d26f8a8a393a1fa"
        LOCAL
)

# project configuration
set(PROJECT_NAME_STR pandora)
project(${PROJECT_NAME_STR} VERSION "0.8.0" LANGUAGES C CXX)
configure_file( include/version.h.in ${CMAKE_BINARY_DIR}/include/version.h )

# add or not feature to print the stack trace
if(PRINT_STACKTRACE)
    message(STATUS "Printing meaningful stacktrace enabled, please have binutils-dev installed")
    add_compile_definitions(BACKWARD_HAS_BFD=1)
    set(BACKWARD_LIBRARIES "-lbfd")
else()
    set(BACKWARD_LIBRARIES "")
endif()


# C++11 required
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else ()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support")
endif ()

# default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED -Wall -Wextra")

if(PRINT_STACKTRACE)
    # -g is to add debug symbols, to make backtraces meaningful.
    # it does not impact performance (see https://stackoverflow.com/a/39223245)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif ()

# compile with openmp only on Linux
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_OPENMP")
else ()
    find_package(OpenMP 4.0 REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# static C and C++ flags
set(STATIC_C_CXX "-static-libgcc -static-libstdc++")

# always link rt
set(RT_LIBRARY "rt")

########################################################################################################################
# EXTERNAL LIBS INSTALLATION
########################################################################################################################
########################################################################################################################
# INSTALL ZLIB
hunter_add_package(ZLIB)
find_package(ZLIB 1.2.11 REQUIRED)
set(ZLIB_LIBRARY ZLIB::ZLIB)
########################################################################################################################

########################################################################################################################
# INSTALL GATB
include(${PROJECT_SOURCE_DIR}/ext/gatb.cmake)
include_directories(SYSTEM
        ${gatb_source_dir}/src
        ${gatb_binary_dir}/include
        )
link_directories(${gatb_binary_dir}/lib)
set(GATB_LIBS
        gatbcore
        hdf5
        hdf5_tools
    )
########################################################################################################################

########################################################################################################################
# INSTALL GTEST
hunter_add_package(GTest)
find_package(GTest CONFIG REQUIRED)
set(Gtest_LIBRARIES GTest::gtest GTest::gmock_main)
########################################################################################################################

########################################################################################################################
# INSTALL BOOST
set(Boost_USE_STATIC_LIBS ON)

if(BIOCONDA)
    find_package(Boost REQUIRED filesystem iostreams log log_setup system thread)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    hunter_add_package(Boost COMPONENTS filesystem iostreams log system thread)
    find_package(Boost CONFIG REQUIRED filesystem iostreams log log_setup system thread)
    set(Boost_LIBRARIES Boost::filesystem Boost::iostreams Boost::log Boost::log_setup Boost::system Boost::thread)
endif()

########################################################################################################################
########################################################################################################################
# END EXTERNAL LIBS INSTALLATION
########################################################################################################################


########################################################################################################################
# PANDORA INSTALLATION
########################################################################################################################
#include directories as SYSTEM includes, thus warnings will be ignored for these
include_directories(SYSTEM
        ${CMAKE_BINARY_DIR}/include
        ${PROJECT_SOURCE_DIR}/thirdparty/cgranges/cpp
        )

# normal includes: warnings will be reported for these
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/thirdparty/include
        ${PROJECT_SOURCE_DIR}/thirdparty/src
)

file(GLOB_RECURSE SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp
        ${PROJECT_SOURCE_DIR}/include/*/*.cpp
        ${PROJECT_SOURCE_DIR}/include/*.hpp
        ${PROJECT_SOURCE_DIR}/include/*.h
        ${PROJECT_SOURCE_DIR}/include/*/*.hpp
        ${PROJECT_SOURCE_DIR}/include/*/*.h
        ${PROJECT_SOURCE_DIR}/thirdparty/src/*.cpp
        ${PROJECT_SOURCE_DIR}/thirdparty/include/*.hpp
        ${PROJECT_SOURCE_DIR}/thirdparty/include/*.h
        )

add_executable(${PROJECT_NAME} ${SRC_FILES})
add_dependencies(${PROJECT_NAME} gatb)

target_link_libraries(${PROJECT_NAME}
        ${GATB_LIBS}
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARY}
        ${CMAKE_DL_LIBS}
        ${STATIC_C_CXX}
        ${BACKWARD_LIBRARIES}
        ${RT_LIBRARY}
)

enable_testing()
add_subdirectory(test)
