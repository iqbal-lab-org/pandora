cmake_minimum_required(VERSION 3.0)


# project configuration
set(PROJECT_NAME_STR pandora)
project(${PROJECT_NAME_STR} VERSION "0.7.0" LANGUAGES C CXX)
configure_file( include/version.h.in ${CMAKE_BINARY_DIR}/include/version.h )


# add a RELEASE_WITH_ASSERTS build type
set(CMAKE_CXX_FLAGS_RELEASE_WITH_ASSERTS "-O3")

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED -Wall -Wextra")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_OPENMP")
else ()
    find_package(OpenMP 4.0 REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()



########################################################################################################################
# FIND PRE-INSTALLED LIBS
########################################################################################################################
find_package(ZLIB REQUIRED)
########################################################################################################################
# END FIND PRE-INSTALLED LIBS
########################################################################################################################




########################################################################################################################
# EXTERNAL LIBS INSTALLATION
########################################################################################################################

########################################################################################################################
# INSTALL GATB
include(${PROJECT_SOURCE_DIR}/ext/gatb.cmake)
include_directories(SYSTEM
        ${gatb_source_dir}/gatb-core/src
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
include(${PROJECT_SOURCE_DIR}/ext/gtest.cmake)
include_directories(
        "${gtest_source_dir}/googletest/include"
        "${gtest_source_dir}/googlemock/include"
)
########################################################################################################################


########################################################################################################################
# INSTALL BOOST
include(${PROJECT_SOURCE_DIR}/ext/boost.cmake)
include_directories(${CMAKE_BINARY_DIR}/boost/include)
link_directories(${CMAKE_BINARY_DIR}/boost/lib)
set(BOOST_LIBS
        boost_atomic
        boost_chrono
        boost_date_time
        boost_filesystem
        boost_iostreams
        boost_log
        boost_log_setup
        boost_regex
        boost_system
        boost_thread)
set(Boost_USE_STATIC_LIBS ON)
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
        ${PROJECT_SOURCE_DIR}/cgranges/cpp
        )

#normal includes: warnings will be reported for these
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)

file(GLOB_RECURSE SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp
        ${PROJECT_SOURCE_DIR}/include/*/*.cpp
        ${PROJECT_SOURCE_DIR}/include/*.hpp
        ${PROJECT_SOURCE_DIR}/include/*.h
        ${PROJECT_SOURCE_DIR}/include/*/*.hpp
        ${PROJECT_SOURCE_DIR}/include/*/*.h
        )

add_executable(${PROJECT_NAME} ${SRC_FILES})

add_dependencies(boost gatb gtest)
add_dependencies(${PROJECT_NAME} gatb gtest boost)

target_link_libraries(${PROJECT_NAME}
        ${ZLIB_LIBRARY}
        ${BOOST_LIBS}
        ${GATB_LIBS}
        ${CMAKE_DL_LIBS}
)


enable_testing()
add_subdirectory(test)
