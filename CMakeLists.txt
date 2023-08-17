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
project(${PROJECT_NAME_STR} VERSION "0.11.0" LANGUAGES C CXX)
set(ADDITIONAL_VERSION_LABELS "")
configure_file( include/version.h.in ${CMAKE_BINARY_DIR}/include/version.h )

# C++14 required
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if (COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else ()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++14 support")
endif ()

# default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED -Wall -Wextra")

# compile with openmp only on Linux
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_OPENMP")
else ()
    find_package(OpenMP 4.0 REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# always link rt
set(RT_LIBRARY "rt")

if(NOT BIOCONDA)
    # static C and C++ flags
    set(STATIC_C_CXX "-static-libgcc -static-libstdc++")
endif()

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
# INSTALL libarchive
SET(ENABLE_TAR OFF CACHE BOOL "Enable tar building" FORCE)
SET(ENABLE_CPIO OFF CACHE BOOL "Enable cpio building" FORCE)
SET(ENABLE_CAT OFF CACHE BOOL "Enable cat building" FORCE)
SET(ENABLE_XATTR OFF CACHE BOOL "Enable extended attribute support" FORCE)
SET(ENABLE_ACL OFF CACHE BOOL "Enable ACL support" FORCE)
SET(ENABLE_ICONV OFF CACHE BOOL "Enable iconv support" FORCE)
SET(ENABLE_TEST OFF CACHE BOOL "Enable unit and regression tests" FORCE)
SET(ENABLE_WERROR OFF CACHE BOOL "Treat warnings as errors - default is ON for Debug, OFF otherwise." FORCE)
add_subdirectory(thirdparty/libarchive)
set(LIBARCHIVE_LIBRARY archive_static)
########################################################################################################################

########################################################################################################################
# INSTALL libzip
SET(ENABLE_COMMONCRYPTO OFF CACHE BOOL "Enable use of CommonCrypto" FORCE)
SET(ENABLE_GNUTLS OFF CACHE BOOL "Enable use of GnuTLS" FORCE)
SET(ENABLE_MBEDTLS OFF CACHE BOOL "Enable use of mbed TLS" FORCE)
SET(ENABLE_OPENSSL OFF CACHE BOOL "Enable use of OpenSSL" FORCE)
SET(ENABLE_WINDOWS_CRYPTO OFF CACHE BOOL "Enable use of Windows cryptography libraries" FORCE)
SET(ENABLE_BZIP2 OFF CACHE BOOL "Enable use of BZip2" FORCE)
SET(ENABLE_LZMA OFF CACHE BOOL "Enable use of LZMA" FORCE)
SET(ENABLE_ZSTD ON CACHE BOOL "Enable use of Zstandard" FORCE)
SET(BUILD_TOOLS OFF CACHE BOOL "Build tools in the src directory (zipcmp, zipmerge, ziptool)" FORCE)
SET(BUILD_REGRESS OFF CACHE BOOL "Build regression tests" FORCE)
SET(BUILD_EXAMPLES OFF CACHE BOOL "Build examples" FORCE)
SET(BUILD_DOC OFF CACHE BOOL "Build documentation" FORCE)
SET(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
SET(LIBZIP_DO_INSTALL OFF CACHE BOOL "Install libzip and the related files" FORCE)
SET(SHARED_LIB_VERSIONNING OFF CACHE BOOL "Add SO version in .so build" FORCE)
add_subdirectory(thirdparty/libzip)
set(LIBZIP_LIBRARY libzip::zip)
########################################################################################################################

########################################################################################################################
# CUSTOM CONFIG MINIMAP2
file(REMOVE_RECURSE ${PROJECT_SOURCE_DIR}/thirdparty/minimap2/lib)
########################################################################################################################

########################################################################################################################
# INSTALL Racon
add_subdirectory(thirdparty/racon)
set(RACON_LIBRARY racon)
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
hunter_add_package(Boost COMPONENTS filesystem iostreams log serialization system thread)
find_package(Boost CONFIG REQUIRED filesystem iostreams log serialization system thread)
set(Boost_LIBRARIES Boost::filesystem Boost::iostreams Boost::log Boost::serialization Boost::system Boost::thread)
########################################################################################################################
########################################################################################################################
# END EXTERNAL LIBS INSTALLATION
########################################################################################################################


########################################################################################################################
# PANDORA INSTALLATION
########################################################################################################################
# allows Seqan to be found
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/thirdparty/seqan/util/cmake")
set(SEQAN_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/thirdparty/seqan/include")

# Load the SeqAn module and fail if not found
find_package (SeqAn REQUIRED)

#include directories as SYSTEM includes, thus warnings will be ignored for these
include_directories(SYSTEM
        ${CMAKE_BINARY_DIR}/include
        ${PROJECT_SOURCE_DIR}/thirdparty/cgranges/cpp
        ${SEQAN_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}/thirdparty/minimap2/
        ${PROJECT_SOURCE_DIR}/thirdparty/libzip/lib
        ${PROJECT_SOURCE_DIR}/thirdparty/libarchive/libarchive
)

# normal includes: warnings will be reported for these
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/thirdparty/include
        ${PROJECT_SOURCE_DIR}/thirdparty/src
        ${PROJECT_SOURCE_DIR}/thirdparty/racon/src
)

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

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
        ${PROJECT_SOURCE_DIR}/thirdparty/minimap2/*.c
        )
list(REMOVE_ITEM SRC_FILES ${PROJECT_SOURCE_DIR}/thirdparty/minimap2/main.c ${PROJECT_SOURCE_DIR}/thirdparty/minimap2/example.c)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${PROJECT_NAME}
        ${LIBARCHIVE_LIBRARY}
        ${LIBZIP_LIBRARY}
        ${RACON_LIBRARY}
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARY}
        ${CMAKE_DL_LIBS}
        ${STATIC_C_CXX}
        ${SEQAN_LIBRARIES}
        ${RT_LIBRARY}
        )

# enable_testing()
# add_subdirectory(test)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
