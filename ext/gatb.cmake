include(ExternalProject)

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_INCLUDE_EXAMPLES  1)

ExternalProject_Add(gatb
        GIT_REPOSITORY https://github.com/GATB/gatb-core.git
        GIT_TAG "v1.4.1"
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/gatb"
        SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/gatb/src/gatb/gatb-core/gatb-core"
        CMAKE_ARGS -DKSIZE_LIST=32
        INSTALL_COMMAND "")

