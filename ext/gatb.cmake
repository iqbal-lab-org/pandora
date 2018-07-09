include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/ext)
execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/bin)
execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib)
execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/include)

set(GATB_DIR ${CMAKE_BINARY_DIR}/ext/gatb-core-1.4.1/gatb-core)
set(GATB_BUILD_DIR ${GATB_DIR}/build)
execute_process(COMMAND mkdir -p ${GATB_BUILD_DIR})

ExternalProject_Add(gatb
        DOWNLOAD_COMMAND  wget https://github.com/GATB/gatb-core/archive/v1.4.1.tar.gz --timestamping -O gatb.tar.gz
        DOWNLOAD_DIR      "${CMAKE_BINARY_DIR}/ext"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     bash -c "cd ${GATB_BUILD_DIR} && cmake ${GATB_DIR} && make"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(gatb extract_tar
        COMMAND tar -xvzf ${CMAKE_BINARY_DIR}/ext/gatb.tar.gz -C ${CMAKE_BINARY_DIR}/ext
        DEPENDEES download
        DEPENDERS configure)

ExternalProject_Add_Step(gatb copy
        COMMAND bash -c "cp -r ${GATB_BUILD_DIR}/bin/* ${CMAKE_BINARY_DIR}/bin/"
        COMMAND bash -c "cp -r ${GATB_BUILD_DIR}/lib/* ${CMAKE_BINARY_DIR}/lib/"
        COMMAND bash -c "cp -r ${GATB_BUILD_DIR}/include/* ${CMAKE_BINARY_DIR}/include/"
        COMMAND bash -c "cp -r ${GATB_DIR}/src/gatb/* ${CMAKE_BINARY_DIR}/include/gatb/"
        DEPENDEES build)
