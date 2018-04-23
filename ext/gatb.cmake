include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/ext)

ExternalProject_Add(gatb
        DOWNLOAD_COMMAND  wget https://github.com/GATB/gatb-core/archive/v1.4.1.tar.gz --timestamping -O gatb.tar.gz
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/ext"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core-1.4.1/gatb-core && mkdir -p build && cd build && cmake .. && make"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(gatb extract_tar
        COMMAND tar -xvzf ${CMAKE_CURRENT_BINARY_DIR}/ext/gatb.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/ext
        DEPENDEES download
        DEPENDERS configure)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core-1.4.1/gatb-core/src
        ${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core-1.4.1/gatb-core/build/include)

link_directories(${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core-1.4.1/gatb-core/build/lib)