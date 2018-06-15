include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_65_1)

ExternalProject_Add(boost
        DOWNLOAD_COMMAND  wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/ext"
        CONFIGURE_COMMAND bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_65_1 && ./bootstrap.sh --with-libraries=filesystem,system,iostreams --prefix=${CMAKE_CURRENT_BINARY_DIR}"
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_65_1 && ./bjam install link=static"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(boost extract_tar
        COMMAND tar -xvzf ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_65_1.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/ext
        DEPENDEES download
        DEPENDERS configure)

