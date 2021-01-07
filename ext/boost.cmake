include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_62_0)

ExternalProject_Add(boost
        DOWNLOAD_COMMAND  wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz --timestamping
        URL_HASH SHA1=34a745901533cef1a5c066199b50fe83368b961b
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/ext"
        CONFIGURE_COMMAND bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_62_0 && ./bootstrap.sh --with-libraries=system,filesystem,iostreams,log,thread,date_time --prefix=${CMAKE_CURRENT_BINARY_DIR}/boost"
        BUILD_COMMAND     ""
        INSTALL_COMMAND   bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_62_0 && ./b2 install link=static"
        TEST_COMMAND      ""
        )

ExternalProject_Add_Step(boost extract_tar
        COMMAND tar -xvzf ${CMAKE_CURRENT_BINARY_DIR}/ext/boost_1_62_0.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/ext
        DEPENDEES download
        DEPENDERS configure
        )
