file(GLOB_RECURSE TEST_SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp
        ${PROJECT_SOURCE_DIR}/test/*.cpp)
list(REMOVE_ITEM TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)
add_executable(pandora_test ${TEST_SRC_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

add_dependencies(pandora_test
        gtest
        gatb
        boost)

target_link_libraries(pandora_test
        libgmock
        libgtest
        ${CMAKE_THREAD_LIBS_INIT}
        ${BOOST_LIBS}
        ${EXTERNAL_LIBS}
        ${CMAKE_DL_LIBS}
        ${ZLIB_LIBRARY})

target_include_directories(pandora_test PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/include
        ${source_dir}/googletest/include #Change here for CLion unit tests inspection
        ${source_dir}/googlemock/include)

add_test(NAME pandora_test COMMAND pandora_test)
