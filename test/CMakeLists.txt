include(${PROJECT_SOURCE_DIR}/ext/gtest.cmake)

file(GLOB TEST_SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp
        ${PROJECT_SOURCE_DIR}/test/*.cpp)
list(REMOVE_ITEM TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)

add_executable(pandora_test
        ${TEST_SRC_FILES})
add_dependencies(pandora_test
        gtest
        boost)

target_link_libraries(pandora_test
        libgmock
        libgtest
        ${ZLIB_LIBRARY})

target_include_directories(pandora_test PUBLIC
        ${CMAKE_BINARY_DIR}/include
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include
        ${PROJECT_SOURCE_DIR}/include)

add_test(NAME pandora_test COMMAND pandora_test)
