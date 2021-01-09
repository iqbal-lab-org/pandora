file(GLOB_RECURSE TEST_SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*/*.cpp
        ${PROJECT_SOURCE_DIR}/test/*.cpp)
list(REMOVE_ITEM TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)
add_executable(${PROJECT_NAME}_test ${TEST_SRC_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

add_dependencies(${PROJECT_NAME}_test gatb gtest)

target_include_directories(${PROJECT_NAME}_test PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/include
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include
        )

target_link_libraries(${PROJECT_NAME}_test
        libgmock
        libgtest
        ${GATB_LIBS}
        ${BOOST_LIBRARIES}
        ${LIBRT}
        ${ZLIB_LIBRARY}
        ${CMAKE_DL_LIBS}
        -static-libgcc -static-libstdc++
        )

add_test(NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME}_test)
