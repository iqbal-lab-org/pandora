file(GLOB_RECURSE TEST_SRC_FILES
        ${SRC_FILES}
        ${PROJECT_SOURCE_DIR}/test/*.cpp)
list(REMOVE_ITEM TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)
add_executable(${PROJECT_NAME}_test ${TEST_SRC_FILES})

string(REPLACE "-fopenmp" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

target_include_directories(${PROJECT_NAME}_test PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/include
        )

target_link_libraries(${PROJECT_NAME}_test
        ${LIBARCHIVE_LIBRARY}
        ${LIBZIP_LIBRARY}
        ${RACON_LIBRARY}
        ${Gtest_LIBRARIES}
        ${GATB_LIBS}
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARY}
        ${CMAKE_DL_LIBS}
        ${STATIC_C_CXX}
        ${BACKWARD_LIBRARIES}
        ${RT_LIBRARY}
)

add_test(NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME}_test)
