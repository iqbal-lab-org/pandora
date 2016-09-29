add_library (Mini minimiser.cpp minimizer.h)
add_executable (demo main.cpp)
target_link_libraries (demo Mini)
