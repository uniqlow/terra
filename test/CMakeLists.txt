add_executable(terra_test main.cpp SphereTest.cpp EllipsoidTest.cpp)
add_test(NAME terra_test COMMAND terra_test)
