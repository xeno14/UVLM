include_directories(../src)

add_gtest(morphing_test morphing_test.cpp morphing)
add_gtest(naca00XX_test naca00XX_test.cpp naca00XXlib)
add_gtest(vortex_test vortex_test.cpp vortex)
