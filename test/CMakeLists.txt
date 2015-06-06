include_directories(../src)

add_gtest(morphing_test morphing_test.cpp morphing)
add_gtest(naca00XX_test naca00XX_test.cpp naca00XXlib)
add_gtest(vortex_test vortex_test.cpp uvlm)
add_gtest(uvlm_vortex_ring_test uvlm_vortex_ring_test.cpp uvlm)
add_gtest(shed_test shed_test.cpp uvlm)
