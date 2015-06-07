include_directories(../src)

# Wing test
add_gtest(naca00XX_test naca00XX_test.cpp wing)

# UVLM test
add_gtest(morphing_test morphing_test.cpp morphing)
add_gtest(vortex_test vortex_test.cpp uvlm)
add_gtest(uvlm_vortex_ring_test uvlm_vortex_ring_test.cpp uvlm)
add_gtest(shed_test shed_test.cpp uvlm)
