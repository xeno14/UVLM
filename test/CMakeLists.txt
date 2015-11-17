include_directories(../src)

# Wing test
add_gtest(rect_test rect_test.cpp wing)
add_gtest(naca00XX_test naca00XX_test.cpp wing)
add_gtest(naca4digit_test naca4digit_test.cpp wing)
add_gtest(wing_test wing_test.cpp wing)

# UVLM test
add_gtest(morphing_test morphing_test.cpp morphing)
add_gtest(proto_adaptor_test proto_adaptor_test.cpp uvlm)
add_gtest(vortex_test vortex_test.cpp uvlm)
add_gtest(vortex_container_test vortex_container_test.cpp uvlm)
add_gtest(uvlm_vortex_ring_test uvlm_vortex_ring_test.cpp uvlm)
add_gtest(shed_test shed_test.cpp uvlm)
add_gtest(wing_builder_test wing_builder_test.cpp uvlm)
add_gtest(calc_load_test calc_load_test.cpp uvlm)
add_gtest(util_test util_test.cpp "")
add_gtest(linear_test linear_test.cpp uvlm)
add_gtest(calc_impulse_test calc_impulse_test.cpp uvlm)
