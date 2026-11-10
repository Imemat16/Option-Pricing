cmake_minimum_required(VERSION 3.16)
project(OptionPricer LANGUAGES CXX)

add_library(option-pricing pricing.cpp)
target_compile_features(option-pricing PUBLIC cxx_std_17)

add_executable(option-pricer main.cpp)
target_link_libraries(option-pricer PRIVATE option-pricing)
install(TARGETS option-pricer RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
    add_executable(pricing-tests pricing_tests.cpp)
    target_link_libraries(pricing-tests PRIVATE option-pricing)
    add_test(NAME pricing COMMAND pricing-tests)
    add_test(NAME cli-example COMMAND option-pricer)
    set_tests_properties(cli-example PROPERTIES PASS_REGULAR_EXPRESSION "Black-Scholes: 10.450584")
endif()
