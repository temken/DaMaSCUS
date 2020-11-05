# Option A: If a static library is built. (see src/CMakeLists.txt)
add_executable(tests
	test_module.cpp)

target_link_libraries(tests 
	PRIVATE
		lib_DaMaSCUS
		gtest_main ) #contains the main function

target_include_directories(tests
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(tests PUBLIC -Wall -pedantic)

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME all_tests COMMAND tests
	WORKING_DIRECTORY ${TESTS_DIR} )