# OPTION A: Build an executable and a static library
	# Executable
	add_executable(DaMaSCUS
		main.cpp)

	target_compile_options(DaMaSCUS PUBLIC -Wall -pedantic)

	target_link_libraries(DaMaSCUS
		PUBLIC
			coverage_config
			lib_DaMaSCUS)

	target_include_directories(DaMaSCUS
		PRIVATE
			${GENERATED_DIR} )

	install(TARGETS DaMaSCUS DESTINATION ${BIN_DIR})

	# Static library
	add_library(lib_DaMaSCUS STATIC
		Data_Generation.cpp
		Earth_Model.cpp 
		Simulation_Trajectory.cpp
		Simulation_Utilities.cpp
		Underground_Distribution.cpp )

	set_target_properties(lib_DaMaSCUS PROPERTIES PREFIX "")

	target_compile_options(lib_DaMaSCUS PUBLIC -Wall -pedantic)

	target_link_libraries(lib_DaMaSCUS
		PUBLIC
			coverage_config
			libobscura
			${MPI_CXX_LIBRARIES})

	target_include_directories(lib_DaMaSCUS
		PRIVATE
			${GENERATED_DIR}
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}
			${INCLUDE_DIR})

	install(TARGETS lib_DaMaSCUS DESTINATION ${LIB_DIR})