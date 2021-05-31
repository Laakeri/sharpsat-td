add_subdirectory (component_types)
add_subdirectory (preprocessor)
add_subdirectory (clhash)

SET(SOURCES alt_component_analyzer.cpp 
	component_cache.cpp
	component_management.cpp
	instance.cpp
	new_component_analyzer.cpp
	solver.cpp
	statistics.cpp
	decomposition.cpp)


ADD_LIBRARY(libsharpSAT ${SOURCES} component_types preprocessor clhash)
target_link_libraries (libsharpSAT component_types preprocessor clhash)
