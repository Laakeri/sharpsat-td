SET(flowcutter_SOURCES 
	cell.cpp 
	greedy_order.cpp
	list_graph.cpp
	tree_decomposition.cpp)


ADD_LIBRARY(libflowcutter ${flowcutter_SOURCES})