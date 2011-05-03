set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
find_package(Boost REQUIRED)
find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

set(CMAKE_CXX_FLAGS "-g -pg -Wall -DUSE_BOOST")

add_executable(test_louvain louvain.cpp)	
target_link_libraries(test_louvain ${graphviz_LIBRARIES})

add_executable(test_kernighanlin kernighanlin.cpp)	
target_link_libraries(test_louvain ${graphviz_LIBRARIES})

add_executable(test_bimap testbimap.cpp)	
target_link_libraries(test_louvain ${graphviz_LIBRARIES})