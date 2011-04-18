set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
find_package(Boost REQUIRED)
find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

include_directories(../../)

set(CMAKE_CXX_FLAGS "-g -pg -Wall -DUSE_BOOST")

add_executable(test_louvain louvain.cpp)	
target_link_libraries(test_louvain gvc)
target_link_libraries(test_louvain graph)

add_library(louvain_lib SHARED louvain.cpp)