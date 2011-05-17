set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
find_package(Boost REQUIRED)
find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})
message(STATUS ${graphviz_LIBRARIES})
set(CMAKE_CXX_COMPILER g++-4.6)
set(CMAKE_CXX_FLAGS "-g -std=gnu++0x -Wall -llapack -DUSE_BOOST")


add_executable(sample_landscape sample_landscape.cpp)

#add_executable(create_landscape create_landscape.cpp)
#target_link_libraries(create_landscape ${graphviz_LIBRARIES})
