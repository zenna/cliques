set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-g -O3 -std=gnu++0x -Wall -llapack -lblas -DUSE_BOOST")

add_executable(sample_landscape sample_landscape.cpp)
target_link_libraries(sample_landscape ${Boost_LIBRARIES})

add_executable(create_landscape create_landscape.cpp)
target_link_libraries(create_landscape ${graphviz_LIBRARIES})