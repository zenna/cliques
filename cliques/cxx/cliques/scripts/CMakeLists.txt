set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
find_package(Boost REQUIRED)
find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

set(CMAKE_CXX_FLAGS "-O3 -Wall -DUSE_BOOST")

add_executable(create_landscape create_landscape.cpp)
target_link_libraries(create_landscape ${graphviz_LIBRARIES})