set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

set(CMAKE_CXX_FLAGS "-g -O3 -std=gnu++0x -Wall -DUSE_BOOST")

#add_executable(compare_robustness compare_robustness.cpp)
#target_link_libraries(compare_robustness ${Boost_LIBRARIES})


ADD_EXECUTABLE(find_communities find_communities.cpp)
target_link_libraries(find_communities expokit)
#target_link_libraries(find_communities gfortran)
target_link_libraries(find_communities blas)
target_link_libraries(find_communities lapack)
target_link_libraries(find_communities ${Boost_LIBRARIES})
