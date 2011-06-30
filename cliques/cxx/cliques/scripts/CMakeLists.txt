set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})


#set(CMAKE_CXX_COMPILER g++-4.6)
set(CMAKE_CXX_FLAGS "-g -O3 -std=gnu++0x -Wall -llapack -lblas -DUSE_BOOST")

#add_executable(compare_robustness compare_robustness.cpp)
#target_link_libraries(compare_robustness ${Boost_LIBRARIES})

add_executable(create_landscape create_landscape.cpp)
target_link_libraries(create_landscape ${Boost_LIBRARIES})



ADD_EXECUTABLE(create_communities create_communities.cpp)
target_link_libraries(create_communities expokit)
target_link_libraries(create_communities gfortran)
target_link_libraries(create_communities blas)
target_link_libraries(create_communities lapack)

target_link_libraries(create_communities ${Boost_LIBRARIES})
