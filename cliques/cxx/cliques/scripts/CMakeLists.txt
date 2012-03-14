set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

if (UseFortran)
  set ( cliques_LIBRARIES expokit gfortran lapack blas ${Boost_LIBRARIES})
else()
  set ( cliques_LIBRARIES expokit lapack blas ${Boost_LIBRARIES})
endif()

set(CMAKE_CXX_FLAGS "-g -std=gnu++0x -Wall -DUSE_BOOST")
#set(CMAKE_CXX_FLAGS "-g -O3 -std=gnu++0x -Wall -DUSE_BOOST")

#add_executable(compare_robustness compare_robustness.cpp)
#target_link_libraries(compare_robustness ${Boost_LIBRARIES})

ADD_EXECUTABLE(create_landscape create_landscape.cpp)
target_link_libraries(create_landscape  ${cliques_LIBRARIES})

ADD_EXECUTABLE(create_communities create_communities.cpp)
target_link_libraries(create_communities  ${cliques_LIBRARIES})

ADD_EXECUTABLE(find_communities find_communities.cpp)
target_link_libraries(find_communities ${cliques_LIBRARIES})

ADD_EXECUTABLE(sample_landscape sample_landscape.cpp)
target_link_libraries(sample_landscape ${cliques_LIBRARIES})

ADD_EXECUTABLE(create_multilandscape create_multilandscape.cpp)
target_link_libraries(create_multilandscape ${cliques_LIBRARIES})