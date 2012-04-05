find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

set(CMAKE_CXX_FLAGS "-g -std=gnu++0x -pg -Wall -DUSE_BOOST")

if (UseFortran)
  set ( cliques_LIBRARIES expokit gfortran lapack blas ${Boost_LIBRARIES})
else()
  set ( cliques_LIBRARIES expokit lapack blas ${Boost_LIBRARIES})
endif()

ADD_EXECUTABLE(test_louvain louvain.cpp)
target_link_libraries(test_louvain ${cliques_LIBRARIES})

add_executable(test_kernighanlin kernighanlin.cpp)
target_link_libraries(test_kernighanlin ${cliques_LIBRARIES})

add_executable(test_optimisation test_optimisation.cpp)
target_link_libraries(test_optimisation ${cliques_LIBRARIES})

#add_executable(test_edit_distance edit_distance.cpp)

#add_executable(test_varofinf varofinf.cpp)