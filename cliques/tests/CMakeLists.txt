include_directories(${BOOST_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-g -O3 -std=gnu++0x  -Wall -DUSE_BOOST")

if (UseFortran)
  set ( cliques_LIBRARIES expokit gfortran lapack blas ${Boost_LIBRARIES})
else()
  set ( cliques_LIBRARIES expokit lapack blas ${Boost_LIBRARIES})
endif()

ADD_EXECUTABLE(test_louvain louvain.cpp)
target_link_libraries(test_louvain ${cliques_LIBRARIES})

ADD_EXECUTABLE(louvain_pairs louvain_pairs.cpp)
target_link_libraries(louvain_pairs ${cliques_LIBRARIES})

add_executable(test_kernighanlin kernighanlin.cpp)
target_link_libraries(test_kernighanlin ${cliques_LIBRARIES})

add_executable(test_optimisation test_optimisation.cpp)
target_link_libraries(test_optimisation ${cliques_LIBRARIES})

#add_executable(test_edit_distance edit_distance.cpp)

#add_executable(test_varofinf varofinf.cpp)