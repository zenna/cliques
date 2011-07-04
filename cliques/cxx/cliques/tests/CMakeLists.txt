find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

set(CMAKE_CXX_FLAGS "-g -std=gnu++0x -pg -Wall -DUSE_BOOST")

ADD_EXECUTABLE(test_louvain louvain.cpp)
target_link_libraries(test_louvain expokit)
target_link_libraries(test_louvain gfortran)
target_link_libraries(test_louvain lapack)
target_link_libraries(test_louvain blas)
target_link_libraries(test_louvain ${Boost_LIBRARIES})

#add_executable(test_louvain louvain.cpp expokit.o)	

#target_link_libraries(test_louvain ${graphviz_LIBRARIES})
#target_link_libraries(test_louvain expokit)

#add_executable(test_kernighanlin kernighanlin.cpp)	
#target_link_libraries(test_louvain ${graphviz_LIBRARIES})

#add_executable(test_edit_distance edit_distance.cpp)

#add_executable(test_varofinf varofinf.cpp)
