find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

include_directories(../../)
include_directories(/usr/include/graphviz/)
link_directories(/usr/lib/graphviz/)

set(CMAKE_CXX_FLAGS "-g -O3 -pg -Wall -DUSE_BOOST")

add_executable(draw_all_partitions draw_all_partitions.cpp)
target_link_libraries(draw_all_partitions gvc)
target_link_libraries(draw_all_partitions graph)

add_executable(graph_to_hsg graph_to_hsg.cpp)
target_link_libraries(graph_to_hsg gvc)
target_link_libraries(graph_to_hsg graph)

add_executable(test_bottleneck bottleneck.cpp)
add_executable(test_louvain louvain.cpp)	
target_link_libraries(test_louvain gvc)
target_link_libraries(test_louvain graph)

add_executable(test_varofinf varofinf.cpp)
add_executable(test_disjoint disjointset.cpp)
add_executable(test_maxima maxima.cpp)
add_executable(test_draw draw.cpp)
target_link_libraries(test_draw gvc)
target_link_libraries(test_draw graph)
