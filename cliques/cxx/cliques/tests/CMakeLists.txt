find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

include_directories(../../)
include_directories(/usr/include/graphviz/)
link_directories(/usr/lib/graphviz/)

set(CMAKE_CXX_FLAGS "-g -pg -Wall -DUSE_BOOST")

add_executable(test_louvain louvain.cpp)	
target_link_libraries(test_louvain gvc)
target_link_libraries(test_louvain graph)

