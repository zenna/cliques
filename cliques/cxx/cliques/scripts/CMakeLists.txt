find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

include_directories(../../)
include_directories(/usr/include/graphviz/)
link_directories(/usr/lib/graphviz/)

set(CMAKE_CXX_FLAGS "-O3 -Wall -DUSE_BOOST")

#add_executable(create_landscape create_landscape.cpp)
#target_link_libraries(create_landscape gvc)
#target_link_libraries(create_landscape graph)