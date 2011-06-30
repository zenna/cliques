find_package(GraphViz REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${graphviz_INCLUDE_DIRECTORIES})

link_directories(/home/mts09/compiler/GCC/libs/install/lib/)
SET (CMAKE_Fortran_COMPILER gfortran)
#SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
#set(CMAKE_CXX_COMPILER g++-4.6)
set(CMAKE_CXX_FLAGS "-g -std=gnu++0x -pg -Wall -DUSE_BOOST")

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
ENABLE_LANGUAGE(Fortran)
SET(F_OBJ expokit.f)

ADD_LIBRARY (expokit STATIC ${F_OBJ})
SET_TARGET_PROPERTIES(expokit PROPERTIES LINKER_LANGUAGE Fortran)

ADD_EXECUTABLE(test_louvain louvain.cpp)
target_link_libraries(test_louvain expokit)
target_link_libraries(test_louvain gfortran)
target_link_libraries(test_louvain blas)
target_link_libraries(test_louvain lapack)

#add_executable(test_louvain louvain.cpp expokit.o)	
#target_link_libraries(test_louvain ${graphviz_LIBRARIES})
#target_link_libraries(test_louvain expokit)

#add_executable(test_kernighanlin kernighanlin.cpp)	
#target_link_libraries(test_louvain ${graphviz_LIBRARIES})

#add_executable(test_edit_distance edit_distance.cpp)

#add_executable(test_varofinf varofinf.cpp)