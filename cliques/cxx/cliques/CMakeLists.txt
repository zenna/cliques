cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
#set(Boost_NO_SYSTEM_PATHS true)
set(BOOST_ROOT /usr/local/Cellar/boost/1.46.1)
set(BOOST_LIBRARYDIR /usr/local/Cellar/boost/1.46.1/lib)

find_package(Boost COMPONENTS program_options)

include_directories (../)
add_subdirectory (scripts)
add_subdirectory (tests)
add_subdirectory (algorithms)
add_subdirectory (structures)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")