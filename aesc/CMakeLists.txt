

cmake_minimum_required(VERSION 3.10)
project(aesc)


set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -DNDEBUG")

#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_DEBUG "-O0")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Boost REQUIRED filesystem system date_time program_options)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )



include_directories("${PROJECT_SOURCE_DIR}")
file(GLOB SOURCES "*.cc")
add_executable(aesc ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})