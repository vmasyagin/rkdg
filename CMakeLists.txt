cmake_minimum_required(VERSION 3.5)
project(rkdg)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

set(SOURCE_FILES rkdg.c)
add_executable(rkdg ${SOURCE_FILES})
target_link_libraries(rkdg m)