cmake_minimum_required(VERSION 3.15)

set(TARGET_NAME DelaunayDemo)

project(TARGET_NAME)

add_subdirectory(thirdparty)

add_subdirectory(delaunay)

add_subdirectory(test)