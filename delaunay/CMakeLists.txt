find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)

file(GLOB_RECURSE HEADER_FILES "*.h")
file(GLOB_RECURSE SOURCE_FILES "*.cpp")
message("headers: "${HEADER_FILES})
message("src: "${SOURCE_FILES})
add_library(delaunay delaunay.cpp edge.cpp triangle.cpp vector2.cpp)

add_executable(${TARGET_NAME} main.cpp)

target_link_libraries(${TARGET_NAME} PRIVATE delaunay sfml-system sfml-network sfml-graphics sfml-window)

