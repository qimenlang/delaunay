find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)

file(GLOB_RECURSE HEADER_FILES "*.h")
file(GLOB_RECURSE SOURCE_FILES "*.cpp")
message("headers: "${HEADER_FILES})
message("src: "${SOURCE_FILES})
add_executable(${TARGET_NAME} ${SOURCE_FILES}) 

target_link_libraries(${TARGET_NAME} PRIVATE sfml-system sfml-network sfml-graphics sfml-window)