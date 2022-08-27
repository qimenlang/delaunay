add_executable(DelaunayTest test.cpp)

target_link_libraries(DelaunayTest PRIVATE delaunay Catch2WithMain)

enable_testing() 

add_test(Name delaunay_test
COMMAND DelaunayTest --success
)