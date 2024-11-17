project(proj5)
set(CMAKE_CXX_STANDARD 20)
add_executable(proj5 main.cpp MATexception.cpp)


add_definitions(-D__NOEXCEPT_TRUE__)
add_definitions(-O3)

