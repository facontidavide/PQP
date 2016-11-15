project(PQP)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(./src SRC_LIST)

include_directories( ./src )

add_executable(demo_1 demos/sample/main.cpp ${SRC_LIST})
