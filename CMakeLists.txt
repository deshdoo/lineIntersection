cmake_minimum_required(VERSION 3.5)
project(lineintersection LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(app
    main.cpp
    geometry.cpp
    geometry.h
    uml.qmodel
)
