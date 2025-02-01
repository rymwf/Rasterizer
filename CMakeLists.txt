cmake_minimum_required(VERSION 3.20)

project(rasterizer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)


file(GLOB commonsrcs common/*.cpp)

#======================
file(GLOB srcs examples/*.cpp)

foreach(example ${srcs}) 

get_filename_component(a ${example} NAME_WLE)

add_executable(${a} ${example})

target_include_directories(${a} PRIVATE 
${CMAKE_SOURCE_DIR}/common
${CMAKE_SOURCE_DIR}/3rdparty
)

endforeach()