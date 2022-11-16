cmake_minimum_required(VERSION 3.3)
project(Fingerprint)

include_directories(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
#add_executable(test_restoration tests/test_restoration.cpp)
find_package(OpenCV CONFIG REQUIRED PATHS "/usr/local/include/opencv4")

include_directories(${OpenCV_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(docs)
add_subdirectory(tests)
IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/data_out)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_out)
ENDIF()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/data_in)
    add_custom_target(files ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_in ${CMAKE_CURRENT_BINARY_DIR}/data_in)
ENDIF()

