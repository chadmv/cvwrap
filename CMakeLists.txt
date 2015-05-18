cmake_minimum_required(VERSION 2.6)
project(cvwrap)

set(PROJECT_PATH ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT})
set(CMAKE_INSTALL_PREFIX ${PROJECT_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cgcmake/modules)

add_subdirectory(src)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/module.txt" "${PROJECT_PATH}/${PROJECT_NAME}.txt")

install(DIRECTORY scripts DESTINATION .)
