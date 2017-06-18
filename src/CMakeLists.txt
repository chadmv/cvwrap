set(SOURCE_FILES
    "pluginMain.cpp"
    "cvWrapCmd.cpp"
    "cvWrapCmd.h"
    "cvWrapDeformer.cpp"
    "cvWrapDeformer.h"
    "bindingio.cpp"
    "bindingio.h"
    "common.cpp"
    "common.h"
    "cvwrap.cl"
)

if (WIN32)
    set(COMPILE_FLAGS "/arch:AVX")
else()
    set(COMPILE_FLAGS "-mavx")
endif()

find_package(Maya REQUIRED)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE Maya::Maya)
target_include_directories(${PROJECT_NAME} PRIVATE Maya::Maya)
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_FLAGS})
MAYA_PLUGIN(${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME} ${MAYA_TARGET_TYPE} DESTINATION plug-ins)
install(FILES "cvwrap.cl" DESTINATION plug-ins)

