set(SOURCE_FILES
    "pluginMain.cpp"
    "cvWrapCmd.cpp"
    "cvWrapCmd.h"
    "cvWrapDeformer.cpp"
    "cvWrapDeformer.h"
    "bindingexporter.cpp"
    "bindingexporter.h"
    #"cvWrapRebindCmd.cpp"
    #"cvWrapRebindCmd.h"
    "common.cpp"
    "common.h"
    "cvwrap.cl"
)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MAYA_CXX_FLAGS} /arch:AVX")
find_package(Maya REQUIRED)

include_directories(${MAYA_INCLUDE_DIR})
link_directories(${MAYA_LIBRARY_DIR})
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${MAYA_LIBRARIES})

MAYA_PLUGIN(${PROJECT_NAME})
install(TARGETS ${PROJECT_NAME} ${MAYA_TARGET_TYPE} DESTINATION plug-ins)
install(FILES "cvwrap.cl" DESTINATION plug-ins)

