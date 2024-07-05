add_executable(ConvertGeoFiles ${PROJECT_MPP_DIR}/tools/geometries/ConvertGeoFiles.cpp)
target_link_libraries(ConvertGeoFiles MPP_LIBRARIES stdc++fs)

if (CONVERT_VTU)
    message(STATUS "Building vtk toolbox")
    find_package(VTK REQUIRED)

    include(${VTK_USE_FILE})

    add_executable(ConvertVtu
            ${PROJECT_MPP_DIR}/tools/geometries/VtuConverter.cpp
            ${PROJECT_MPP_DIR}/tools/geometries/ConvertVtuFiles.cpp)
    target_link_libraries(ConvertVtu MPP_LIBRARIES ${VTK_LIBRARIES})
endif ()
