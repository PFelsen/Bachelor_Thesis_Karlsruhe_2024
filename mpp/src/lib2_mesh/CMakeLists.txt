set(cell
        cells/Celltype.cpp
        cells/Intval.cpp
        cells/ICell.cpp
        cells/Cell.cpp
        )

if (USE_SPACETIME) # todo rmv hack
    set(cell
            ${cell}
            cells/TCell.cpp
            )
endif()

if (SPACE_DIM MATCHES "^[2-3]$")
    set(cell
            ${cell}
            cells/Triangle.cpp
            cells/Quadrilateral.cpp
            )
endif ()
if (SPACE_DIM MATCHES "3")
    set(cell
            ${cell}
            cells/Tetrahedron.cpp
            cells/Hexahedron.cpp
            )
endif ()

set(mesh_src
        ${cell}
        distribution/Distribution.cpp
        mesh/Mesh.cpp
        mesh/Meshes.cpp
        mesh/MeshInfo.cpp
        mesh/MeshSettings.cpp
        mesh/CheckOrientation.cpp
        mesh/CoarseGeometry.cpp
        mesh/MeshesCreator.cpp
        meshparts/Face.cpp
        meshparts/Edge.cpp
        meshparts/Vertex.cpp
        meshparts/Identify.cpp
        meshparts/MeshPart.cpp
        meshparts/ProcSet.cpp
        )

add_library(LIB_MESH SHARED ${mesh_src})
target_link_libraries(LIB_MESH LIB_MATH)