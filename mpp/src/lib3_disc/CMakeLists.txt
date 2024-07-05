set(disc_src
        dof/NodalPointProvider.cpp
        dof/IDoF.cpp
        dof/MixedDoF.cpp
        dof/RTDoF.cpp
        dof/DGDoF.cpp
        dof/EGDoF.cpp
        dof/SerendipityDoF.cpp
        dof/MulitpartDoF.cpp
        dof/ArgyrisDoF.cpp
        dof/CurlDoF.cpp
        dof/WCDoF.cpp
        shapes/Shapes.cpp
        quadrature/QuadratureSym.cpp
        quadrature/QuadratureGL.cpp
        quadrature/Quadrature.cpp
        matrixgraph/Row.cpp
        matrixgraph/IMatrixGraph.cpp
        matrixgraph/MatrixGraph.cpp
        matrixgraph/EGMatrixGraph.cpp
        matrixgraph/DGMatrixGraph.cpp
        matrixgraph/HybridMatrixGraph.cpp
        matrixgraph/FullMatrixGraph.cpp
        )

if (USE_SPACETIME)
    set (disc_src ${disc_src}
            matrixgraph/SpaceTimeMatrixGraph.cpp
            quadrature/SpaceTimeQuadrature.cpp
    )
endif()

add_library(LIB_DISC SHARED ${disc_src})
target_link_libraries(LIB_DISC LIB_MESH)
