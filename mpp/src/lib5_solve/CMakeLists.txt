set(solve_src
        Assemble.cpp
        solver/preconditioner/PreconditionerBase.cpp
        solver/preconditioner/Jacobi.cpp
        solver/preconditioner/Richardson.cpp
        solver/preconditioner/GaussSeidel.cpp
        solver/preconditioner/SSOR.cpp
        solver/preconditioner/Multigrid.cpp
        solver/preconditioner/Preconditioner.cpp
        solver/preconditioner/CyclicPreconditioner.cpp
        solver/linearsolver/GMRES.cpp
        solver/linearsolver/LinearSolver.cpp
        solver/linearsolver/FullMultigrid.cpp
        solver/Newton.cpp
        timeintegrator/nonlinear/ImplicitEuler.cpp
        timeintegrator/TimeIntegrator.cpp
        timeintegrator/linear/RungeKutta.cpp
        timeintegrator/linear/ExponentialIntegrator.cpp
        timeintegrator/TimeSeries.cpp
        spectrum/eigensolver/IEigenSolver.cpp
        spectrum/eigensolver/RitzGalerkin.cpp
        spectrum/eigensolver/LOBPCG.cpp
        spectrum/eigensolver/LOBPCGSelective.cpp
        spectrum/Eigenpair.cpp
        )

if (USE_SUPERLU)
    set(solve_src ${solve_src} solver/preconditioner/SuperLU.cpp)
    set(solve_src ${solve_src} solver/preconditioner/ParallelSolverMatrix.cpp)
    set(solve_src ${solve_src} solver/preconditioner/ParallelSolverSteps.cpp)
    set(solve_src ${solve_src} solver/preconditioner/ParallelSolver.cpp)
endif ()


if (USE_CXSC)
    set(solve_src ${solve_src}
            spectrum/IAEigenvalueMethods.cpp
            )
endif ()

add_library(LIB_SOLVE SHARED ${solve_src})

target_link_libraries(LIB_SOLVE LIB_FEM)
