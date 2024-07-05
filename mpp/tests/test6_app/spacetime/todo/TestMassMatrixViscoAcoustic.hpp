#ifndef TESTMASSMATRIXVISCOACOUSTIC_HPP
#define TESTMASSMATRIXVISCOACOUSTIC_HPP

#include "DGViscoAcoustic.hpp"
#include "DebuggingTools.hpp"
#include "IMatrixGraph.hpp"
#include "TestEnvironment.hpp"
#include "random"

class TestMassMatrix : public Test {
public:
  Matrix *M = nullptr;

  void SetUp() override {

    std::map<string, string> configMap({
        {"Overlap", "dG2"},
        {"LinearSolver", "GMRES"},
        {"MinMu", "1e-3"},
        {"ProblemLevel", "2"},
        {"ModelImageKappa", "marmousi2-vp"},
        {"LinearPrintSteps", "100"},
        {"MaxRho", "2626.999855041504"},
        {"pml", "1.0"},
        {"doMeasure", "0"},
        {"Smoother", "PointBlockGaussSeidel"},
        {"ModelImageRho", "marmousi2-density"},
        {"LinearVerbose", "1"},
        {"plevel", "2"},
        {"Epsilon", "1e-10"},
        {"numL", "3"},
        {"LinearEpsilon", "1e-10"},
        {"LinearReduction", "1e-20"},
        {"BaseSolverEpsilon", "1e-6"},
        {"BaseSolverSteps", "10000"},
        {"LinearSteps", "1000"},
        {"Verbose", "10"},
        {"roi_min", "4.75, 0.1, 4.0"},
        {"BaseSolverReduction", "1e-2"},
        {"BaseSolver", "BiCGStab"},
        {"MaxMu", "2.802"},
        {"BasePreconditioner", "PointBlockGaussSeidel"},
        {"truncateSTMesh", "1"},
        {"source_duration", "0.1"},
        {"ProblemMid", "3.0, 0.25, 0.15"},
        {"MinRho", "1009.9992752075195"},
        {"linear_weight", "1"},
        {"MinKappa", "1.0279998779296875"},
        {"ConfigVerbose", "4"},
        {"MaxKappa", "4.7"},
        {"ModelImageMu", "marmousi2-vs"},
        {"DebugLevel", "0"},
        {"mlLength", "1.0"},
        {"useL2Proj", "0"},
    });

    Config::Initialize(configMap);


    SpaceTimeMeshes STM("ST_marmousi2_squares_abgeschnitten", 2, 3);
    Plot plot(STM.space_fine());
    Marmousi2_VA_rhs problem{};
    STDiscretization disc{1, 1, problem.problemDim(), STM.fine(), false};
    TDGViscoAcousticAssemble assemble{disc, problem, plot};
    /*SpaceTimeMatrixGraph cellmg(STM.fine(), disc,false, BlockDiagonal(), false);
    VectorMatrixBase mg(cellmg);
    M = new Matrix{mg};
    assemble.MassMatrix(*M);

    printFirstBlock(*M);


    auto S = Solver();
    S(*M);
    std::uniform_real_distribution<double> unif(-1,1);
    std::default_random_engine re(5);
    Vector v{mg};
    for (double* vp = v(); vp!= v()+v.size(); vp++){
        *vp = unif(re);
    }
    mout << norm(v) << endl;

    Vector Mv{(*M)*v};
    //mout << norm(vv) << endl;
    Vector SMv{S*Mv};
    mout << norm(SMv) << endl;
    Vector difference{v - SMv};
    mout << norm(difference) << endl;
    Vector Mtv{mg};
    M->multiply_transpose_plus(Mtv, v);
    difference = Mtv - Mv;
    mout << norm(difference) << endl;*/
    /*{
        SpaceTimeMatrixGraph cmg(STM.fine(), disc,false, false);
        VectorMatrixBase mg{cmg};
        Vector rhs{mg};
        Matrix M{mg};
        Matrix M0{mg};
        M = 0.0;
        M0 = 0.0;
        assemble.System(M0, rhs);
        Vector v{1.0,mg};
        Vector v2 = v;
        assemble.ApplySystemMatrixFree(M,v2);
        Vector tm = M*v;
        Vector ttm = tm + v2;
        Vector Mv = M0 * v;
        double diff = (Mv - ttm).norm();
        mout << diff << endl;


    }*/
  }

  void TearDown() override {}
};

TEST_F(TestMassMatrix, MassMatrix) {}


#endif // TESTMASSMATRIXVISCOACOUSTIC_HPP
