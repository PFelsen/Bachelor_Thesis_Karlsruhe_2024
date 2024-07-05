#ifndef TESTMASSMATRIXVISCOACOUSTIC_HPP
#define TESTMASSMATRIXVISCOACOUSTIC_HPP

#include "DGViscoAcoustic.hpp"
#include "SpaceTimeMatrixGraph.hpp"
#include "TestEnvironment.hpp"

#include "DebuggingTools.hpp"

#include <random>

struct Data {
  std::shared_ptr<SpaceTimeMesh> stm;
  VectorMatrixBase *mg;
  TDGViscoAcousticAssemble *assemble;
  std::shared_ptr<STDiscretization> disc;
};

class TestRestriction : public Test {
public:
  Data getData(int sl, int tl, SpaceTimeMeshes &STM, TProblemViscoAcoustic &problem) {
    auto mesh = STM.getPlusPLevel(sl, tl);
    mout << "getData tcellcount: " << mesh->TCellCount() << endl;
    Plot plot(STM.space_fine());
    auto disc = std::make_shared<STDiscretization>(1, 1, problem.problemDim(), *mesh, false);

    SpaceTimeMatrixGraph cellmg(*mesh, *disc, false, false);
    auto *mg = new VectorMatrixBase(cellmg);
    auto *assemble = new TDGViscoAcousticAssemble{*disc, problem, plot};
    return {mesh, mg, assemble, disc};
  }

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


    SpaceTimeMeshes STM("ST_marmousi2_squares_abgeschnitten", 1, 2);

    Marmousi2_VA_rhs problem{};
    mout << "haha" << endl;
    Data d2 = getData(1, 0, STM, problem);
    Data d1 = getData(0, 0, STM, problem);

    mout << "haha2" << endl;
    Vector v(1.0, *d1.mg);

    {
      mout << "#####" << endl;
      mout << d1.stm->TCellCount() << " " << d1.mg->rowsize() << endl;
      mout << d2.stm->TCellCount() << " " << d2.mg->rowsize() << endl;
    }


    Vector v2{*d2.mg};
    auto ct = NewCellTransfer_2D(d2.disc, d1.disc);
    auto a = PointT<double, 2, 1>{-50, -50, 2};
    auto b = Point{50, 50, 2};

    mout << v.norm() << " " << endl;
    v = v * ct;
    mout << v2.norm() << " " << endl;


    // TDGViscoAcousticAssemble assemble{disc, problem, plot};
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

TEST_F(TestRestriction, TestRestriction) {}


#endif // TESTMASSMATRIXVISCOACOUSTIC_HPP
