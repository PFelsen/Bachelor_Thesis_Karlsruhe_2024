#include "SpaceTimePreconditioner.hpp"


int STMultiGridPC::cycleNameToNumber(const string &cycle_name) {
  if (cycle_name == "V") return 1;
  else if (cycle_name == "W") return 2;
  else if (cycle_name == "VW") return 3;
  else if (cycle_name == "0") return 0;
  return 0;
}


STMultiGridPC::STMultiGridPC(std::unique_ptr<PathStrategy> pathStrategy,
                             const STAssemble &assemble,
                             const string &prefix,
                             std::shared_ptr<Preconditioner> precond,
                             std::shared_ptr<LinearSolver> base_solver)
    :
    Preconditioner(),
    pathStrategy(std::move(pathStrategy)),
    disc(assemble.GetSharedDisc()),
    assemble(assemble),
    BP(std::move(precond)),
    BS(std::move(base_solver)),
    dual(false) {


  Config::Get(prefix + "Smoother", name_smoother);
  Config::Get(prefix + "SmootherDamp_space", theta_space);
  Config::Get(prefix + "Presmoothing_space", pre_space);
  Config::Get(prefix + "Postsmoothing_space", pos_space);
  Config::Get(prefix + "SmootherDamp_time", theta_time);
  Config::Get(prefix + "Presmoothing_time", pre_time);
  Config::Get(prefix + "Postsmoothing_time", pos_time);
  Config::Get(prefix + "cycle_time", name_time);
  Config::Get(prefix + "cycle_space", name_space);
  Config::Get(prefix + "UseSparseMatrix", useSparseMatrix);
  Config::Get(prefix + "Transfer", name_transfer);

  Config::Get(prefix + "SkipBaseSolver", skipbasesolver);

  if (!BP) {
    string base_pc{"PointBlockGaussSeidel"};
    Config::Get(prefix + "BasePreconditioner", base_pc);
    BP = std::shared_ptr<Preconditioner>(GetPC(base_pc));
  }
  if (!BS) {
    string solver{"gmres"};
    Config::Get(prefix + "BaseSolver", solver);
    string basePreconditioner = "PointBlockGaussSeidel";
    Config::Get(prefix + "BasePreconditioner", basePreconditioner);
    BS = std::shared_ptr<LinearSolver>(
        GetLinearSolver(solver, GetPC(basePreconditioner), prefix + "BaseSolver"));
  }
  space_cycle = cycleNameToNumber(name_space);
  time_cycle = cycleNameToNumber(name_time);
}

void STMultiGridPC::Construct(const Matrix &_A) {
  this->Destruct();
  data.clear();
  Path path = pathStrategy->createPath(_A.Level(), _A.GetDisc().GetMeshes().PLevel());
  data.resize(path.size());
  for (int level = 0; level < path.size(); ++level) {
    int rlevel = path.size() - 1 - level;
    data[level].u = std::make_unique<Vector>(disc, path[rlevel]);
  }
  for (int level = 1; level < path.size(); ++level) {
    data[level].CT = std::make_unique<NewCellTransfer_2D>(*disc);
    data[level].transfer = GetSpaceTimeTransfer(*disc, name_transfer);
    int rlevel = path.size() - 1 - level;
    data[level].SM = std::unique_ptr<Preconditioner>(GetPC(name_smoother));
    if (path[rlevel + 1].time == path[rlevel].time) {
      data[level].pre = pre_space;
      data[level].pos = pos_space;
      data[level].theta = theta_space;
      data[level].cycle = space_cycle;
    } else if (path[rlevel + 1].space == path[rlevel].space) {
      data[level].pre = pre_time;
      data[level].pos = pos_time;
      data[level].theta = theta_time;
      data[level].cycle = time_cycle;
    }else {
        data[level].pre = pre_space;
        data[level].pos = pos_space;
        data[level].theta = theta_space;
        data[level].cycle = space_cycle;
    }
  }

  vout(1) << "Constructing MGPC at level" << _A.Level() << endl;
  for (int i = 0; i < data.size() - 1; ++i) {
    vout(3) << "level " << data[i].u->Level() << " Matrix start assemble " << endl;
    if (!data[i].A) {
      Matrix *A = new Matrix(*data[i].u);
      Vector v(0.0, *data[i].u);
      assemble.System(*A, v);
      data[i].A = A;
      if (dual) assemble.DualMatrix(*A, *data[i].A);
      if (useSparseMatrix) A->CreateSparse();
      vout(1) << "level " << data[i].u->Level() << " Matrix finish assemble " << endl;
    }
  }
  // TODO: Check this assignment
  data[data.size() - 1].A = &_A; //new Matrix(_A);
  (*BS)(*data[0].A);

  if (data.size() == 1) {
    (*BP).Construct(_A);
  }

  for (int level = 1; level < data.size(); ++level) {
    vout(10) << "     level " << data[level].u->Level() << " SM destruct " << Date() << endl;
    data[level].SM->Destruct();
    vout(10) << "     level " << data[level].u->Level() << " SM start " << Date() << endl;
    data[level].SM->Construct(*data[level].A);
    vout(10) << "     level " << data[level].u->Level() << " SM finish " << Date() << endl;
  }
}

void STMultiGridPC::Destruct() {
  for (int level = 0; level + 1 < data.size(); ++level) {
    if (data[level].A) delete data[level].A;
    data[level].A = nullptr;
  }
  for (int level = 0; level < data.size(); ++level) {
    if (data[level].SM) data[level].SM->Destruct();
  }
  data.clear();
}

void STMultiGridPC::Cycle(int level, Vector &u, Vector &r) const {
  if (level == 0) {
    if(!skipbasesolver) {
      vout(3) << "Solve on level: " << u.Level() << endl;
      u = (*BS) * r;

      //u.GetMesh().procSets.CheckConsistency();
      //u.GetMatrixGraph().GetProcSets().CheckConsistency();
      /*
      pout << DOUT(u) << endl;
      double l2BaseError = assemble.L2Error(u);
      mout << DOUT(l2BaseError) << endl;
      Vector utemp(0.0, u);
      assemble.get_exact_solution(utemp);
      double l2BaseExactError = assemble.L2Error(utemp);
      mout << DOUT(l2BaseExactError) << endl;

      Matrix temp(u);
      Vector rhs(0.0, u);
      assemble.System(temp, rhs);
      pout << r-rhs << endl;
      mout << norm(rhs) << endl;*/

    }
    return;
  }
  mout.StartBlock();
  Vector w(0.0, u);
  double norm_r_start = norm(r);
  for (int i = 0; i < data[level].pre; ++i) {
    w = (*data[level].SM) * r;
    w *= data[level].theta;
    r -= (*data[level].A) * w;
    u += w;
  }
  double norm_r_after_pre = norm(r);
  data[level].smoothing_steps += data[level].pre;

  vout(1) << "on level " << u.Level() << ": "
          << data[level].pre << " presmoothing (" << data[level].SM->Name()
          << ", damp = " << data[level].theta << ") steps done " << endl;
  Vector d(0.0, *data[level - 1].u);
  d.Clear();
  Vector c(0.0, *data[level - 1].u);

  mout.StartBlock("Restriction");
  data[level].transfer->Restrict(r, d);

  mout.EndBlock(verbose < 5);

  for (int i = 0; i < data[level].cycle; ++i) {
    Cycle(level - 1, c, d);
  }
  data[level].transfer->Prolongate(w, c);
  u += w;
  Vector temp = (*data[level].A) * w;
  r -= temp;

  double norm_r_after_mg = norm(r);

  mout.StartBlock("Postsmoothing");
  for (int i = 0; i < data[level].pos; ++i) {
    w = (*data[level].SM) * r;
    w *= data[level].theta;
    r -= (*data[level].A) * w;
    u += w;
  }
  mout.EndBlock(verbose < 5);
  vout(100) << std::fixed;
  vout(100) << "PRE -Correction: " << norm_r_start << " => " << norm_r_after_pre
            << " == " << norm_r_after_pre / norm_r_start << " %" << endl;
  vout(100) << "MG  -Correction: " << norm_r_after_pre << " => " << norm_r_after_mg
            << " == " << norm_r_after_mg / norm_r_after_pre << " %" << endl;
  vout(100) << "POST-Correction: " << norm_r_after_mg << " => " << norm(r)
            << " == " << norm(r) / norm_r_after_mg << " %" << endl;
  vout(100) << "All -Correction: " << norm_r_start << " => " << norm(r)
            << " == " << norm(r) / norm_r_start << " %" << endl;
  vout(100) << std::defaultfloat;
  data[level].smoothing_steps += data[level].pos;
  vout(1) << "on level " << u.Level() << ": "
          << data[level].pos << " postsmoothing (" << data[level].SM->Name()
          << ", damp = " << data[level].theta << ") steps done" << endl;

  mout.EndBlock(verbose < 100);
}

void STMultiGridPC::multiply(Vector &u, const Vector &b) const {
  if (data.size() == 1) {
    u = (*BP) * b;
    return;
  }

  Vector r(b);
  u = 0.0;
  Cycle(data.size() - 1, u, r);

  vout(5) << "cycle done" << endl;
}

void STMultiGridPC::transpose() const {
  Exit("Solve const Matrix")
  /*for (int level = 0; level < A.size() - 1; ++level) {
      assemble[level]->SystemAddDoubleD(*data[level].A);
      assemble[level]->TransposeMatrix(*data[level].A);
      (*data[level].A) *= -1.;
      //assemble[level]->SystemDual(*data[level].A);
      if (useSparseMatrix) data[level].A->CreateSparse();
  }
  (*BS)(*A[0]);
  for (int level = 1; level < A.size(); ++level) {
      CT[level]->dual = true;
      data[level].SM->Destruct();
      data[level].SM->Construct(*data[level].A);
  }*/
}

string STMultiGridPC::Name() const {
  return "MultilevelPreconditioner";
}

STMultiGridPC::~STMultiGridPC() {
  if (BP) BP->Destruct();
  Destruct();
}
