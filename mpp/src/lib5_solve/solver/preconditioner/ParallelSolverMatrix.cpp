#include "ParallelSolverMatrix.hpp"
#include "ParallelSolverSteps.hpp"

using namespace std;

double ps_abs(double &d) {
  if (d < 0) return -d;
  return d;
}

ParallelSolverMatrix::ParallelSolverMatrix(const ParallelSolverOneStep &one_step, int matrixsize,
                                           int maxP) :
    cluster(one_step.get_own_cluster()), step(one_step), C(cluster.getCommunicationModule()),
    C_loc(0), vps(one_step.get_vps()), A(0), Right(0), Left(0), Schur(0), rhs(0), size_rhs(0),
    nrhs(1), info(0) {
  Schur = new ParallelMatrix(cluster.matrix_size_Schur(), cluster.matrix_size_Schur(), C,
                             matrixsize, maxP);
  Right = new ParallelMatrix(cluster.matrix_size_Sol(), *Schur, matrixsize, maxP);

  if (matrixsize == -1) {
    A = new ParallelMatrix(cluster.matrix_size_Sol(), 0, C, -1);      // no cols, but rows
    Left = new ParallelMatrix(cluster.matrix_size_Schur(), 0, C, -1); // no cols, but rows
                                                                      //         Left = new
    //         ParallelMatrix(cluster.matrix_size_Schur(),cluster.matrix_size_Sol(),C,matrixsize);
  } else {
    if (C.Size() < 3)
      A = new ParallelMatrix(cluster.matrix_size_Sol(), cluster.matrix_size_Sol(), C, -1);
    else
      A = new ParallelMatrix(cluster.matrix_size_Sol(), cluster.matrix_size_Sol(), C, matrixsize,
                             maxP);
    Left = new ParallelMatrix(cluster.matrix_size_Schur(), *A, matrixsize, maxP);
  }


  C_loc = A->get_PC_loc();
  if (Schur->get_dec_size() > A->get_dec_size()) C_loc = Schur->get_PC_loc();

  if (!C_loc) C_loc = &C;

  size_rhs = A->rows() + Schur->rows();
  //     Left = new ParallelMatrix(cluster.matrix_size_Schur(),cluster.matrix_size_Sol(),C);
}

void ParallelSolverMatrix::Set(const Scalar *a, const int *d, const int *col, int n) {
  if (cluster.matrix_size_Sol() + cluster.matrix_size_Schur() == 0) return;
  vector<int> searchindex_Sol;
  searchindex_Sol.resize(n);
  vector<int> searchindex_Schur;
  searchindex_Schur.resize(n);
  int akt = 0;
  for (unsigned int i = 0; i < searchindex_Sol.size(); ++i)
    searchindex_Sol[i] = -1;
  for (unsigned int i = 0; i < searchindex_Schur.size(); ++i)
    searchindex_Schur[i] = -1;
  for (int i = 0; i < cluster.size_Sol(); ++i) { // should be only i==0 in the first step...
    for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j, ++akt) {
      searchindex_Sol[(*vps.getvps(cluster.get_Sol_element(i)))[j]] = akt;
    }
  }

  akt = 0;
  for (int i = 0; i < cluster.size_Schur(); ++i) {
    for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j, ++akt) {
      searchindex_Schur[(*vps.getvps(cluster.get_Schur_element(i)))[j]] = akt;
    }
  }

  akt = 0; // set upper part
  for (int i = 0; i < cluster.size_Sol(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j, ++akt)
      for (int k = d[(*vps.getvps(cluster.get_Sol_element(i)))[j]];
           k < d[(*vps.getvps(cluster.get_Sol_element(i)))[j] + 1]; ++k) {
        int search_k = searchindex_Sol[col[k]];
        if (search_k != -1) {
          A->set(akt, search_k, a[k]);
        } else Right->set(akt, searchindex_Schur[col[k]], a[k]);
      }

  akt = 0; // set lower part
  for (int i = 0; i < cluster.size_Schur(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j, ++akt)
      for (int k = d[(*vps.getvps(cluster.get_Schur_element(i)))[j]];
           k < d[(*vps.getvps(cluster.get_Schur_element(i)))[j] + 1]; ++k) {
        int search_k = searchindex_Schur[col[k]];
        if (search_k != -1) {
          Schur->set(akt, search_k, a[k]);
        } else Left->set(akt, searchindex_Sol[col[k]], a[k]);
      }
}

void ParallelSolverMatrix::Communicate_other_dec(const ParallelSolverMatrix &next,
                                                 vector<pardec> &other_dec, int shift_own_proc,
                                                 int shift_other_proc, int &other_cluster_number) {
  // kommuniziere mit proc(0)<->other_proc(0)
  // danach Bcast mit aktuellem Communicator

  int own_size = 0;
  if (Schur->rows() != 0) own_size = Schur->get_dec_size();
  int cluster_number = step.cluster_number();
  int other_size = -1;

  if (shift_own_proc == 0)
    if (C.Proc() == 0) {
      next.C.Send(&own_size, 1, shift_other_proc);
      next.C.Receive(&other_size, 1, shift_other_proc);
    }
  if (shift_own_proc > 0)
    if (C.Proc() == 0) {
      next.C.Receive(&other_size, 1, shift_other_proc);
      next.C.Send(&own_size, 1, shift_other_proc);
    }

  C.Broadcast(&other_size, 1, 0);

  int *own_tmp = new int[3 * own_size + 1];
  int *other_tmp = new int[3 * other_size + 1];

  int j = 0;
  for (int i = 0; i < own_size; ++i, j += 3) {
    own_tmp[j] = Schur->get_dec(i).size;
    own_tmp[j + 1] = Schur->get_dec(i).proc;
    own_tmp[j + 2] = Schur->get_dec(i).shift;
  }
  own_tmp[j] = cluster_number;

  if (shift_own_proc == 0)
    if (C.Proc() == 0) {
      next.C.Send(own_tmp, 3 * own_size + 1, shift_other_proc);
      next.C.Receive(other_tmp, 3 * other_size + 1, shift_other_proc);
    }

  if (shift_own_proc > 0)
    if (C.Proc() == 0) {
      next.C.Receive(other_tmp, 3 * other_size + 1, shift_other_proc);
      next.C.Send(own_tmp, 3 * own_size + 1, shift_other_proc);
    }

  C.Broadcast(other_tmp, 3 * other_size + 1, 0);

  other_dec.resize(other_size);

  j = 0;
  for (unsigned int i = 0; i < other_dec.size(); ++i, j += 3) {
    other_dec[i].size = other_tmp[j];
    other_dec[i].proc = other_tmp[j + 1]; // !!!!!!
    other_dec[i].shift = other_tmp[j + 2];
  }
  other_cluster_number = other_tmp[j];

  delete[] own_tmp;
  delete[] other_tmp;
}

void ParallelSolverMatrix::Copy_Schur_Matrix(ParallelSolverMatrix &next) {
  next.Schur->Copy(*Schur);
  Schur->Clear();
}

void ParallelSolverMatrix::COPY(Scalar *data, size_t size, size_t numrow, Scalar *tmp) {
  memcpy(tmp, data, size * numrow * sizeof(Scalar));
}

void ParallelSolverMatrix::SEND(Communicator &Comm, Scalar *data, size_t size, size_t numrow,
                                int dest) {
  Comm.Send(data, size * numrow, dest);
}

void ParallelSolverMatrix::RECEIVE(Communicator &Comm, Scalar *tmp, size_t size, size_t numrow,
                                   int source) {
  Comm.Receive(tmp, size * numrow, source);
}

// current_shift: proc-shift!

void ParallelSolverMatrix::Set_LEFT(ParallelSolverMatrix &next,
                                    const ParallelSolverClusterStep *current_cluster,
                                    const Scalar *tmp, const int numrows, const int numcols,
                                    const int shift) {
  int r = 0;

  for (int ne = 0; ne < current_cluster->size_Schur(); ++ne) {
    int element = current_cluster->get_Schur_element(ne);
    int size = current_cluster->get_Schur_size(ne);
    int rshift = -1;
    for (int j = 0; j < next.cluster.size_Sol(); ++j)
      if (element == next.cluster.get_Sol_element(j)) {
        rshift = next.cluster.get_Sol_shift(j);
        for (int i = 0; i < size; ++i, ++r)
          for (int c = 0; c < numcols; ++c)
            next.A->loc_add(i + rshift, c + shift, *(tmp + r + c * numrows));
        break;
      }
    if (rshift == -1)
      for (int j = 0; j < next.cluster.size_Schur(); ++j)
        if (element == next.cluster.get_Schur_element(j)) {
          rshift = next.cluster.get_Schur_shift(j);
          for (int i = 0; i < size; ++i, ++r)
            for (int c = 0; c < numcols; ++c)
              next.Left->loc_add(i + rshift, c + shift, *(tmp + r + c * numrows));
          break;
        }
  }
}

void ParallelSolverMatrix::set_sr_dec_LEFT(ParallelSolverMatrix &next,
                                           const ParallelSolverClusterStep *current_cluster,
                                           const vector<pardec> *current_dec, int &proc_shift,
                                           int &shift_send, int &size, int &shift_receive) {
  if (size == 0) { return; }
  if ((*current_dec).size() == 0) { return; }

  int cur_i = 0;
  int cur_shift = 0;
  int next_i = 0;
  int next_shift = 0;

  vector<pardec> sending_dec;
  vector<pardec> receiving_dec;
  int sr_size = 0;

  while (size > 0) {
    while (cur_shift + (*current_dec)[cur_i].size <= shift_send) {
      cur_shift += (*current_dec)[cur_i].size;
      ++cur_i;
    }
    while (next_shift + (*next.A->get_dec())[next_i].size <= shift_receive) {
      next_shift += (*next.A->get_dec())[next_i].size;
      ++next_i;
    }

    int startshift_send = shift_send - cur_shift;
    int startshift_receive = shift_receive - next_shift;

    int sendsize = min((*next.A->get_dec())[next_i].size - startshift_receive,
                       (*current_dec)[cur_i].size - startshift_send);
    sendsize = min(sendsize, size);

    ++sr_size;

    sending_dec.resize(sr_size);
    receiving_dec.resize(sr_size);

    sending_dec[sr_size - 1].proc = (*current_dec)[cur_i].proc + proc_shift;
    sending_dec[sr_size - 1].size = sendsize;
    sending_dec[sr_size - 1].shift = startshift_send;

    receiving_dec[sr_size - 1].proc = (*next.A->get_dec())[next_i].proc;
    receiving_dec[sr_size - 1].size = sendsize;
    receiving_dec[sr_size - 1].shift = startshift_receive;


    size -= sendsize;
    shift_send += sendsize;
    shift_receive += sendsize;
  }

  size_t tmpsize = 0;
  tmpsize += next.A->rows() * next.A->cols_loc();
  tmpsize += next.Left->rows() * next.Left->cols_loc();
  Scalar *tmp = new Scalar[tmpsize];

  int rows = current_cluster->matrix_size_Schur();

  for (int i = 0; i < sr_size; ++i) {
    if (sending_dec[i].proc == receiving_dec[i].proc && sending_dec[i].proc == next.C.Proc()) {
      COPY(Schur->ref() + rows * sending_dec[i].shift, sending_dec[i].size, rows, tmp);
      Set_LEFT(next, current_cluster, tmp, rows, receiving_dec[i].size, receiving_dec[i].shift);
    } else if (next.C.Proc() == sending_dec[i].proc)
      SEND(next.C, Schur->ref() + Schur->rows() * sending_dec[i].shift, sending_dec[i].size, rows,
           receiving_dec[i].proc);
    else if (next.C.Proc() == receiving_dec[i].proc) {
      RECEIVE(next.C, tmp, sending_dec[i].size, rows, sending_dec[i].proc);
      Set_LEFT(next, current_cluster, tmp, rows, receiving_dec[i].size, receiving_dec[i].shift);
    }
  }
  delete[] tmp;
}

void ParallelSolverMatrix::SetNext_Matrix_LEFT_CLUSTERED(
    ParallelSolverMatrix &next, const ParallelSolverClusterStep *current_cluster,
    const vector<pardec> *current_dec, int &current_proc_shift) {
  for (int next_element = 0; next_element < next.cluster.size_Sol(); ++next_element) {
    int next_sol_element = next.cluster.get_Sol_element(next_element);
    int current_element = -1;
    for (int j = 0; j < current_cluster->size_Schur(); ++j) {
      if (next_sol_element == current_cluster->get_Schur_element(j)) {
        current_element = j;
        break;
      }
    }
    if (current_element == -1) continue;

    int shift_receive = next.cluster.get_Sol_shift(next_element);
    int shift_send = current_cluster->get_Schur_shift(current_element);
    int size = current_cluster->get_Schur_size(current_element);

    while (next.cluster.get_Sol_element(next_element + 1) != -1
           && next.cluster.get_Sol_element(next_element + 1)
                  == current_cluster->get_Schur_element(current_element + 1)) {
      next_element++;
      current_element++;
      size += current_cluster->get_Schur_size(current_element);
    }
    set_sr_dec_LEFT(next, current_cluster, current_dec, current_proc_shift, shift_send, size,
                    shift_receive);
  }
}

void ParallelSolverMatrix::Set_RIGHT(ParallelSolverMatrix &next,
                                     const ParallelSolverClusterStep *current_cluster,
                                     const Scalar *tmp, const int numrows, const int numcols,
                                     const int shift) {
  int r = 0;

  for (int ne = 0; ne < current_cluster->size_Schur(); ++ne) {
    int element = current_cluster->get_Schur_element(ne);
    int size = current_cluster->get_Schur_size(ne);
    int rshift = -1;
    for (int j = 0; j < next.cluster.size_Sol(); ++j)
      if (element == next.cluster.get_Sol_element(j)) {
        rshift = next.cluster.get_Sol_shift(j);
        for (int i = 0; i < size; ++i, ++r)
          for (int c = 0; c < numcols; ++c)
            next.Right->loc_add(i + rshift, c + shift, *(tmp + r + c * numrows));
        break;
      }
    if (rshift == -1)
      for (int j = 0; j < next.cluster.size_Schur(); ++j)
        if (element == next.cluster.get_Schur_element(j)) {
          rshift = next.cluster.get_Schur_shift(j);
          for (int i = 0; i < size; ++i, ++r)
            for (int c = 0; c < numcols; ++c)
              next.Schur->loc_add(i + rshift, c + shift, *(tmp + r + c * numrows));
          break;
        }
  }
}

void ParallelSolverMatrix::set_sr_dec_RIGHT(ParallelSolverMatrix &next,
                                            const ParallelSolverClusterStep *current_cluster,
                                            const vector<pardec> *current_dec, int &proc_shift,
                                            int &shift_send, int &size, int &shift_receive) {
  if (size == 0) { return; }
  if ((*current_dec).size() == 0) { return; }

  int cur_i = 0;
  int cur_shift = 0;
  int next_i = 0;
  int next_shift = 0;

  vector<pardec> sending_dec;
  vector<pardec> receiving_dec;
  int sr_size = 0;

  while (size > 0) {
    while (cur_shift + (*current_dec)[cur_i].size <= shift_send) {
      cur_shift += (*current_dec)[cur_i].size;
      ++cur_i;
    }
    while (next_shift + (*next.Schur->get_dec())[next_i].size <= shift_receive) {
      next_shift += (*next.Schur->get_dec())[next_i].size;
      ++next_i;
    }

    int startshift_send = shift_send - cur_shift;
    int startshift_receive = shift_receive - next_shift;

    int sendsize = min((*next.Schur->get_dec())[next_i].size - startshift_receive,
                       (*current_dec)[cur_i].size - startshift_send);
    sendsize = min(sendsize, size);

    ++sr_size;

    sending_dec.resize(sr_size);
    receiving_dec.resize(sr_size);

    sending_dec[sr_size - 1].proc = (*current_dec)[cur_i].proc + proc_shift;
    sending_dec[sr_size - 1].size = sendsize;
    sending_dec[sr_size - 1].shift = startshift_send;

    receiving_dec[sr_size - 1].proc = (*next.Schur->get_dec())[next_i].proc;
    receiving_dec[sr_size - 1].size = sendsize;
    receiving_dec[sr_size - 1].shift = startshift_receive;


    size -= sendsize;
    shift_send += sendsize;
    shift_receive += sendsize;
  }

  size_t tmpsize = 0;
  tmpsize += next.Schur->rows() * next.Schur->cols_loc();
  tmpsize += next.Right->rows() * next.Right->cols_loc();
  Scalar *tmp = new Scalar[tmpsize];

  int rows = current_cluster->matrix_size_Schur();

  for (int i = 0; i < sr_size; ++i) {
    if (sending_dec[i].proc == receiving_dec[i].proc && sending_dec[i].proc == next.C.Proc()) {
      COPY(Schur->ref() + rows * sending_dec[i].shift, sending_dec[i].size, rows, tmp);
      Set_RIGHT(next, current_cluster, tmp, rows, receiving_dec[i].size, receiving_dec[i].shift);
    } else if (next.C.Proc() == sending_dec[i].proc)
      SEND(next.C, Schur->ref() + rows * sending_dec[i].shift, sending_dec[i].size, rows,
           receiving_dec[i].proc);
    else if (next.C.Proc() == receiving_dec[i].proc) {
      RECEIVE(next.C, tmp, sending_dec[i].size, rows, sending_dec[i].proc);
      Set_RIGHT(next, current_cluster, tmp, rows, receiving_dec[i].size, receiving_dec[i].shift);
    }
  }
  delete[] tmp;
}

void ParallelSolverMatrix::SetNext_Matrix_RIGHT_CLUSTERED(
    ParallelSolverMatrix &next, const ParallelSolverClusterStep *current_cluster,
    const vector<pardec> *current_dec, int &current_proc_shift) {
  for (int next_element = 0; next_element < next.cluster.size_Schur(); ++next_element) {
    int next_sol_element = next.cluster.get_Schur_element(next_element);
    int current_element = -1;
    for (int j = 0; j < current_cluster->size_Schur(); ++j) {
      if (next_sol_element == current_cluster->get_Schur_element(j)) {
        current_element = j;
        break;
      }
    }
    if (current_element == -1) continue;

    int shift_receive = next.cluster.get_Schur_shift(next_element);
    int shift_send = current_cluster->get_Schur_shift(current_element);
    int size = current_cluster->get_Schur_size(current_element);

    while (next.cluster.get_Schur_element(next_element + 1) != -1
           && next.cluster.get_Schur_element(next_element + 1)
                  == current_cluster->get_Schur_element(current_element + 1)) {
      next_element++;
      current_element++;
      size += current_cluster->get_Schur_size(current_element);
    }
    set_sr_dec_RIGHT(next, current_cluster, current_dec, current_proc_shift, shift_send, size,
                     shift_receive);
  }
}

void ParallelSolverMatrix::SetNext_Matrix(ParallelSolverMatrix &next) {
  if (next.C.Size() == C.Size()) {
    Copy_Schur_Matrix(next);
    return;
  }

  int shift_own_proc = 0;
  if (next.C.Proc() > C.Proc()) shift_own_proc = next.C.Size() - C.Size();
  int shift_other_proc = 0;
  if (next.C.Proc() == C.Proc()) shift_other_proc = C.Size();

  vector<pardec> other_dec;
  int other_cluster_number = -1;

  Communicate_other_dec(next, other_dec, shift_own_proc, shift_other_proc, other_cluster_number);


  const ParallelSolverClusterStep *other_cluster = step.get_cluster(other_cluster_number);

  const ParallelSolverClusterStep *current_cluster;
  int current_proc_shift;
  const vector<pardec> *current_dec = NULL;
  const vector<pardec> pseudodec(0);

  if (shift_own_proc == 0) {
    current_cluster = &cluster;
    current_proc_shift = shift_own_proc;
    if (Schur->rows()) current_dec = Schur->get_dec();
    else current_dec = &pseudodec;
  } else {
    current_cluster = other_cluster;
    current_proc_shift = shift_other_proc;
    current_dec = &other_dec;
  }

  SetNext_Matrix_LEFT_CLUSTERED(next, current_cluster, current_dec, current_proc_shift);
  SetNext_Matrix_RIGHT_CLUSTERED(next, current_cluster, current_dec, current_proc_shift);

  //     do the same with the other cluster...
  if (shift_own_proc != 0) {
    current_cluster = &cluster;
    current_proc_shift = shift_own_proc;
    if (Schur->rows()) current_dec = Schur->get_dec();
    else current_dec = &pseudodec;
  } else {
    current_cluster = other_cluster;
    current_proc_shift = shift_other_proc;
    current_dec = &other_dec;
  }

  SetNext_Matrix_LEFT_CLUSTERED(next, current_cluster, current_dec, current_proc_shift);
  SetNext_Matrix_RIGHT_CLUSTERED(next, current_cluster, current_dec, current_proc_shift);

  Schur->Clear();
}

void ParallelSolverMatrix::makeLU() {
  if (A->rows() == 0) return;
  if (A->ref() == 0 && Schur->ref() == 0) return;
  int decsize = A->get_dec_size();

  for (int i = 0; i < decsize; ++i) {
    A->makeLU(i);
    if (C_loc->Proc() == C.Proc()) A->Communicate_Column_intern(i, C_loc);

    pardec dec = A->get_dec(i);

    if (A->ref() && C.Proc() > i) {
      A->Solve_own(i);
      A->MMM_own(i);
    }

    if (Schur->rows() != 0) {
      if (C_loc->Proc() == C.Proc()) Left->Communicate_Column_intern(i, C_loc, false);
      if (Left && C.Proc() > i) Left->MMM_left(A->ref() + dec.shift, A->rows(), dec);
      if (Right) {
        A->Solve_right(Right->ref(), Right->cols_loc(), dec);
        Right->MMM_right(A->get_other_column(), dec);
        Schur->MMM_schur(Left->get_other_column(), Right->ref(), Right->rows(), dec);
      }
    }
  }
}

void ParallelSolverMatrix::Set_rhs(const Scalar *u, int shift) {
  Scalar *t = rhs + shift * size_rhs;
  for (int i = 0; i < cluster.size_Sol(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j)
      *(t++) = u[(*vps.getvps(cluster.get_Sol_element(i)))[j]];
  for (int i = 0; i < cluster.size_Schur(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j)
      *(t++) = u[(*vps.getvps(cluster.get_Schur_element(i)))[j]];
}

void ParallelSolverMatrix::Create_rhs(int NRHS) {
  if (rhs) delete[] rhs;
  rhs = 0;
  nrhs = NRHS;
  rhs = new Scalar[size_rhs * nrhs];
  Scalar *t = rhs;
  for (int i = 0; i < size_rhs * nrhs; ++i)
    *(t++) = 0;
}

void ParallelSolverMatrix::Write_rhs(Scalar *u, int shift) {
  Scalar *t = rhs + shift * size_rhs;
  for (int i = 0; i < cluster.size_Sol(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j)
      u[(*vps.getvps(cluster.get_Sol_element(i)))[j]] = *(t++);
  for (int i = 0; i < cluster.size_Schur(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j)
      u[(*vps.getvps(cluster.get_Schur_element(i)))[j]] = *(t++);
}

void ParallelSolverMatrix::SetNext_rhs_LEFT(ParallelSolverMatrix &next) {
  Scalar *t = rhs;
  int shift_A = next.A->rows();
  int next_size = next.size_rhs;

  for (int n = 0; n < nrhs; ++n) {
    t += A->rows();
    for (int i = 0; i < cluster.size_Schur(); ++i) {
      int element = cluster.get_Schur_element(i);
      int size = cluster.get_Schur_size(i);
      int rshift = -1;
      for (int j = 0; j < next.cluster.size_Sol(); ++j)
        if (element == next.cluster.get_Sol_element(j)) {
          rshift = next.cluster.get_Sol_shift(j);
          for (int m = 0; m < size; ++m) {
            next.rhs[m + rshift + n * next_size] = *(t++);
          }
        }
      if (rshift == -1)
        for (int j = 0; j < next.cluster.size_Schur(); ++j)
          if (element == next.cluster.get_Schur_element(j)) {
            rshift = next.cluster.get_Schur_shift(j);
            for (int m = 0; m < size; ++m) {
              next.rhs[m + rshift + shift_A + n * next_size] = *(t++);
            }
          }
    }
  }
  if (C.Size() == next.C.Size()) return;
  size_t size_data = next.size_rhs * nrhs;
  Scalar *data = new Scalar[size_data];
  if (next.C.Proc() == 0) {
    SEND(next.C, next.rhs, size_data, 1, C.Size());
    RECEIVE(next.C, data, size_data, 1, C.Size());
  }

  if (C.Proc() == 0 && next.C.Proc() > 0) {
    RECEIVE(next.C, data, size_data, 1, 0);
    SEND(next.C, next.rhs, size_data, 1, 0);
  }
  C.Broadcast(data, size_data, 0);
  for (unsigned int i = 0; i < size_data; ++i)
    next.rhs[i] += data[i];

  delete[] data;
  // addiere von verschiedenen proc-mengen!
  // finde c-menge raus
  // addiere die beiden jeweils auf dem 0er
  // kommuniziere an alle...
}

void ParallelSolverMatrix::SetNext_rhs_RIGHT(ParallelSolverMatrix &next) {
  Scalar *t = next.rhs;
  int shift_A = A->rows();

  for (int n = 0; n < nrhs; ++n) {
    t += next.A->rows();
    for (int i = 0; i < next.cluster.size_Schur(); ++i) {
      int element = next.cluster.get_Schur_element(i);
      int size = next.cluster.get_Schur_size(i);
      int rshift = -1;
      for (int j = 0; j < cluster.size_Sol(); ++j)
        if (element == cluster.get_Sol_element(j)) {
          rshift = cluster.get_Sol_shift(j);
          for (int m = 0; m < size; ++m) {
            *(t++) = rhs[m + rshift + n * size_rhs];
          }
        }
      if (rshift == -1)
        for (int j = 0; j < cluster.size_Schur(); ++j)
          if (element == cluster.get_Schur_element(j)) {
            rshift = cluster.get_Schur_shift(j);
            for (int m = 0; m < size; ++m) {
              *(t++) = rhs[m + rshift + shift_A + n * size_rhs];
            }
          }
    }
  }
}

void ParallelSolverMatrix::SolveL() {
  if (A->rows() == 0) return;
  int decsize = A->get_dec_size();

  for (int i = 0; i < decsize; ++i) {
    if (C.Proc() == i) {
      A->Solve(i, rhs, nrhs, size_rhs);
      if (i != decsize - 1) A->MV_rhs(i, rhs, nrhs, size_rhs);
      if (Left->rows() != 0) { Left->MV_rhs_left(i, A->rows(), rhs, nrhs, size_rhs); }
      // SEND to i+1
      if (i < decsize - 1) SEND(C, rhs, nrhs, size_rhs, i + 1);
    }
    if (i < decsize - 1 && C.Proc() == i + 1) {
      // RECEIVE from i
      RECEIVE(C, rhs, nrhs, size_rhs, i);
    }
  }
  if (C.Proc() == C_loc->Proc()) C_loc->Broadcast(rhs, size_rhs * nrhs, decsize - 1);
}

void ParallelSolverMatrix::SolveU() {
  if (A->rows() == 0) return;
  int decsize = Right->get_dec_size();
  if (decsize > 0) {
    if (C.Proc() != decsize - 1)
      for (int K = 0; K < nrhs; ++K)
        for (int i = 0; i < Right->rows(); ++i)
          rhs[K * size_rhs + i] = 0;
    Right->MV_rhs_right(C.Proc(), rhs, nrhs, size_rhs);
    for (int K = 0; K < nrhs; ++K)
      C.Sum(rhs + K * size_rhs, Right->rows());
  }

  decsize = A->get_dec_size();
  if (decsize <= 1) {
    C.Broadcast(rhs, size_rhs * nrhs, 0);
    return;
  }
  for (int i = decsize - 1; i >= 1; --i) {
    if (C.Proc() == i) {
      A->MV_rhs_Z(i, rhs, nrhs, size_rhs);
      if (i > 1) SEND(C, rhs, nrhs, size_rhs, i - 1);
    }
    if (i > 1 && C.Proc() == i - 1) RECEIVE(C, rhs, nrhs, size_rhs, i);
  }
  if (decsize > 1) C.Broadcast(rhs, size_rhs * nrhs, 1);
}

ParallelSolverMatrix::~ParallelSolverMatrix() {
  C_loc = 0;
  Destruct();
  if (A) delete A;
  A = NULL;
  if (Schur) delete Schur;
  Schur = NULL;
  if (Left) delete Left;
  Left = NULL;
  if (Right) delete Right;
  Right = NULL;
  if (rhs) delete[] rhs;
  rhs = NULL;
}

int ParallelSolverMatrix::get_LU_procs() {
  if (A->cols() > 0) {
    if (A->get_PC_loc()) return A->get_PC_loc()->Size();
    else return A->get_PC().Size();
  } else return 0;
}

int ParallelSolverMatrix::get_Schur_procs() {
  if (Right->cols() > 0) {
    if (Right->get_PC_loc()) return Right->get_PC_loc()->Size();
    else return Right->get_PC().Size();
  } else return 0;
}

void ParallelSolverMatrix_Sparse::Set(const Scalar *Sa, const int *Sd, const int *Scol, int Sn) {
  if (cluster.matrix_size_Sol() + cluster.matrix_size_Schur() == 0) return;
  vector<int> searchindex_Sol;
  searchindex_Sol.resize(Sn);
  vector<int> searchindex_Schur;
  searchindex_Schur.resize(Sn);
  int akt = 0;
  for (unsigned int i = 0; i < searchindex_Sol.size(); ++i)
    searchindex_Sol[i] = -1;
  for (unsigned int i = 0; i < searchindex_Schur.size(); ++i)
    searchindex_Schur[i] = -1;
  for (int i = 0; i < cluster.size_Sol(); ++i) { // should be only i==0 in the first step...
    for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j, ++akt) {
      searchindex_Sol[(*vps.getvps(cluster.get_Sol_element(i)))[j]] = akt;
    }
  }

  akt = 0;
  for (int i = 0; i < cluster.size_Schur(); ++i) {
    for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j, ++akt) {
      searchindex_Schur[(*vps.getvps(cluster.get_Schur_element(i)))[j]] = akt;
    }
  }

  int nonzeroSparseA = 0;
  for (int i = 0; i < cluster.size_Sol(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j)
      for (int k = Sd[(*vps.getvps(cluster.get_Sol_element(i)))[j]];
           k < Sd[(*vps.getvps(cluster.get_Sol_element(i)))[j] + 1]; ++k) {
        if (searchindex_Sol[Scol[k]] != -1 && Sa[k] != 0.0) nonzeroSparseA++;
      }

  int nonzeroSparseLeft = 0;

  for (int i = 0; i < cluster.size_Schur(); ++i)
    for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j)
      for (int k = Sd[(*vps.getvps(cluster.get_Schur_element(i)))[j]];
           k < Sd[(*vps.getvps(cluster.get_Schur_element(i)))[j] + 1]; ++k) {
        if (searchindex_Sol[Scol[k]] != -1 && Sa[k] != 0.0) nonzeroSparseLeft++;
      }


  if (cluster.matrix_size_Sol() > 0) {
    SparseA = new SparseRMatrix(cluster.matrix_size_Sol(), nonzeroSparseA);

    int *d = (*SparseA).rowind();
    Scalar *nzval = (*SparseA).nzval();
    int *colptr = (*SparseA).colptr();
    int mi = 0;
    d[0] = 0;

    akt = 0; // set upper part
    for (int i = 0; i < cluster.size_Sol(); ++i)
      for (int j = 0; j < vps.getvps(cluster.get_Sol_element(i))->get_global_size(); ++j, ++akt) {
        d[j + 1] = d[j];
        for (int k = Sd[(*vps.getvps(cluster.get_Sol_element(i)))[j]];
             k < Sd[(*vps.getvps(cluster.get_Sol_element(i)))[j] + 1]; ++k) {
          int search_k = searchindex_Sol[Scol[k]];
          if (search_k != -1) {
            Scalar val = Sa[k];
            if (ps_abs(val) != 0) {
              colptr[mi] = search_k;
              nzval[mi] = val;
              ++d[j + 1];
              mi++;
            }
          } else Right->set(akt, searchindex_Schur[Scol[k]], Sa[k]);
        }
      }
  }
  if (cluster.matrix_size_Schur() > 0) {
    SparseLeft = new SparseRMatrix(cluster.matrix_size_Schur(), nonzeroSparseLeft);

    int *d = (*SparseLeft).rowind();
    Scalar *nzval = (*SparseLeft).nzval();
    int *colptr = (*SparseLeft).colptr();
    int mi = 0;
    d[0] = 0;

    akt = 0; // set lower part
    for (int i = 0; i < cluster.size_Schur(); ++i)
      for (int j = 0; j < vps.getvps(cluster.get_Schur_element(i))->get_global_size(); ++j, ++akt) {
        d[akt + 1] = d[akt];
        for (int k = Sd[(*vps.getvps(cluster.get_Schur_element(i)))[j]];
             k < Sd[(*vps.getvps(cluster.get_Schur_element(i)))[j] + 1]; ++k) {
          int search_k = searchindex_Schur[Scol[k]];
          if (search_k != -1) {
            Schur->set(akt, search_k, Sa[k]);
          } else {
            Scalar val = Sa[k];
            if (ps_abs(val) != 0) {
              colptr[mi] = searchindex_Sol[Scol[k]];
              nzval[mi] = val;
              ++d[akt + 1];
              mi++;
            }
          }
        }
      }
  }
}

void ParallelSolverMatrix_Sparse::makeLU() {
  if (!SparseA) return;
  S = GetSparseSolver(*SparseA);
  delete SparseA;
  SparseA = 0;

  Scalar *right = Right->ref();
  int local_cols = Right->cols_loc();
  int leftrows = Left->rows();

  S->Solve(right, local_cols);

  if (right && SparseLeft) {
    for (int i = 0; i < leftrows; ++i)
      for (int j = 0; j < leftrows; ++j) {
        Scalar s = 0;
        for (int k = SparseLeft->rowind(i); k < SparseLeft->rowind(i + 1); ++k)
          s += SparseLeft->nzval(k) * (*(*Right)(SparseLeft->colptr(k), j));
        *(*Schur)(i, j) -= s;
      }
  }
}

void ParallelSolverMatrix_Sparse::SolveL() {
  if (A->rows() == 0) return;
  for (int i = 0; i < nrhs; ++i)
    S->Solve(rhs + size_rhs * i);

  int row = A->rows();

  int leftrows = Left->rows();
  if (SparseLeft) {
    for (int n = 0; n < nrhs; ++n) {
      for (int i = 0; i < leftrows; ++i) {
        Scalar s = 0;
        for (int k = SparseLeft->rowind(i); k < SparseLeft->rowind(i + 1); ++k)
          s += SparseLeft->nzval(k) * rhs[SparseLeft->colptr(k) + n * size_rhs];
        rhs[i + row + n * size_rhs] -= s;
      }
    }
  }
}

void ParallelSolverMatrix_Sparse::SolveU() { ParallelSolverMatrix::SolveU(); }

void ParallelSolverMatrix_Sparse::Destruct() {
  if (SparseA) delete SparseA;
  SparseA = NULL;
  if (S) delete S;
  S = NULL;
  if (SparseLeft) delete SparseLeft;
  SparseLeft = NULL;
}
