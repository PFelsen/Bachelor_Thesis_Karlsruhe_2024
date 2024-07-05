#ifndef SPACETIME_DGVISCOACOUSTICDGDG_GLGL_HPP
#define SPACETIME_DGVISCOACOUSTICDGDG_GLGL_HPP

#include "STDGViscoAcousticAssemble.hpp"

class TDGViscoAcoustic_DGT_SystemOperator_GLGL : public Operator {
private:
  STDGViscoAcousticAssemble &assemble;
  const VectorMatrixBase &mg;
  const STDiscretization &disc;
  const AcousticProblem &problem;
  std::unique_ptr<SpaceTimeMatrixGraph> MG;


public:
  TDGViscoAcoustic_DGT_SystemOperator_GLGL(STDGViscoAcousticAssemble &assemble,
                                           const VectorMatrixBase &mg, bool assembleDiagonal = false)
      :
      assemble(assemble),
      mg(mg),
      disc(assemble.GetDisc()),
      problem(assemble.GetProblem()) {

  }

  void multiply_plus(Vector &b, const Vector &x) const override;

  std::string Info() const {
    return "Memory: 0 MB";
  }

};


#endif //SPACETIME_DGVISCOACOUSTICDGDG_GLGL_HPP
