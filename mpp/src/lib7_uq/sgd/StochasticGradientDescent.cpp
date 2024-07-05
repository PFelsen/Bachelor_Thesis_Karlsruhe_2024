#include "StochasticGradientDescent.hpp"

double StochasticGradientDescent::Compute() {
  Config::PrintInfo();

  std::vector<double> u_k_l2_norm{};
  std::vector<double> j_k_l2_norm{};

  if (alpha < 1.0) {
    sum_of_stepsizes = std::make_unique<double>(0.0);
  }

  SampleID id(level, 0, false);

  MultiSampleFEMConfig conf(protocolName, Quantity);
  conf = conf.WithCoarseLevel(level - 1).WithFineLevel(level);

  auto msFEM = CreateSharedMSFEM(conf);
  SampleSolution gradient = msFEM->RunProtocol(id);

  SetStepsize();

  Vector vector = gradient.solution.vector;
  vector *= Stepsize();
  SampleSolution control(vector, id);

  if (StochasticGradientDescent::descentType == "ADAM") {
    m_k = std::make_unique<Vector>(0.0, control.solution.vector);
    v_k = std::make_unique<Vector>(0.0, control.solution.vector);
  }

  // start of Main Algorithm
  mout.StartBlock("SGD");
  while(!Breakpoint()) {
    vout(1) << "Iteration: " << Iteration() + 1 << endl;  //

    gradient = msFEM->RunProtocol(id, control);
    gradient = DescentEstimator(gradient);

    //    vout(6) << "D J(u_k, xi_k) = " << gradient_k.solution.vector << endl;

    control = Update(control, gradient);

    UpdateFunctionQOI(j_k_l2_norm, id, msFEM, control);

    //    vout(6) << "u_k+1 = " << u_k.solution.vector << endl;

    // Todo: Refactor the block after this comment.
    //  consider how to calculate the averaged coefficients in the loop.

    // AveragingSolution(control, u_T, Iteration());//for alpha==1 no effect
    // Redo with aggregate object

    if (plotting == 1) { PlotControl(control); }

    UpdateIteration();
    SetStepsize();
  }
  mout.EndBlock(verbose == 0);
  // end Main Algorithm

  // output

  mout.PrintInfo("SGD", 1, PrintInfoEntry("u_k ", vec2str(u_k_l2_norm), 1));

  mout.PrintInfo("SGD", 1,PrintInfoEntry("j(u_k) ", vec2str(j_k_l2_norm), 1));

  return j_k_l2_norm.back();
}

void StochasticGradientDescent::UpdateFunctionQOI(vector<double> &j_k_l2_norm, const SampleID &id, std::shared_ptr<MultiSampleFEM> &msFEM,
                                               const SampleSolution &control) const {
  if (overkill or Iteration() == maxSteps - 1) {
    double j_function_value;
    j_function_value = msFEM->RunProtocolQOI(id, control);
    j_k_l2_norm.push_back(j_function_value);
    vout(1) << "Estimated function value j: " << j_function_value << endl;
  }
}

void StochasticGradientDescent::PlotControl(const SampleSolution &control) const {
  //mpp::plot("control") << control.solution.vector << mpp::endp;
  mpp::PlotData("Control", "control." + control.IdString() + ".iteration." + to_string(Iteration()+1), control.solution.vector);
}

SampleSolution StochasticGradientDescent::AdmissibleMapping(
    SampleSolution &control) const {
  for (row r = control.solution.vector.rows();
       r != control.solution.vector.rows_end(); ++r) {
    if (control.solution.vector(r, 0) < u_a) {
      control.solution.vector(r, 0) = u_a;
    } else if (control.solution.vector(r, 0) > u_b) {
      control.solution.vector(r, 0) = u_b;
    }
  }
  return control;
}

SampleSolution StochasticGradientDescent::Update(
    SampleSolution &control, const SampleSolution &descent) const {  // Todo private machen
  control.solution.vector -= Stepsize() * descent.solution.vector;
  return AdmissibleMapping(control);
}

SampleSolution StochasticGradientDescent::DescentEstimator(
    SampleSolution &gradient) const {
  if (StochasticGradientDescent::descentType == "ADAM") {
    Vector g_hadarmard(gradient.solution.vector);
    Vector m_hat(gradient.solution.vector);
    Vector v_hat(gradient.solution.vector);
    Vector epsilon_vector(epsilonADAM, gradient.solution.vector);

    // calculation of ADAM step
    (*m_k) = beta1 * (*m_k);
    (*m_k) += (1 - beta1) * gradient.solution.vector;
    (*v_k) = beta2 * (*v_k);
    g_hadarmard = HadamardProductVector(gradient.solution.vector,
                                        gradient.solution.vector);
    (*v_k) += (1 - beta2) * g_hadarmard;

    m_hat = (1 / (1 - std::pow(beta1, Iteration() + 1))) * (*m_k);
    v_hat = (1 / (1 - std::pow(beta2, Iteration() + 1))) * (*v_k);
    v_hat = ComponentSqrt(v_hat) + epsilon_vector;
    gradient.solution.vector = ComponentDivide(m_hat, v_hat);
  }
  return gradient;
}

void StochasticGradientDescent::SetStepsize() {
  //Todo: Schrittweite direkt übergeben oder aus dem Kontext des Models auslesen (Geiersbach)
  //      Better to return the value directly?
  if (StochasticGradientDescent::descentType == "ADAM") {
    if (StochasticGradientDescent::stepsizeRule == "decreasing") {
      stepsize = gammaADAM / sqrt(Iteration() + 1);
      vout(4) << "stepsize (ADAM - d): " << stepsize << endl;
    } else {
      stepsize = gammaADAM;
      vout(4) << "stepsize (ADAM - c): " << stepsize << endl;
    }
  } else if (targetCost > 10e-5) {
    stepsize = theta / (Iteration() + 1.0 + nu);
    vout(4) << "stepsize (µ conv): " << stepsize << endl;
  } else if (StochasticGradientDescent::stepsizeRule == "decreasing") {
    stepsize = gamma / sqrt(Iteration() + 1.0);
    vout(4) << "stepsize (conv - d): " << stepsize << endl;
  } else if (StochasticGradientDescent::stepsizeRule == "constant") {
    stepsize = gamma / sqrt(maxSteps);
    vout(4) << "stepsize (conv - c): " << stepsize << endl;
  }
}


void StochasticGradientDescent::AveragingSolution(const SampleSolution &u_k,
                                                  SampleSolution &u_T) const {
  if (alpha < 1.0 && Iteration() >= maxSteps * alpha) {
    double sz = Stepsize();
    u_T.solution.vector += sz * u_k.solution.vector;
    (*sum_of_stepsizes) += sz;
  }
}

bool StochasticGradientDescent::Breakpoint() const {//Todo: implement for different breakpoint criteria
  if (Iteration()<maxSteps) return false;
  else return true;
}
