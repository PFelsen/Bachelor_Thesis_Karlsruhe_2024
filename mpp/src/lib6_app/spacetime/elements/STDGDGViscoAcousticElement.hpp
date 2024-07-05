#ifndef SPACETIME_STDGDGVISCOACOUSTICELEMENT_HPP
#define SPACETIME_STDGDGVISCOACOUSTICELEMENT_HPP

#include "SpaceTimeDiscretization.hpp"
#include "SpaceTimeQuadrature.hpp"
#include "SpaceTimeShape.hpp"
#include "Vector.hpp"
#include  "AcousticProblems.hpp"

const SpaceTimeShape &GetSpaceTimeShape(const VectorMatrixBase &base,
                                        const Cell &c);

const SpaceTimeShape &GetSpaceTimeShapeHighOrder(const VectorMatrixBase &base,
                                                 const Cell &c);

const SpaceTimeQuadrature &GetSpaceTimeQuadrature(const VectorMatrixBase &base,
                                                  const Cell &c);

const SpaceTimeQuadrature &GetSpaceTimeQuadratureHighOrder(
    const VectorMatrixBase &base);

class STDGDGViscoAcousticElement {
public:
  double det_space;
  double det_time;

  std::vector<double> space_qWeight;
  std::vector<double> time_qWeight;
  std::vector<double> qWeight;
  std::vector<Point> qPoint;

  const int dimension;
  const SpaceTimeShape &shape;
  const SpaceTimeQuadrature &quad;
  const Quadrature &space_quad;
  const cell &TC;
  const row r;

  ShapeValues<VectorFieldT<double, SpaceDimension, 1>> gradients;
  std::vector<VectorField> velocityST;
  std::vector<VectorField> dtVelocityST;
  std::vector<double> divVelocityST;

  std::vector<double> pressureST;
  std::vector<double> dtPressureST;
  std::vector<VectorField> gradPressureST;

  int dampingComponentCount;

  std::vector<std::vector<int>> component_to_index;

  STDGDGViscoAcousticElement(const VectorMatrixBase &, const cell &, int nL = 0,
                             bool max_quad_order = false);

  inline size_t GetComponentCount() const {
    return (1 + dimension + dampingComponentCount);
  }

  inline size_t GetFullDimensions() const {
    return shape.size() * GetComponentCount();
  }

  inline size_t GetQuadratureSize() const { return quad.size(); }

  int nQ() const;

  int nSpaceQ() const;

  int nTimeQ() const;

  double QWeight(int q) const;

  double SpaceQWeight(int qq);

  Point SpaceQPoint(int q) const;

  Point TimeQPoint(int q) const;

  double TimeQWeight(int qq);

  Point LocalTimeQPoint(int q) const;

  const Point &QPoint(int q) const;

  double Area() const;

  int i_dimension() const;

  int j_dimension() const;

  int variable(int j) const;

  int GetComponent(int j) const;

  std::vector<int> GetIndicesForComponent(int component) const {
    return component_to_index[component];
  }

  double EvaluateComponentLocal(const Point &p, const Vector &u,
                                int component) const;

  double EvaluateComponentGlobal(const Point &p, const Vector &u,
                                 int component) const;

  double EvaluateComponent(int nq, const Vector &u, int component) const;

  VectorField Velocity(int nq, int j) const;

  VectorField Velocity(int nq, const Vector &U) const;

  VectorField VelocityLocal(const Point &localPoint, int j) const;

  VectorField VelocityLocal(const Point &localPoint, const Vector &u) const;

  VectorField VelocityGlobal(const Point &globalPoint, const Vector &u) const;

  VectorField DtVelocity(int nq, int j) const;

  VectorField DtVelocity(int nq, const Vector &u) const;

  double DivVelocity(int nq, int j) const;

  double DivVelocity(int nq, const Vector &U) const;

  double Pressure(int nq, int j) const;

  double Pressure(int nq, const Vector &U) const;

  DampingVector DampingPressure(int nq, const Vector &u) const;

  double PressureLocal(const Point &localPoint, int j) const;

  double EvaluateTimeLocal(int q, int i) const;

  int GetTimeDimension() const;

  double PressureLocal(const Point &localPoint, const Vector &u) const;

  double PressureGlobal(const Point &globalPoint, int i) const;

  double PressureGlobal(const Point &globalPoint, const Vector &u) const;

  double DtPressure(int nq, int j) const;

  double DtPressure(int nq, const Vector &u) const;

  VectorField GradPressure(int nq, int j) const;

  VectorField GradPressure(int nq, const Vector &u) const;


  inline DampingVector DampingLocal(const Point &localPoint,
                                    const Vector &vector) const {
    DampingVector dampingVector(dampingComponentCount);
    const auto spaceDimension = TC.SpaceCell().dim();
    for (int i = 0; i < dampingComponentCount; ++i) {
      for (int j : component_to_index[spaceDimension + 1 + i]) {
        dampingVector[i] += vector(r, j) * PressureLocal(localPoint, j);
      }
    }
    return dampingVector;
  }

  inline DampingVector DampingGlobal(const Point &globalPoint,
                                     const Vector &vector) const {
    Point localPoint = TC.GlobalToLocal(globalPoint);
    return DampingLocal(localPoint, vector);
  }
};

using SpaceTimeViscoAcousticDGTElement = STDGDGViscoAcousticElement;

#endif  // SPACETIME_STDGDGVISCOACOUSTICELEMENT_HPP
