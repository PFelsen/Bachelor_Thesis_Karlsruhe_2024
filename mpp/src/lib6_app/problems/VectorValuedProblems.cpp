//
// Created by lstengel on 21.07.22.
//

#include "EllipticProblems.hpp"

IVectorValuedProblem *CreateVectorValuedProblem(const string &problemName) {
  if (problemName == "P0Test1D") return new VectorValuedPolynomialTest<0, INTERVAL>();
  else if (problemName == "P1Test1D") return new VectorValuedPolynomialTest<1, INTERVAL>();
  else if (problemName == "P2Test1D") return new VectorValuedPolynomialTest<2, INTERVAL>();
  else if (problemName == "P3Test1D") return new VectorValuedPolynomialTest<3, INTERVAL>();
  else if (problemName == "P4Test1D") return new VectorValuedPolynomialTest<4, INTERVAL>();
  else if (problemName == "P0Test2D") return new VectorValuedPolynomialTest<0, QUADRILATERAL>();
  else if (problemName == "P1Test2D") return new VectorValuedPolynomialTest<1, QUADRILATERAL>();
  else if (problemName == "P2Test2D") return new VectorValuedPolynomialTest<2, QUADRILATERAL>();
  else if (problemName == "P3Test2D") return new VectorValuedPolynomialTest<3, QUADRILATERAL>();
  else if (problemName == "P4Test2D") return new VectorValuedPolynomialTest<4, QUADRILATERAL>();
  else if (problemName == "P0Test3D") return new VectorValuedPolynomialTest<0, HEXAHEDRON>();
  else if (problemName == "P1Test3D") return new VectorValuedPolynomialTest<1, HEXAHEDRON>();
  else if (problemName == "P2Test3D") return new VectorValuedPolynomialTest<2, HEXAHEDRON>();
  else if (problemName == "P3Test3D") return new VectorValuedPolynomialTest<3, HEXAHEDRON>();
  else if (problemName == "P4Test3D") return new VectorValuedPolynomialTest<4, HEXAHEDRON>();
  else if (problemName == "P0Test2DTet") return new VectorValuedPolynomialTest<0, TRIANGLE>();
  else if (problemName == "P1Test2DTet") return new VectorValuedPolynomialTest<1, TRIANGLE>();
  else if (problemName == "P2Test2DTet") return new VectorValuedPolynomialTest<2, TRIANGLE>();
  else if (problemName == "P3Test2DTet") return new VectorValuedPolynomialTest<3, TRIANGLE>();
  else if (problemName == "P4Test2DTet") return new VectorValuedPolynomialTest<4, TRIANGLE>();
  else if (problemName == "P0Test3DTet") return new VectorValuedPolynomialTest<0, TETRAHEDRON>();
  else if (problemName == "P1Test3DTet") return new VectorValuedPolynomialTest<1, TETRAHEDRON>();
  else if (problemName == "P2Test3DTet") return new VectorValuedPolynomialTest<2, TETRAHEDRON>();
  else if (problemName == "P3Test3DTet") return new VectorValuedPolynomialTest<3, TETRAHEDRON>();
  else if (problemName == "P4Test3DTet") return new VectorValuedPolynomialTest<4, TETRAHEDRON>();
  else Exit(problemName + " not found");
}

std::unique_ptr<IVectorValuedProblem>
CreateVectorValuedProblemUnique(const std::string &problemName) {
  return std::unique_ptr<IVectorValuedProblem>(CreateVectorValuedProblem(problemName));
}

std::shared_ptr<IVectorValuedProblem>
CreateVectorValuedProblemShared(const std::string &problemName) {
  return std::shared_ptr<IVectorValuedProblem>(CreateVectorValuedProblem(problemName));
}