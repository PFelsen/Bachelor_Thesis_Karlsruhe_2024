#include "SymmetricCovariance.hpp"


RVector MOCE(RVector &toepRow, RVector &toepCol) {
  RVector circRow(toepRow.size() + toepCol.size() - 1);
  circRow.insert(toepRow, 0);
  std::vector<double> temp(toepCol.rbegin(), toepCol.rend() - 1);
  circRow.insert(RVector(std::move(temp)), toepRow.size());
  return circRow;
}

RMatrix BMOCE(RMatrix &toepRows, RMatrix &toepCols) {
  int circSize = toepRows.rows() + toepCols.rows() - 1;
  RMatrix circRows(circSize, circSize);
  for (int i = 0; i < toepRows.rows(); i++) {
    RVector toepRow = toepRows.row(i);
    RVector toepCol = toepCols.row(i);
    RVector circRow = MOCE(toepRow, toepCol);
    circRows.InsertRow(circRow, i);
  }
  for (int i = 0; i < toepRows.rows() - 1; i++) {
    RVector circRow = circRows.row(toepRows.rows() - 1 - i);
    circRows.InsertRow(circRow, toepRows.rows() + i);
  }
  return circRows;
}

RTensor BBMOCE(RTensor &toepRows, RTensor &toepCols) {
    int circSize = toepRows.ThirdDimension() + toepCols.ThirdDimension() - 1;
    RTensor circRows(circSize,circSize,circSize);
    for (int i = 0; i < toepRows.ThirdDimension(); i++) {
        RMatrix toepRow = toepRows.ThirdComponent(i);
        RMatrix toepCol = toepCols.ThirdComponent(i);
        RMatrix circRow = BMOCE(toepRow, toepCol);
        circRows.InsertMatrixThirdDimension(circRow, i, 0, 0);
    }
    for (int i = 0; i < toepRows.ThirdDimension() - 1; i++) {
        RMatrix circRow = circRows.ThirdComponent(toepRows.ThirdDimension() - 1 - i);
        circRows.InsertMatrixThirdDimension(circRow, toepRows.ThirdDimension() + i);
    }
    return circRows;
}

std::vector<double> linspace(const double &start, const double &end, int num) {
  std::vector<double> linspaced;
  if (num == 0) { return linspaced; }
  if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }
  double delta = (end - start) / (num - 1);
  for (int i = 0; i < num - 1; i++) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end);
  return linspaced;
}

void CovarianceFunction1D::ToeplitzMatrix(CovarianceFunction1D::T &toepRow,
                                          CovarianceFunction1D::T &toepCol) {
  RVector firstCoord(linspace(0, 1, toepRow.size()));
  double y_rows = firstCoord.front();

  for (int index_c = 0; index_c < toepRow.size(); index_c++) {
    double c1 = firstCoord[index_c];
    double x_rows = c1;
    double tau_rows = x_rows - y_rows;

    double x_columns = firstCoord.front();
    double y_columns = c1;
    double tau_columns = x_columns - y_columns;

    toepRow[index_c] = Function(&tau_rows);
    toepCol[index_c] = Function(&tau_columns);
  }
}

void CovarianceFunction2D::ToeplitzMatrix(CovarianceFunction2D::T &toepRows,
                                          CovarianceFunction2D::T &toepCols) {
  RVector coordinate1(linspace(0, 1, toepRows.rows()));
  RVector coordinate2(linspace(0, 1, toepCols.rows()));

  double y_rows[2] = {coordinate1.front(), coordinate2.front()};
  for (int index_c2 = 0; index_c2 < toepRows.rows(); index_c2++) {
    for (int index_c1 = 0; index_c1 < toepRows.rows(); index_c1++) {
      double c1 = coordinate1[index_c1];
      double c2 = coordinate2[index_c2];
      double x_rows[2] = {c1, c2};
      double tau_rows[2] = {x_rows[0] - y_rows[0], x_rows[1] - y_rows[1]};

      double x_columns[2] = {coordinate1.front(), c2};
      double y_columns[2] = {c1, coordinate2.front()};
      double tau_columns[2] = {x_columns[0] - y_columns[0],
                               x_columns[1] - y_columns[1]};

      toepCols[index_c2][index_c1] = Function(tau_columns);
      toepRows[index_c2][index_c1] = Function(tau_rows);
    }
  }
}

void CovarianceFunction3D::ToeplitzMatrix(CovarianceFunction3D::T &toepRows,
                                          CovarianceFunction3D::T &toepCols) {
    RVector coordinate1(linspace(0, 1, toepRows.FirstDimension()));
    RVector coordinate2(linspace(0, 1, toepCols.FirstDimension()));
    RVector coordinate3(linspace(0, 1, toepCols.FirstDimension()));

    double y_rows[3] = {coordinate1.front(), coordinate2.front(), coordinate3.front()};
    for (int index_c3 = 0; index_c3 < toepRows.FirstDimension(); index_c3++) {
        for (int index_c2 = 0; index_c2 < toepRows.FirstDimension(); index_c2++) {
            for (int index_c1 = 0; index_c1 < toepRows.FirstDimension(); index_c1++) {
                double c1 = coordinate1[index_c1];
                double c2 = coordinate2[index_c2];
                double c3 = coordinate3[index_c3];
                double x_rows[3] = {c1, c2, c3};
                double tau_rows[3] = {x_rows[0] - y_rows[0], x_rows[1] - y_rows[1], x_rows[2] - y_rows[2]};

                //Todo: Check what columns need to be. either c3 or coordinate3.front()
                double x_columns[3] = {coordinate1.front(), c2, coordinate3.front()};
                double y_columns[3] = {c1, coordinate2.front(), c3};
                double tau_columns[3] = {x_columns[0] - y_columns[0],
                                         x_columns[1] - y_columns[1],
                                         x_columns[2] - y_columns[2]};

                //Todo: Check indices for toepCols and Rows

                toepCols(index_c3, index_c2, index_c1) = Function(tau_columns);
                toepRows(index_c3, index_c2, index_c1) = Function(tau_rows);
            }
        }
    }
}
