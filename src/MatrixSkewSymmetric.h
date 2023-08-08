#ifndef MATRIXSKEWSYMMETRIC_H
#define MATRIXSKEWSYMMETRIC_H

#pragma once

#include "MatrixSymmetry.h"

namespace pmat {

class MatrixSkewSymmetric : public pmat::MatrixSymmetry {

   public:
      MatrixSkewSymmetric() = default;
      explicit MatrixSkewSymmetric(const unsigned &size) : MatrixSymmetry::MatrixSymmetry(size){};
      MatrixSkewSymmetric(const MatrixSkewSymmetric &matrix)
          : MatrixSymmetry::MatrixSymmetry{std::move(matrix)} {}
      MatrixSkewSymmetric(MatrixSkewSymmetric &&matrix)
          : MatrixSymmetry::MatrixSymmetry{std::move(matrix)} {};
      MatrixSkewSymmetric &operator=(const MatrixSkewSymmetric &matrix) = default;
      MatrixSkewSymmetric &operator=(MatrixSkewSymmetric &&matrix) = default;
      ~MatrixSkewSymmetric() override = default;
      double operator()(const unsigned &row, const unsigned &column) const override;
      MatrixSkewSymmetric operator+(const MatrixSkewSymmetric &matrix) const;
      MatrixSquare operator+(const MatrixSymmetry &matrix) const;
      virtual void addBy(const MatrixSkewSymmetric &matrix);
      MatrixSkewSymmetric operator-(const MatrixSkewSymmetric &matrix) const;
      MatrixSquare operator-(const MatrixSymmetry &matrix) const;
      virtual void subtractBy(const MatrixSkewSymmetric &matrix);
      MatrixSkewSymmetric operator*(const double &scalar) const;
      MatrixSquare operator*(const MatrixSkewSymmetric &matrix) const;
      Vector operator*(const Vector &vector) const override;
      Matrix operator*(const Matrix &matrix) const override;
      void multiplyBy(const double &scalar) override;
      void transpose() override;
      void fillWithRandomValues(const double &min, const double &max) override;
};

} // namespace pmat
#endif