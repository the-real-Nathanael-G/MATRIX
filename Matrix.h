/*
Nominal numbers (i.e. NxM) are index 1
Ordinal numbers (i.e. when retrieving Matrix(row, col)) are index 0
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <Arduino.h>
#include <math.h>

class Matrix {
private:
  int rows;
  int cols;
  float* data;

public:
  //Constructor: Creates object but no data is stored
  Matrix(){
    this->rows = 0;
    this->cols = 0;
    data = nullptr;
  }

  //Constructor: Allocates memory for a rows x cols matrix. Matrix is empty
  Matrix(int rows, int cols) : rows(rows), cols(cols) {
    data = new float[rows * cols];
  }

  //Constructor: Allocates memory and fils with initial values
  Matrix(int rows, int cols, const float* initVals) : rows(rows), cols(cols) {
    data = new float[rows * cols];
    for (int i = 0; i < rows * cols; i++) {
      data[i] = initVals[i];
    }
  }

  //Copy constructor
  Matrix(const Matrix &other) : rows(other.rows), cols(other.cols) {
    data = new float[rows * cols];
    for (int i = 0; i < rows * cols; i++) {
      data[i] = other.data[i];
    }
  }

  //Overloaded '=' assignment operator
  Matrix& operator=(const Matrix& other){
    if (this == &other){
      return *this;
    }
    float* newData = new (std::nothrow) float[other.rows * other.cols];
    //if (!newData) { /* handle error (return *this or flag) */ }
    for (int i=0; i<other.rows*other.cols; ++i){
      newData[i] = other.data[i];
    }
    delete[] data;
    data = newData; rows = other.rows; cols = other.cols;
    return *this;
  }

  //Overloaded multiplication operator for matrix multiplication
  Matrix operator*(const Matrix &other) const {
    if (cols != other.rows) {
      Serial.println("Error: Matrix dimensions are incompatible for multiplication.");
      return Matrix(0, 0);
    }
    Matrix result(rows, other.cols);

    const int Arows = rows;
    const int Acols = cols;          // == other.rows
    const int Bcols = other.cols;

    for (int i = 0; i < Arows; ++i) {
      const int arow = i * Acols;
      const int rrow = i * Bcols;
      for (int j = 0; j < Bcols; ++j) {
        float sum = 0.0f;
        // unrolled inner indexing: other.data[k*Bcols + j]
        for (int k = 0; k < Acols; ++k) {
          sum += data[arow + k] * other.data[k * Bcols + j];
        }
        result.data[rrow + j] = sum;
      }
    }
    return result;
  }

  //Overloaded multiplication operator for scalar multiplication
  Matrix operator*(float scalar) const {
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; i++) {
      result.data[i] = data[i] * scalar;
    }
    return result;
  }

  //Overloaded addition operator for matrix addition
  Matrix operator+(const Matrix &other) const {
    if (rows != other.rows || cols != other.cols) {
      Serial.println("Error: Matrices must have the same dimensions for addition.");
      return Matrix(0, 0); // Return an empty matrix in case of error
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        result(i, j) = (*this)(i, j) + other(i, j);
      }
    }
    return result;
  }

  //Overloaded subtraction operator for matrix subtraction
  Matrix operator-(const Matrix &other) const {
    if (rows != other.rows || cols != other.cols) {
      Serial.println("Error: Matrices must have the same dimensions for subtraction.");
      return Matrix(0, 0); // Return an empty matrix in case of error
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        result(i, j) = (*this)(i, j) - other(i, j);
      }
    }
    return result;
  }

  //Access operator (non-const) to access matrix element at (r, c)
  float& operator()(int r, int c) {
    return data[r * cols + c];
  }

  //Access operator (const) for read-only access
  const float& operator()(int r, int c) const {
    return data[r * cols + c];
  }

  //Function to set all elements to a given value
  void fill(float value) {
    for (int i = 0; i < rows * cols; i++) {
      data[i] = value;
    }
  }

  //Function to horizontally join this matrix to another (appending columns)
  void joinHorizontal(const Matrix &other) {
    if (rows != other.rows) {
      Serial.println("Error: joinHorizontal - Matrices must have the same number of rows.");
      return;
    }
    const int newCols = cols + other.cols;
    float* newData = new float[rows * newCols];
    if (!newData) {
      Serial.println("Error: Out of memory in joinHorizontal.");
      return;
    }

    for (int i = 0; i < rows; ++i) {
      // left block (this row)
      memcpy(newData + i * newCols, data + i * cols, cols * sizeof(float));
      // right block (other row)
      memcpy(newData + i * newCols + cols, other.data + i * other.cols, other.cols * sizeof(float));
    }

    delete[] data;
    data = newData;
    cols = newCols;
  }

  //Function to vertically join this matrix to another (appending rows)
  void joinVertical(const Matrix &other) {
    if (cols != other.cols) {
      Serial.println("Error: joinVertical - Matrices must have the same number of columns.");
      return;
    }
    const int newRows = rows + other.rows;
    float* newData = new float[newRows * cols];
    if (!newData) {
      Serial.println("Error: Out of memory in joinVertical.");
      return;
    }

    // top block: this
    memcpy(newData, data, rows * cols * sizeof(float));
    // bottom block: other
    memcpy(newData + rows * cols, other.data, other.rows * cols * sizeof(float));

    delete[] data;
    data = newData;
    rows = newRows;
  }

  //Return a 1xcols Matrix of required row
  Matrix getRow(int rowIndex) const {
    if (rowIndex < 0 || rowIndex >= rows) {
      Serial.println("Error: Row index out of range.");
      return Matrix(0, 0);
    }
    Matrix out(1, cols);
    for (int j = 0; j < cols; ++j) out(0, j) = data[rowIndex * cols + j];
    return out;
  }

  //Copy required row into buffer provided. Return success/failure
  bool copyRowInto(int rowIndex, float* buf, int buf_len) const {
    if (rowIndex < 0 || rowIndex >= rows) {
      Serial.println("Error: Row index out of range.");
      return false;
    }
    if (buf_len < cols) {
      Serial.println("Error: Destination buffer too small.");
      return false;
    }
    // contiguous copy
    memcpy(buf, &data[rowIndex * cols], cols * sizeof(float));
    return true;
  }

  //When we initialise an empty Matrix, we need to give it a size before assigning values to it
  void reSize(int newRows, int newCols) {
    delete[] data; //Just to be sure the Matrix is empty
    rows = newRows;
    cols = newCols;
    data = new float[rows * cols];
  }

  //Determinant: 1x1/2x2/3x3 fast paths, otherwise Gaussian elimination with partial pivot
  float determinant() const {
    if (rows != cols) {
      Serial.println("Error: Determinant undefined for non-square matrices.");
      return 0.0f;
    }
    const int n = rows;
    if (n == 0) return 0.0f;

    // 1x1
    if (n == 1) return data[0];

    // 2x2
    if (n == 2) {
      const float a = data[0], b = data[1];
      const float c = data[2], d = data[3];
      return a*d - b*c;
    }

    // 3x3 (Sarrus)
    if (n == 3) {
      const float a = data[0], b = data[1], c = data[2];
      const float d = data[3], e = data[4], f = data[5];
      const float g = data[6], h = data[7], i = data[8];
      return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    }

    // n >= 4: Gaussian elimination with partial pivoting
    float* a = new float[n * n];
    if (!a) {
      Serial.println("Error: Out of memory in determinant().");
      return 0.0f;
    }
    memcpy(a, data, n * n * sizeof(float));

    int sign = 1;
    float det = 1.0f;
    for (int i = 0; i < n; ++i) {
      // pivot search
      int pivot = i;
      float maxAbs = fabsf(a[i*n + i]);
      for (int r = i + 1; r < n; ++r) {
        float v = fabsf(a[r*n + i]);
        if (v > maxAbs) { maxAbs = v; pivot = r; }
      }
      if (maxAbs < 1e-12f) { det = 0.0f; break; }

      // row swap if needed
      if (pivot != i) {
        for (int c = i; c < n; ++c) {
          float tmp = a[i*n + c];
          a[i*n + c] = a[pivot*n + c];
          a[pivot*n + c] = tmp;
        }
        sign = -sign;
      }

      const float piv = a[i*n + i];
      det *= piv;

      // eliminate below
      for (int r = i + 1; r < n; ++r) {
        float f = a[r*n + i] / piv;
        // a[r*n + i] becomes zero
        for (int c = i + 1; c < n; ++c) {
          a[r*n + c] -= f * a[i*n + c];
        }
      }
    }

    delete[] a;
    return det * sign;
  }

  //Inverse: 1x1/2x2/3x3 fast paths, otherwise Gauss-Jordan on [A|I] with pivoting
  Matrix inverse() const {
    if (rows != cols) {
      Serial.println("Error: Inverse undefined for non-square matrices.");
      return Matrix(0, 0);
    }
    const int n = rows;
    if (n == 0) return Matrix(0, 0);

    // 1x1
    if (n == 1) {
      if (fabsf(data[0]) < 1e-12f) {
        Serial.println("Error: Matrix is singular (1x1).");
        return Matrix(0, 0);
      }
      Matrix inv(1,1);
      inv(0,0) = 1.0f / data[0];
      return inv;
    }

    // 2x2 closed form
    if (n == 2) {
      const float a = data[0], b = data[1];
      const float c = data[2], d = data[3];
      const float det = a*d - b*c;
      if (fabsf(det) < 1e-12f) {
        Serial.println("Error: Matrix is singular (2x2).");
        return Matrix(0, 0);
      }
      Matrix inv(2,2);
      inv(0,0) =  d / det; inv(0,1) = -b / det;
      inv(1,0) = -c / det; inv(1,1) =  a / det;
      return inv;
    }

    // 3x3 closed form via adjugate / det
    if (n == 3) {
      const float a = data[0], b = data[1], c = data[2];
      const float d = data[3], e = data[4], f = data[5];
      const float g = data[6], h = data[7], i = data[8];

      const float A =  (e*i - f*h);
      const float B = -(d*i - f*g);
      const float C =  (d*h - e*g);
      const float D = -(b*i - c*h);
      const float E =  (a*i - c*g);
      const float F = -(a*h - b*g);
      const float G =  (b*f - c*e);
      const float H = -(a*f - c*d);
      const float I =  (a*e - b*d);

      const float det = a*A + b*B + c*C;
      if (fabsf(det) < 1e-12f) {
        Serial.println("Error: Matrix is singular (3x3).");
        return Matrix(0, 0);
      }

      Matrix inv(3,3);
      // adjugate (transpose of cofactor matrix), then / det
      inv(0,0) = A / det; inv(0,1) = D / det; inv(0,2) = G / det;
      inv(1,0) = B / det; inv(1,1) = E / det; inv(1,2) = H / det;
      inv(2,0) = C / det; inv(2,1) = F / det; inv(2,2) = I / det;
      return inv;
    }

    // n >= 4: Gauss-Jordan on augmented [A|I] with partial pivoting
    const int W = 2 * n;
    float* aug = new float[n * W];
    if (!aug) {
      Serial.println("Error: Out of memory in inverse().");
      return Matrix(0, 0);
    }

    // build [A | I]
    for (int r = 0; r < n; ++r) {
      for (int c = 0; c < n; ++c)      aug[r*W + c]     = data[r*cols + c];
      for (int c = 0; c < n; ++c)      aug[r*W + n + c] = (r == c) ? 1.0f : 0.0f;
    }

    for (int i = 0; i < n; ++i) {
      // pivot
      int pivot = i;
      float maxAbs = fabsf(aug[i*W + i]);
      for (int r = i + 1; r < n; ++r) {
        float v = fabsf(aug[r*W + i]);
        if (v > maxAbs) { maxAbs = v; pivot = r; }
      }
      if (maxAbs < 1e-12f) {
        Serial.println("Error: Matrix is singular (Gauss-Jordan).");
        delete[] aug;
        return Matrix(0, 0);
      }
      if (pivot != i) {
        for (int c = 0; c < W; ++c) {
          float tmp = aug[i*W + c];
          aug[i*W + c] = aug[pivot*W + c];
          aug[pivot*W + c] = tmp;
        }
      }

      // normalize pivot row
      const float piv = aug[i*W + i];
      const float invPiv = 1.0f / piv;
      for (int c = 0; c < W; ++c) aug[i*W + c] *= invPiv;

      // eliminate other rows
      for (int r = 0; r < n; ++r) {
        if (r == i) continue;
        const float f = aug[r*W + i];
        if (f != 0.0f) {
          for (int c = 0; c < W; ++c) {
            aug[r*W + c] -= f * aug[i*W + c];
          }
        }
      }
    }

    Matrix inv(n, n);
    for (int r = 0; r < n; ++r)
      for (int c = 0; c < n; ++c)
        inv(r, c) = aug[r*W + (n + c)];

    delete[] aug;
    return inv;
  }

  //Read-only accessor for the number of rows
  int getRows() const {
    return rows;
  }

  //Read-only accessor for the number of columns
  int getCols() const {
    return cols;
  }

  //Read-only accessor for the data pointer
  const float* getData() const {
    return data;
  }

  // Function to print the matrix to Serial
  void print() const {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        Serial.print((*this)(i, j));
        Serial.print("\t");
      }
      Serial.println();
    }
  }

  //Function to delete a row at the specified index from the matrix
  void deleteRow(int rowIndex) {
    if (rowIndex < 0 || rowIndex >= rows) {
      Serial.println("Error: Row index out of range for deletion.");
      return;
    }
    int newRows = rows - 1;
    //Allocate new array with one fewer row
    float* newData = new float[newRows * cols];

    //Copy rows before the deleted row
    for (int i = 0; i < rowIndex; i++) {
      for (int j = 0; j < cols; j++) {
        newData[i * cols + j] = (*this)(i, j);
      }
    }
    //Copy rows after the deleted row
    for (int i = rowIndex + 1; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        newData[(i - 1) * cols + j] = (*this)(i, j);
      }
    }

    //Delete old data and update matrix dimensions
    delete[] data;
    data = newData;
    rows = newRows;
  }

  //Function to delete a column at the specified index from the matrix
  void deleteColumn(int colIndex) {
    if (colIndex < 0 || colIndex >= cols) {
      Serial.println("Error: Column index out of range for deletion.");
      return;
    }
    int newCols = cols - 1;
    float* newData = new float[rows * newCols];
    
    //For each row, copy over every column except the one to delete
    for (int i = 0; i < rows; i++) {
      int newCol = 0; //New column index
      for (int j = 0; j < cols; j++) {
        if (j == colIndex) continue;  //Skip the column to delete
        newData[i * newCols + newCol] = (*this)(i, j);
        newCol++;
      }
    }
    
    //Delete old data and update matrix dimensions
    delete[] data;
    data = newData;
    cols = newCols;
  }

  //Destructor
  ~Matrix() {
    delete[] data;
  }
};

#endif
