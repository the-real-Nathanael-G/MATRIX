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
  Matrix();

  //Constructor: Allocates memory for a rows x cols matrix. Matrix is empty
  Matrix(int rows, int cols);

  //Constructor: Allocates memory and fils with initial values
  Matrix(int rows, int cols, const float* initVals);

  //Copy constructor
  Matrix(const Matrix &other);

  //Overloaded '=' assignment operator
  Matrix& operator=(const Matrix& other);

  //Overloaded multiplication operator for matrix multiplication
  Matrix operator*(const Matrix &other) const;

  //Overloaded multiplication operator for scalar multiplication
  Matrix operator*(float scalar) const;

  //Overloaded addition operator for matrix addition
  Matrix operator+(const Matrix &other) const;

  //Overloaded subtraction operator for matrix subtraction
  Matrix operator-(const Matrix &other) const;

  //Access operator (non-const) to access matrix element at (r, c)
  float& operator()(int r, int c);

  //Access operator (const) for read-only access
  const float& operator()(int r, int c) const;

  //Function to set all elements to a given value
  void fill(float value);

  //Function to horizontally join this matrix to another (appending columns)
  void joinHorizontal(const Matrix &other);

  //Function to vertically join this matrix to another (appending rows)
  void joinVertical(const Matrix &other);

  //Return a 1xcols Matrix of required row
  Matrix getRow(int rowIndex) const;

  //Copy required row into buffer provided. Return success/failure
  bool copyRowInto(int rowIndex, float* buf, int buf_len) const;

  //When we initialise an empty Matrix, we need to give it a size before assigning values to it
  void reSize(int newRows, int newCols);

  //Determinant: 1x1/2x2/3x3 fast paths, otherwise Gaussian elimination with partial pivot
  float determinant() const;

  //Inverse: 1x1/2x2/3x3 fast paths, otherwise Gauss-Jordan on [A|I] with pivoting
  Matrix inverse() const;

  //Read-only accessor for the number of rows
  int getRows() const;

  //Read-only accessor for the number of columns
  int getCols() const;

  //Read-only accessor for the data pointer
  const float* getData() const;

  // Function to print the matrix to Serial
  void print() const;

  //Function to delete a row at the specified index from the matrix
  void deleteRow(int rowIndex);

  //Function to delete a column at the specified index from the matrix
  void deleteColumn(int colIndex);

  //Destructor
  ~Matrix();
};

#endif
