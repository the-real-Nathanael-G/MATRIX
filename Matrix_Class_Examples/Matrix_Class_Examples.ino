#include <Matrix.h>

void setup() {
  Serial.begin(115200);
  
  //3x3 Matrix filled with given data
  float data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Matrix Example1(3, 3, data);
  Example1.print();
  Serial.println("");
  /*
  Example1 Data:
    | 1 , 2 , 3 |
    | 4 , 5 , 6 |
    | 7 , 8 , 9 |
  */

  //Empty 2x2 Matrix
  Matrix Example2(2, 2);
  //Fills Example2 with number '3' in every position
  Example2.fill(3);
  Example2.print();
  Serial.println("");
  //Assign each element of the Example2 matrix individually
  Example2(0, 0) = 6;
  Example2(0, 1) = 23;
  Example2(1, 0) = 18;
  Example2(1, 1) = 2;
  /*
  Example2 Data:
    | 6 , 23 |
    | 18,  2 |
  */
  
  //MATRIX CLASS TESTBENCHES
  test_joins();
  test_matrix_multiplier();
  test_copyRow_safe();
  test_determinant_inverse();

}

void loop() {
  // put your main code here, to run repeatedly:

}

void test_joins() {
  Serial.println("\n[TEST 3b] joinHorizontal / joinVertical");

  // A (2x2)
  float Av[] = {1,2, 3,4};
  Matrix A(2,2, Av);

  // B (2x1) for horizontal join
  float Bv[] = {5, 6};
  Matrix B(2,1, Bv);

  // C (1x2) for vertical join needs same cols as result after H-join test; but first test H-join:
  Serial.println("A:");
  A.print();
  Serial.println("B:");
  B.print();

  Matrix AH = A; // copy
  AH.joinHorizontal(B);
  // Expected:
  // [1 2 5]
  // [3 4 6]
  Serial.println("A | B (expected rows: [1 2 5], [3 4 6])");
  AH.print();

  // For vertical join, use same cols: A (2x2) and D (1x2)
  float Dv[] = {7, 8};
  Matrix D(1,2, Dv);
  Matrix AV = A; // copy
  AV.joinVertical(D);
  // Expected:
  // [1 2]
  // [3 4]
  // [7 8]
  Serial.println("A over D (expected rows: [1 2], [3 4], [7 8])");
  AV.print();
}

void test_matrix_multiplier() {
  Serial.println("\n[TEST 3a] fast matmul");

  // A: 3x2
  float Av[] = {1,2, 3,4, 5,6};
  Matrix A(3,2, Av);
  // B: 2x3
  float Bv[] = {7,8,9, 10,11,12};
  Matrix B(2,3, Bv);

  // Expected C = A*B:
  // [1*7+2*10, 1*8+2*11, 1*9+2*12] = [27, 30, 33]
  // [3*7+4*10, 3*8+4*11, 3*9+4*12] = [61, 68, 75]
  // [5*7+6*10, 5*8+6*11, 5*9+6*12] = [95,106,117]
  Matrix C = A * B;
  Serial.println("C = A*B (expected rows: [27 30 33], [61 68 75], [95 106 117])");
  C.print();
}

void test_copyRow_safe() {
  Serial.println("\n[TEST 2] getRow / copyRowInto");

  float vals[] = {1,2,3, 4,5,6, 7,8,9};
  Matrix M(3,3, vals);

  Matrix r1 = M.getRow(1); // expected [4 5 6]
  Serial.println("Row 1 (Matrix 1x3):");
  r1.print();

  float buf[3] = {0,0,0};
  bool ok = M.copyRowInto(2, buf, 3);
  Serial.print("copyRowInto ok? ");
  Serial.println(ok ? "true":"false");
  Serial.print("buf: "); Serial.print(buf[0]);
  Serial.print(" ");
  Serial.print(buf[1]);
  Serial.print(" ");
  Serial.println(buf[2]); // expected [7 8 9]
}

void test_determinant_inverse() {
  Serial.println("\n[TEST 1] determinant & inverse");

  // 2x2
  float a2v[] = {4, 7, 2, 6};
  Matrix A2(2,2,a2v);
  float det2 = A2.determinant();
  Serial.print("det(A2) expected 10, got ");
  Serial.println(det2, 6);
  Matrix I2 = A2 * A2.inverse();
  Serial.println("A2 * inv(A2) ~ I:");
  I2.print();

  // 3x3
  float a3v[] = {
    1, 2, 3,
    0, 1, 4,
    5, 6, 0
  };
  Matrix A3(3,3,a3v);
  float det3 = A3.determinant();   // expected 1*(1*0-4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5) = -24 + 40 - 15 = 1
  Serial.print("det(A3) expected 1, got "); Serial.println(det3, 6);
  Matrix I3 = A3 * A3.inverse();
  Serial.println("A3 * inv(A3) ~ I:");
  I3.print();

  // 4x4 (general path)
  Matrix T(4,4); T.fill(0);
  // make a well-conditioned upper-triangular matrix with diag [2,3,4,5]
  T(0,0)=2; T(0,1)=1; T(0,2)=0; T(0,3)=0;
  T(1,1)=3; T(1,2)=1; T(1,3)=0;
  T(2,2)=4; T(2,3)=1;
  T(3,3)=5;
  float det4 = T.determinant();    // product of diag = 2*3*4*5 = 120
  Serial.print("det(T) expected 120, got "); Serial.println(det4, 6);
  Matrix I4 = T * T.inverse();
  Serial.println("T * inv(T) ~ I:");
  I4.print();
}
