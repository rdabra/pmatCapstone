#include "../src/DecompositionPLU.h"
#include "../src/MatrixSkewSymmetric.h"
#include "../src/MatrixSymmetric.h"
#include "../src/utils.h"
#include "gtest/gtest.h"

using namespace pmat;

TEST(TestMatrixSkewSymmetric, TestEqualityOperator) {
   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   MatrixSkewSymmetric v(4);
   v.setValue(1.0, 0, 0);
   v.setValue(4.0, 1, 0);
   v.setValue(5.0, 1, 1);
   v.setValue(7.0, 2, 0);
   v.setValue(8.0, 2, 1);
   v.setValue(9.0, 2, 2);
   v.setValue(10.0, 3, 0);
   v.setValue(11.0, 3, 1);
   v.setValue(12.0, 3, 2);
   v.setValue(13.0, 3, 3);

   EXPECT_TRUE(v == z);
}

TEST(TestMatrixSkewSymmetric, TestDotProduct) {
   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   MatrixSkewSymmetric v(4);
   v.setValue(1.5, 0, 0);
   v.setValue(4.5, 1, 0);
   v.setValue(5.5, 1, 1);
   v.setValue(7.5, 2, 0);
   v.setValue(8.5, 2, 1);
   v.setValue(9.5, 2, 2);
   v.setValue(10.5, 3, 0);
   v.setValue(11.5, 3, 1);
   v.setValue(12.5, 3, 2);
   v.setValue(13.5, 3, 3);

   double resp = z.dotProduct(v);

   EXPECT_TRUE(utils::areEqual(resp, 1330.0));
}

TEST(TestMatrixSkewSymmetric, TestPlus) {
   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   MatrixSkewSymmetric v(4);
   v.setValue(1.0, 0, 0);
   v.setValue(4.0, 1, 0);
   v.setValue(5.0, 1, 1);
   v.setValue(7.0, 2, 0);
   v.setValue(8.0, 2, 1);
   v.setValue(9.0, 2, 2);
   v.setValue(10.0, 3, 0);
   v.setValue(11.0, 3, 1);
   v.setValue(12.0, 3, 2);
   v.setValue(13.0, 3, 3);

   MatrixSkewSymmetric resp(4);
   resp.setValue(2.0, 0, 0);
   resp.setValue(8.0, 1, 0);
   resp.setValue(10.0, 1, 1);
   resp.setValue(14.0, 2, 0);
   resp.setValue(16.0, 2, 1);
   resp.setValue(18.0, 2, 2);
   resp.setValue(20.0, 3, 0);
   resp.setValue(22.0, 3, 1);
   resp.setValue(24.0, 3, 2);
   resp.setValue(26.0, 3, 3);

   MatrixSquare resp1(4);
   resp1.setValue(2.0, 0, 0);
   resp1.setValue(8.0, 1, 0);
   resp1.setValue(-8.0, 0, 1);
   resp1.setValue(10.0, 1, 1);
   resp1.setValue(14.0, 2, 0);
   resp1.setValue(-14.0, 0, 2);
   resp1.setValue(16.0, 2, 1);
   resp1.setValue(-16.0, 1, 2);
   resp1.setValue(18.0, 2, 2);
   resp1.setValue(20.0, 3, 0);
   resp1.setValue(-20.0, 0, 3);
   resp1.setValue(22.0, 3, 1);
   resp1.setValue(-22.0, 1, 3);
   resp1.setValue(24.0, 3, 2);
   resp1.setValue(-24.0, 2, 3);
   resp1.setValue(26.0, 3, 3);

   MatrixSkewSymmetric r1(4);
   r1 = z + v;
   MatrixSkewSymmetric r2(z + v);
   z.addBy(v);

   MatrixSymmetric s{4};
   MatrixSquare rr{s + z};

   EXPECT_TRUE(resp == r1);
   EXPECT_TRUE(resp1 == r2);
   EXPECT_TRUE(resp == z);
   EXPECT_TRUE(rr == z);
}

TEST(TestMatrixSkewSymmetric, TestMinus) {
   MatrixSkewSymmetric z(4);
   z.setValue(2.0, 0, 0);
   z.setValue(8.0, 1, 0);
   z.setValue(10.0, 1, 1);
   z.setValue(14.0, 2, 0);
   z.setValue(16.0, 2, 1);
   z.setValue(18.0, 2, 2);
   z.setValue(20.0, 3, 0);
   z.setValue(22.0, 3, 1);
   z.setValue(24.0, 3, 2);
   z.setValue(26.0, 3, 3);

   MatrixSkewSymmetric v(4);
   v.setValue(1.0, 0, 0);
   v.setValue(4.0, 1, 0);
   v.setValue(5.0, 1, 1);
   v.setValue(7.0, 2, 0);
   v.setValue(8.0, 2, 1);
   v.setValue(9.0, 2, 2);
   v.setValue(10.0, 3, 0);
   v.setValue(11.0, 3, 1);
   v.setValue(12.0, 3, 2);
   v.setValue(13.0, 3, 3);

   MatrixSkewSymmetric resp(4);
   resp.setValue(1.0, 0, 0);
   resp.setValue(4.0, 1, 0);
   resp.setValue(5.0, 1, 1);
   resp.setValue(7.0, 2, 0);
   resp.setValue(8.0, 2, 1);
   resp.setValue(9.0, 2, 2);
   resp.setValue(10.0, 3, 0);
   resp.setValue(11.0, 3, 1);
   resp.setValue(12.0, 3, 2);
   resp.setValue(13.0, 3, 3);

   MatrixSquare resp1(4);
   resp1.setValue(1.0, 0, 0);
   resp1.setValue(4.0, 1, 0);
   resp1.setValue(-4.0, 0, 1);
   resp1.setValue(5.0, 1, 1);
   resp1.setValue(7.0, 2, 0);
   resp1.setValue(-7.0, 0, 2);
   resp1.setValue(8.0, 2, 1);
   resp1.setValue(-8.0, 1, 2);
   resp1.setValue(9.0, 2, 2);
   resp1.setValue(10.0, 3, 0);
   resp1.setValue(-10.0, 0, 3);
   resp1.setValue(11.0, 3, 1);
   resp1.setValue(-11.0, 1, 3);
   resp1.setValue(12.0, 3, 2);
   resp1.setValue(-12.0, 2, 3);
   resp1.setValue(13.0, 3, 3);

   MatrixSkewSymmetric r1(4);
   r1 = z - v;
   MatrixSkewSymmetric r2(z - v);
   z.subtractBy(v);

   MatrixSymmetric s{4};
   MatrixSquare rr{z - s};

   EXPECT_TRUE(resp == r1);
   EXPECT_TRUE(resp1 == r2);
   EXPECT_TRUE(resp == z);
   EXPECT_TRUE(rr == z);
}

TEST(TestMatrixSkewSymmetric, TestTimes) {

   MatrixSkewSymmetric z1(2);
   z1.setValue(1.0, 0, 0);
   z1.setValue(2.0, 1, 0);
   z1.setValue(3.0, 1, 1);

   MatrixSkewSymmetric v1(2);
   z1.setValue(1.0, 0, 0);
   z1.setValue(2.0, 1, 0);
   z1.setValue(3.0, 1, 1);

   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   MatrixSkewSymmetric v(4);
   v.setValue(2.0, 0, 0);
   v.setValue(5.0, 1, 0);
   v.setValue(6.0, 1, 1);
   v.setValue(8.0, 2, 0);
   v.setValue(9.0, 2, 1);
   v.setValue(10.0, 2, 2);
   v.setValue(11.0, 3, 0);
   v.setValue(12.0, 3, 1);
   v.setValue(13.0, 3, 2);
   v.setValue(14.0, 3, 3);

   MatrixSquare resp(4);

   resp.setValue(-184.0, 0, 0);
   resp.setValue(-212.0, 0, 1);
   resp.setValue(-172.0, 0, 2);
   resp.setValue(-12.0, 0, 3);

   resp.setValue(-152.0, 1, 0);
   resp.setValue(-194.0, 1, 1);
   resp.setValue(-300.0, 1, 2);
   resp.setValue(-154.0, 1, 3);

   resp.setValue(-6.0, 2, 0);
   resp.setValue(-50.0, 2, 1);
   resp.setValue(-194.0, 2, 2);
   resp.setValue(-458.0, 2, 3);

   resp.setValue(314.0, 3, 0);
   resp.setValue(280.0, 3, 1);
   resp.setValue(110.0, 3, 2);
   resp.setValue(-216.0, 3, 3);

   MatrixSkewSymmetric resp1(4);
   resp1.setValue(2.0, 0, 0);
   resp1.setValue(8.0, 1, 0);
   resp1.setValue(10.0, 1, 1);
   resp1.setValue(14.0, 2, 0);
   resp1.setValue(16.0, 2, 1);
   resp1.setValue(18.0, 2, 2);
   resp1.setValue(20.0, 3, 0);
   resp1.setValue(22.0, 3, 1);
   resp1.setValue(24.0, 3, 2);
   resp1.setValue(26.0, 3, 3);

   MatrixSquare r1(4);
   r1 = z * v;
   MatrixSquare r2(z * v);
   MatrixSkewSymmetric r3;
   r3 = z;

   r3.multiplyBy(2.0);

   EXPECT_TRUE(resp == r1);
   EXPECT_TRUE(resp == r2);
   EXPECT_TRUE(z * 2.0 == resp1);
   EXPECT_TRUE(r3 == resp1);
}

TEST(TestMatrixSkewSymmetric, TestFrobenius) {
   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   EXPECT_TRUE(utils::areEqual(z.getFrobeniusNorm(), 35.55277767));
}

TEST(TestMatrixSkewSymmetric, TestDeterminant) {
   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   EXPECT_TRUE(utils::areEqual(z.decomposeToPLU().determinant(), 15384.0000000000));
}

TEST(TestMatrixSkewSymmetric, TestTranspose) {
   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   MatrixSkewSymmetric v(z * -1.0);

   z.transpose();

   EXPECT_TRUE(v == z);
}

TEST(TestMatrixSkewSymmetric, TestMisc) {

   MatrixSkewSymmetric z(4);
   z.setValue(1.0, 0, 0);
   z.setValue(4.0, 1, 0);
   z.setValue(5.0, 1, 1);
   z.setValue(7.0, 2, 0);
   z.setValue(8.0, 2, 1);
   z.setValue(9.0, 2, 2);
   z.setValue(10.0, 3, 0);
   z.setValue(11.0, 3, 1);
   z.setValue(12.0, 3, 2);
   z.setValue(13.0, 3, 3);

   MatrixSkewSymmetric v(4);
   v.setValue(1.0, 0, 0);
   v.setValue(4.0, 1, 0);
   v.setValue(5.0, 1, 1);
   v.setValue(7.0, 2, 0);
   v.setValue(8.0, 2, 1);
   v.setValue(9.0, 2, 2);
   v.setValue(10.0, 3, 0);
   v.setValue(11.0, 3, 1);
   v.setValue(12.0, 3, 2);
   v.setValue(13.0, 3, 3);

   MatrixSkewSymmetric resp(4);
   resp.setValue(1.0, 0, 0);
   resp.setValue(4.0, 1, 0);
   resp.setValue(5.0, 1, 1);
   resp.setValue(7.0, 2, 0);
   resp.setValue(8.0, 2, 1);
   resp.setValue(9.0, 2, 2);
   resp.setValue(10.0, 3, 0);
   resp.setValue(11.0, 3, 1);
   resp.setValue(12.0, 3, 2);
   resp.setValue(13.0, 3, 3);

   MatrixSkewSymmetric resp1(2);
   resp1.setValue(0.0, 0, 0);
   resp1.setValue(0.0, 1, 0);
   resp1.setValue(0.0, 1, 1);

   pmat::Vector vet{4};
   vet.setValue(2.0, 0);
   vet.setValue(8.0, 1);
   vet.setValue(10.0, 2);
   vet.setValue(14.0, 3);

   pmat::Vector respv{4};
   respv.setValue(-240.0, 0);
   respv.setValue(-186.0, 1);
   respv.setValue(0.0, 2);
   respv.setValue(410.0, 3);

   Vector vv{z * vet};

   MatrixSkewSymmetric a;
   a = (std::move(z));

   v.resize(2);

   EXPECT_TRUE(a == resp);
   EXPECT_TRUE(v == resp1);
   EXPECT_TRUE(vv == respv);
}
