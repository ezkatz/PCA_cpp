//
//  main.cpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/2/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//

/*
 1. Standardize the data. (with mean =0 and variance = 1)
 2. Compute the Covariance matrix of dimensions.
 3. Obtain the Eigenvectors and Eigenvalues from the covariance matrix (we can also use correlation matrix or even Single
  value decomposition, however in this post will focus on covariance matrix).
 4. Sort eigenvalues in descending order and choose the top k Eigenvectors that correspond to the k largest eigenvalues
  (k will become the number of dimensions of the new feature subspace k≤d, d is the number of original dimensions).
 5. Construct the projection matrix W from the selected k Eigenvectors.
 6. Transform the original data set X via W to obtain the new k-dimensional feature subspace Y.
 */

#include <iostream>
#include <cmath>
#include "stats_util/stats_util.hpp"
#include "matrix_algebra_util/matrix_algebra_util.hpp"
#include "matrix/matrix.hpp"

using std::cout;
using std::endl;

void printStats(vector<double> X) {
    double mean = computeMean(X);
    cout << "mean: " << mean << endl;

    double variance = computeVariance(X);
    cout << "variance: " << variance << endl;

    double std_dev = computeStandardDeviation(X);
    cout << "std_dev: " << std_dev << endl;
}

void printVector(vector<double> X) {
    cout << endl;
    for (double i : X) {
        cout << i << endl;
    }
}

void testEigenvector(vector<double> X, Matrix M) {
    vector<double> res = M.rightMultiply(X);
    double eigenval = eigenvalue(X, res);
    if (!isnan(eigenval)) {
        printVector(X);
        printVector(res);
        cout << "Eigenvector! eigenvalue: " << eigenval << endl;
    }
}

int main(int argc, const char * argv[]) {

    vector<double> X1{ 1, 2, 5, 6, 12, 15, 25, 45, 68, 67, 65, 98 };
    vector<double> X2{ 0, 8, 12, 20 };
    vector<double> X3{ 8, 9, 11, 12 };
    vector<double> X4{ 2, 2, -1 };
    vector<double> X5{ -1, 0, 2 };
    vector<double> X6{ -1, 1, 3 };
    vector<double> X7{ 0, 1, 0 };
    vector<double> X8{ 3, 2, 1 };

    vector<double> Y1{ 1, 3 };
    vector<double> Y2{ 3, 2 };

    Matrix M1 = Matrix(2, 2);
    M1.setRow(0, { 2, 3 });
    M1.setRow(1, { 1, -5 });

    Matrix M2 = Matrix(2, 3);
    M2.setRow(0, { 4, 3, 6 });
    M2.setRow(1, { 1, -2, 3 });

    Matrix M3 = Matrix(3, 3);
    M3.setRow(0, { 3, 0, 1 });
    M3.setRow(1, { -4, 1, 2 });
    M3.setRow(2, { -6, 0, -2 });

    Matrix M4 = Matrix(3, 3);
    M4.setRow(0, { 1, 1, 1 });
    M4.setRow(1, { 0, 1, 1 });
    M4.setRow(2, { 0, 0, 1 });

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    printStats(X1);
//    cout << endl;
//
//    printStats(X2);
//    cout << endl;
//
//    printStats(X3);
//    cout << endl;

//    double covariance = computeCovariance(X2, X3);
//    cout << "covariance: " << covariance << endl;
//
//    double norm = euclideanNorm(Y2);
//    cout << "norm of Y2: " << norm << endl;
//
//    vector<double> res = M1.rightMultiply(Y1);
//    printVector(res);
//
//    Matrix matrixMult = M1.rightMultiply(M2);
//    matrixMult.print();
//
//    testEigenvector(X4, M3);
//    testEigenvector(X5, M3);
//    testEigenvector(X6, M3);
//    testEigenvector(X7, M3);
//    testEigenvector(X8, M3);
//
//    vector<double> zzz = zeroMean(Y2);
//    printVector(Y2);
//    printVector(zzz);
//    cout << endl;
//
////    vector<vector<double>> cov{
////        { 1, -1, 4 },
////        { 2, 1, 3 },
////        { 1, 3, -1 },
////    };
////    vector<vector<double>> cov2{
////        { 90, 90, 60, 60, 30 },
////        { 60, 90, 60, 60, 30 },
////        { 90, 30, 60, 90, 30 }
////    };
//
////    Matrix cov_res = computeCovarianceMatrix(cov);
////    cov_res.print();
////
////    Matrix cov2_res = computeCovarianceMatrix(cov2);
////    cov2_res.print();
//
//    vector<vector<double>> cov3{
//        { .69, -1.31, .39, .09, 1.29, .49, .19, -.81, -.31, -.71 },
//        { .49, -1.21, .99, .29, 1.09, .79, -.31, -.81, -.31, -1.01 }
//    };
//    Matrix cov3_res = computeCovarianceMatrix(cov3);
//    cov3_res.print();
//
//    Matrix M5 = Matrix(4, 4);
//    M5.setRow(0, { 1, -3, 1, -2 });
//    M5.setRow(1, { 2, -5, -1, -2 });
//    M5.setRow(2, { 0, -4, 5, 1 });
//    M5.setRow(3, { -3, 10, -6, 8 });
//
//    M5.print();
//    cout << endl;
//    Matrix echelon = M5.echelonForm();
//    echelon.print();

    Matrix A = Matrix(4, 4);
    A.setRow(0, { 1, 2, 1, 1 });
    A.setRow(1, { 1, 3, 1, 1 });
    A.setRow(2, { 1, 1, 4, 1 });
    A.setRow(3, { 5, 1, 1, 1 });

    Matrix B = Matrix(4, 4);
    B.setRow(0, { 1, 2, 3, 4 });
	B.setRow(1, { 5, 6, 7, 8 });
	B.setRow(2, { 9, 0, 1, 2 });
	B.setRow(3, { 8, 9, 3, 4 });

//	Matrix multip = A.strassenMultiply(B);
//	multip.print();
//
//	Matrix multip2 = A.rightMultiply(B);
//	multip2.print();

	Matrix C = Matrix(2, 2);
	C.setRow(0, { 5, -3 });
	C.setRow(1, { -6, 2 });
	vector<double> eigenvalues = C.computeEigenvalues();

	for (double d : eigenvalues) {
		cout << "Eigenvalue: " << d << endl;
	}

    return 0;
}
