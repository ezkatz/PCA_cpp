//
//  matrix.cpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/3/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <exception>
#include "../matrix.hpp"

using std::cout;
using std::endl;
using std::exception;

// M x N matrix
// rows = M, cols = N
Matrix::Matrix(unsigned long rows, unsigned long cols) {
    M_ = rows;
    N_ = cols;
    for (int i = 0; i < rows; ++i) {
        vector<double> v(cols, 0);
        data.push_back(v);
    }
}

void Matrix::setRow(unsigned long row, vector<double> v) {
    // TODO -- check dimensions? Or let runtime error
    data[row] = v;
}

void Matrix::setCol(unsigned long col, vector<double> v) {
	// TODO -- check dimensions? Or let runtime error
    for (int i = 0; i < M_; ++i) {
        data[i][col] = v[i];
    }
}

void Matrix::set(unsigned long row, unsigned long col, double val) {
	// TODO -- check dimensions? Or let runtime error
    data[row][col] = val;
}

// V * M
// [1 x M][M x N] = [1 x N]
//vector<int> Matrix::leftMultiply(vector<int> v) {
//    vector<int> result;
//
//    return result;
//}

class InvalidDimensions : public exception {
	private:
		const char* msg = "Invalid dimensions for matrix multiplication";
	public:
		const virtual char* what() const throw () {
			return msg;
		}
};

// M * V
// [M x N][N x 1] = [M x 1]
vector<double> Matrix::rightMultiply(vector<double> v) {
    int size = static_cast<int>(v.size());

    if (size != N_) {
        throw InvalidDimensions();
    }

    vector<double> result;
    double sum;
    for (int i = 0; i < M_; ++i) {
        sum = 0;
        vector<double> x = data[i];
        for (int j = 0; j < N_; ++j) {
            sum += x[j] * v[j];
        }
        result.push_back(sum);
    }
    return result;
}

//Matrix Matrix::leftMultiply(Matrix m) {
//    return m;
//}

vector<double> Matrix::getColumn(unsigned long col) {
    vector<double> v;
    for (int i = 0; i < M_; ++i) {
        v.push_back(data[i][col]);
    }
    return v;
}

double dotProduct(vector<double> v, vector<double> y) {
    auto size = v.size();

    double sum = 0;
    for (decltype(v.size()) i = 0; i < size; ++i) {
        sum += v[i] * y[i];
    }
    return sum;
}

// [Low, high)
Matrix Matrix::copySubmatrix(int low_y, int high_y, int low_x, int high_x) {
	int rows = high_y - low_y;
	int cols = high_x - low_x;
	Matrix res = Matrix(rows, cols);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			res.data[i][j] = data[i + low_y][j + low_x];
		}
	}

	return res;
}

Matrix Matrix::add(Matrix M, bool positive) {
	Matrix res = Matrix(M_, N_);
	for (int i = 0; i < M_; ++i) {
		for (int j = 0; j < N_; ++j) {
			if (positive) {
				res.data[i][j] = data[i][j] + M.data[i][j];
			} else {
				res.data[i][j] = data[i][j] - M.data[i][j];
			}
		}
	}
	return res;
}

// TODO -- use new and add delete calls? or keep as is?
// TODO -- make this generic for Matrix multiplication other than [M x M][M x M]
Matrix Matrix::strassenMultiply(Matrix B) {
	int half = M_/2;

	if (M_ == 1) { // base case // TODO -- need to change when no longer generic
		Matrix o = Matrix(1, 1);
		o.set(0, 0, data[0][0] * B.data[0][0]);
		return o;
	}

	// A = 	A_11 A_12	B = B_11 B_12
	// 		A_21 A_22       B_21 B_22
	Matrix A_11 = copySubmatrix(0, half, 0, half);
	Matrix A_12 = copySubmatrix(0, half, half, M_);
	Matrix A_21 = copySubmatrix(half, M_, 0, half);
	Matrix A_22 = copySubmatrix(half, M_, half, M_);

	Matrix B_11 = B.copySubmatrix(0, half, 0, half);
	Matrix B_12 = B.copySubmatrix(0, half, half, M_);
	Matrix B_21 = B.copySubmatrix(half, M_, 0, half);
	Matrix B_22 = B.copySubmatrix(half, M_, half, M_);

	Matrix S_1 = B_12.add(B_22, false); //  B_12 - B_22
	Matrix S_2 = A_11.add(A_12, true);  //  A_11 + A_12
	Matrix S_3 = A_21.add(A_22, true);  //  A_21 + A_22
	Matrix S_4 = B_21.add(B_11, false); //  B_21 - B_11
	Matrix S_5 = A_11.add(A_22, true);  //  A_11 + A_22
	Matrix S_6 = B_11.add(B_22, true);  //  B_11 + B_22
	Matrix S_7 = A_12.add(A_22, false); //  A_12 - A_22
	Matrix S_8 = B_21.add(B_22, true);  //  B_21 + B_22
	Matrix S_9 = A_11.add(A_21, false); //  A_11 - A_21
    Matrix S_10 = B_11.add(B_12, true); //  B_11 + B_12

	// use strassen recursively
    Matrix P_1 = A_11.strassenMultiply(S_1); // A_11 * S_1
    Matrix P_2 = S_2.strassenMultiply(B_22); // S_2 * B_22
    Matrix P_3 = S_3.strassenMultiply(B_11); // S_3 * B_11
    Matrix P_4 = A_22.strassenMultiply(S_4); // A_22 * S_4
    Matrix P_5 = S_5.strassenMultiply(S_6);  // S_5 * S_6
    Matrix P_6 = S_7.strassenMultiply(S_8);  // S_7 * S_8
    Matrix P_7 = S_9.strassenMultiply(S_10); // S_9 * S_10

	Matrix res = Matrix(M_, M_);

	for (int i = 0; i < M_; ++i) {
		for (int j = 0; j < M_; ++j) {
			double val;

			if (i < half) {
				if (j < half) { //	C_11 = P_5 + P_4 - P_2 + P_6
					val = P_5.data[i][j]
						+ P_4.data[i][j]
						- P_2.data[i][j]
						+ P_6.data[i][j];
				} else { //	C_12 = P_1 + P_2
					val = P_1.data[i][j - half]
						+ P_2.data[i][j - half];
				}
			} else {
				if (j < half) { //	C_21 = P_3 + P_4
					val = P_3.data[i - half][j]
						+ P_4.data[i - half][j];
				} else { //	C_22 = P_5 + P_1 - P_3 - P_7
					val = P_5.data[i - half][j - half]
						+ P_1.data[i - half][j - half]
						- P_3.data[i - half][j - half]
						- P_7.data[i - half][j - half];
				}
			}

			res.data[i][j] = val;
		}
	}

	return res;
}

// [M x N][N x P] = [M x P]
Matrix Matrix::rightMultiply(Matrix m) {
    Matrix res = Matrix(M_, m.N_);
    for (int i = 0; i < N_; ++i) {
        vector<double> v = data[i];
        vector<double> y; // rename
        for (int j = 0; j < m.N_; ++j) {
            y.push_back(dotProduct(v, m.getColumn(j)));
        }
        res.setRow(i, y);
    }
    return res;
}

double Matrix::get(unsigned long row, unsigned long col) {
    return data[row][col];
}

void Matrix::print() {
	cout << "__";
	for (int i = 1; i < M_; ++i) {
		cout << "\t";
	}
	cout << "  __" << endl;
    for (int i = 0; i < M_; ++i) {
    	cout << "| ";
        vector<double> d = data[i];
        for (int j = 0; j < N_; ++j) {
            cout << d[j];
            if (j < N_ - 1) cout << "\t";
        }
        cout << " |";
        cout << endl;
    }
    cout << "￣";
	for (int i = 1; i < M_; ++i) {
		cout << "\t";
	}
	cout << "  ￣" << endl;
}

double computeDeterminant() {
	// make upper triangle
	// multiply diagonal entries
	// divide by scaling factors used in transformation
    return 0;
}

bool Matrix::isUpperTriangle() {
    for (int i = 0; i < N_; ++i) {
        if (data[i][i] == 0) {
            return false;
        }
        for (int j = 0; j < i; ++j) {
            if (data[i][j] != 0) {
                return false;
            }
        }
    }

    return true;
}

void Matrix::swapRows(unsigned long row1, unsigned long row2) {
    vector<double> temp = data[row1];
    data[row1] = data[row2];
    data[row2] = temp;
}

void Matrix::scaleRow(unsigned long row, double scalar) {
    vector<double> x = data[row];
    for (double& x_i : x) {
        x_i *= scalar;
    }
    setRow(row, x);
}

// subtract row2 from row1
void Matrix::subtractRow(unsigned long row1, unsigned long row2) {
    for (int i = 0; i < N_; ++i) {
        data[row1][i] -= data[row2][i];
    }
}

Matrix Matrix::echelonForm() {
    Matrix copy = Matrix(*this); // how does this work? Haven't defined copy constructor
    // this doesn't account for having to swap if a row doesn't have the proposed pivot column
    for (int z = 0; z < M_; ++z) {
        for (int i = z + 1; i < M_; ++i) {
            if (copy.data[i][z] != 0) { // TODO


//                copy.scaleRow(i, copy.data[z][z]/copy.data[i][z]);
//                copy.subtractRow(i, z);
            }
            copy.print();
            cout << endl;
        }
    }
    return copy;
}

Matrix Matrix::rowReduce() {

    return *this;
}

vector<double> Matrix::computeEigenvalues() {
    vector<double> v;
    if (M_ == N_) { // only if square matrix
    	// TODO -- generalize for larger matrices
    	if (2 == M_) {
    		// det ( w-l x   ) = 0
    		//     ( y	 z-l )
    		//
    		// (w-l)(z-l) - xy = 0
    		// (wz -lz -wl +l^2 - xy) = 0
    		// l^2 - (w+z)l + (wd-xy) = 0
    		int w = data[0][0];
    		int x = data[0][1];
    		int y = data[1][0];
    		int z = data[1][1];

    		// (-b +/- sqrt(b^2 - 4ac))/2a
    		// double a = 1;
    		// (-b +/- sqrt(b^2 - 4c))/2
    		double negB = (w + z); // b = -(w + z), negB = w + z
    		double c = (w * z) - (x * y); // (wz - xy)

    		double d = sqrt((negB * negB) - (4 * c));
    		v.push_back((negB + d)/2);
    		v.push_back((negB - d)/2);
    	}
    }
    return v;
}
