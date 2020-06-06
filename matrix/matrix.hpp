//
//  matrix.hpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/3/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdio.h>

#include <vector>
using std::vector;

class Matrix {
    private:
        unsigned long M_;
        unsigned long N_;
        vector<vector<double>> data;

        vector<double> getColumn(unsigned long col);

        void swapRows(unsigned long row1, unsigned long row2);
        void scaleRow(unsigned long row, double scalar);
        void subtractRow(unsigned long row1, unsigned long row2);

//        void add(Matrix orig, Matrix& new_, int low_x, int high_x, int low_y, int high_y, bool positive);
        Matrix copySubmatrix(int low_x, int high_x, int low_y, int high_y);
        Matrix add(Matrix M, bool positive);
    public:
        Matrix(unsigned long rows, unsigned long cols);
//        ~Matrix();

        void setRow(unsigned long row, vector<double> v);
        void setCol(unsigned long col, vector<double> v);
        void set(unsigned long row, unsigned long col, double val);

//        vector<double> leftMultiply(vector<double> v);
        vector<double> rightMultiply(vector<double> v);

//        Matrix leftMultiply(Matrix m);
        Matrix rightMultiply(Matrix m);

        double get(unsigned long row, unsigned long col);
    // TODO -- getRow
    // TODO -- getCol

        vector<double> computeEigenvalues();

        bool isUpperTriangle(); // TODO -- keep public?
        Matrix echelonForm();
        Matrix rowReduce();
        Matrix strassenMultiply(Matrix m);

        void print();
};

#endif /* matrix_hpp */
