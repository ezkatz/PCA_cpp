//
//  stats_util.cpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/2/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//

#include <iostream>
#include "math.h"
#include "../stats_util.hpp"

using std::cout; // TODO -- remove
using std::endl;

double computeMean(vector<double> X) {
    double sum = 0;
    for (double x_i :X) {
        sum += x_i;
    }
    return sum / X.size();
}

double computeCovariance(vector<double> X, vector<double> Y) {
    double mean_x = computeMean(X);
    double mean_y = computeMean(Y);

    double sum = 0;
    for (decltype(X.size()) i = 0; i < X.size(); ++i) {
        sum += (X[i] - mean_x) * (Y[i] - mean_y);
    }

    return sum / (X.size() - 1);
}

double computeVariance(vector<double> X) {
    double mean = computeMean(X);

    double variance = 0;
    for (double x_i : X) {
        double d = x_i - mean; // TODO -- rename?
        variance += d * d;
    }
    variance /= (X.size() - 1);

    return variance;
}

double computeStandardDeviation(vector<double> X) {
    return sqrt(computeVariance(X));
}

Matrix computeCovarianceMatrix(vector<vector<double>> v) {
    decltype(v.size()) i;
    decltype(v.size()) j;
    decltype(v.size()) size = v.size();

    Matrix res = Matrix(size, size);

    // upper triangle
    for (i = 0; i < size; ++i) {
        vector<double> v_i = v[i];
        for (j = i + 1; j < size; ++j) {
            res.set(i, j, computeCovariance(v_i, v[j]));
        }
    }

    // TODO -- decide if we should keep the diagonal separate from the
    // upper triangle or if we should merge
    // diagonal
    for (i = 0; i < size; ++i) {
        vector<double> v_i = v[i];
        res.set(i, i, computeCovariance(v_i, v_i));
    }

    // lower triangle
    for (i = 1; i < size; ++i) {
        for (j = 0; j < i; ++j) {
            res.set(i, j, res.get(j, i)); // copy reflection
        }
    }

    return res;
}
