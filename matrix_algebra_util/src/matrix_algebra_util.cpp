//
//  matrix_algebra_util.cpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/3/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//


#include <iostream>
#include "math.h"
#include "../matrix_algebra_util.hpp"

double euclideanNorm(vector<double> X) {
    double sum = 0;
    for (int x_i : X) {
        sum += (x_i * x_i);
    }
    return sqrt(sum);
}

double eigenvalue(vector<double> original, vector<double> multiple) {
    double value = NAN;
    int size = static_cast<int>(original.size());
    for (int i = 0; i < size; ++i) {
        if (isnan(value)) {
            if (original[i] != 0) {
                value = multiple[i] / original[i];
            }
        } else {
            // TODO -- switch from division to multiplication? requires epsilon?
            if (original[i] == 0) {
                if (multiple[i] != 0) {
                    return NAN;
                }
            } else if (multiple[i]/original[i] != value) {
                return NAN;
            }
        }
    }

    return value;
}

vector<double> zeroMean(vector<double> v) {
    vector<double> res(v);
    double mean = computeMean(v);

    for (double& i : res) {
        i -= mean;
    }

    return res;
}
