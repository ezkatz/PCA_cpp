//
//  matrix_algebra_util.hpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/3/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//

#ifndef matrix_algebra_util_hpp
#define matrix_algebra_util_hpp

#include <stdio.h>
#include <vector>
#include "../stats_util/stats_util.hpp"

using std::vector;

double euclideanNorm(vector<double> X);
double eigenvalue(vector<double> original, vector<double> multiple);
vector<double> zeroMean(vector<double> v);

#endif /* matrix_algebra_util_hpp */
