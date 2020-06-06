//
//  stats_util.hpp
//  principal_component_analysis
//
//  Created by Eli Katz on 6/2/20.
//  Copyright © 2020 Eli Katz. All rights reserved.
//

#ifndef stats_util_hpp
#define stats_util_hpp

#include <stdio.h>
#include <vector>
#include "../matrix/matrix.hpp"

using std::vector;

double computeMean(vector<double> X);
double computeCovariance(vector<double> X, vector<double> Y);
double computeVariance(vector<double> X);
double computeStandardDeviation(vector<double> X);
Matrix computeCovarianceMatrix(vector<vector<double>> v);

#endif /* stats_util_hpp */
