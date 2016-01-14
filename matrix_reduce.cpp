/*
 * matrix_reduce.cpp
 * 
 * Copyright 2016 Christian Diener <mail[at]cdiener.com>
 * 
 * MIT license. See LICENSE for more information.
 */

#include <Rcpp.h>
#include <vector>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

NumericVector mean_reduce(NumericMatrix m, std::vector<int> idx) {
    NumericVector res(m.ncol(), 0.0);
    
    for(int i=0; i<idx.size(); ++i) res += m(idx[i], _);
    
    return res/idx.size();
}

// [[Rcpp::export]]
NumericVector matrix_reduce(NumericMatrix m, IntegerVector idx, IntegerVector factors) {
    int nf = max(factors)+1;
    int minf = min(factors);
    std::vector< std::vector<int> > map(nf);
    
    //Rcout<<"Building factor map..."<<std::endl;
    int fi;
    for(int i=0; i<factors.size(); ++i) {
        fi = factors[i] - minf;
        map[fi].push_back(idx[i]-1);
    }
    
    //Rcout<<"Reducing matrix..."<<std::endl;
    Progress p(nf, true);
    NumericMatrix res(nf, m.ncol());
    for(int fi=0; fi<nf; ++fi) {
        p.increment();
        res(fi, _) = mean_reduce(m, map[fi]);
    }
    
    return res;
}
