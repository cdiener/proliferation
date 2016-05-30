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

template <typename C>
class index_sorter {
    public:
        index_sorter(C const& c) : c(c) {}
    bool operator()(std::size_t const& lhs, std::size_t const& rhs) const {
        return c[lhs] < c[rhs];
    }
    private:
        C const& c;
};

NumericVector max_reduce(NumericMatrix m, std::vector<int> idx) {
    NumericVector res(m.ncol(), 0.0);
    double val;
    
    int maxidx = 0;
    double maxval = mean(m(idx[0], _));
    for(int i=0; i<idx.size(); ++i) {
        val = mean(m(idx[i], _));
        if(val > maxval) { 
            maxval = val;
            maxidx = idx[i];
        }
    }
    res = m(maxidx, _);
    
    return res;
}

NumericVector mean_reduce(NumericMatrix m, std::vector<int> idx) {
    NumericVector res(m.ncol(), 0.0);
    
    for(int i=0; i<idx.size(); ++i) res += m(idx[i], _);
    
    return res/idx.size();
}

NumericVector median_reduce(NumericMatrix m, std::vector<int> idx) {
    NumericVector res(m.ncol(), 0.0);
    std::vector<double> means(idx.size());
    
    for(int i=0; i<idx.size(); ++i) means[i] = mean(m(idx[i], _));
    
    int i = idx.size()/2;
    if(idx.size() % 2 == 0) {
        std::nth_element(idx.begin(), idx.begin()+i-1, idx.end(), 
            index_sorter<std::vector<double> >(means));
        res += m(idx[i-1], _);
        std::nth_element(idx.begin(), idx.begin()+i, idx.end(), 
            index_sorter<std::vector<double> >(means));
        res += m(idx[i-1], _);
        res = 0.5*res;
    } else {
        std::nth_element(idx.begin(), idx.begin()+i-1, idx.end(), 
             index_sorter<std::vector<double> >(means));
        res += m(idx[i], _);
    }
    
    return res;
}

// [[Rcpp::export]]
NumericVector matrix_reduce(NumericMatrix m, IntegerVector idx, 
    IntegerVector factors, std::string method) {
    int minf = min(factors);
    int nf = max(factors) - min(factors) + 1;
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
        
        if(method == "max") res(fi, _) = max_reduce(m, map[fi]);
        else if(method == "median") res(fi, _) = median_reduce(m, map[fi]);
        else res(fi, _) = mean_reduce(m, map[fi]);
    }
    
    return res;
}
