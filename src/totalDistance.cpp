#include <Rcpp.h>
using namespace Rcpp;

// O(1) distance index calculation for triangular matrix storage
// Previous version used O(n) loop, now uses direct formula
double getDistance(NumericVector dist_vec, int nrows, int point1, int point2) {

    if (point1 == point2) {
        return 0;
    }
    
    // Ensure point1 < point2
    if (point1 > point2) {
        std::swap(point1, point2);
    }

    // O(1) direct formula for triangular matrix index
    // Replaces O(n) loop: for (j=0; j<point1; j++) total += (j==0?0:1)*(nrows-j)
    // Formula: (point1-1) * (2*nrows - point1) / 2 + (point2 - point1)
    int total = (point1 - 1) * (2 * nrows - point1) / 2 + (point2 - point1);

    // Return distance value (1-indexed in R, so subtract 1)
    return dist_vec[total - 1];
}

// [[Rcpp::export]]
double totalDistance(NumericVector dist_vec, int nrows, IntegerVector medoids, IntegerVector clustering) {
    double totaldis = 0.0;
    for (int i = 0; i < nrows; ++i) {
        int medoid_index = clustering[i];
        int medoid = medoids[medoid_index - 1];
        int point = i + 1;
        totaldis += getDistance(dist_vec, nrows, point, medoid);
    }
    return totaldis;
}
