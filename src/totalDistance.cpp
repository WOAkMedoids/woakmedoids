#include <Rcpp.h>
using namespace Rcpp;

double getDistance(NumericVector dist_vec, int nrows, int point1, int point2) {

    if (point1 == point2) {
        return 0;
    }
    // point1 < point2
    if (point1 > point2) {
        std::swap(point1, point2);
    }

    // Calculate the index in a one-dimensional array
    int total = 0;
    for (int j = 0; j < point1; ++j) {
        total += (j == 0 ? 0 : 1) * (nrows - j);
    }
    total += (point2 - point1);

    // Return distance value
    return dist_vec(total - 1);
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
