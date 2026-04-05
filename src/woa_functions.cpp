#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// O(1) distance index calculation for triangular matrix storage
static inline double getDistanceHelper(const double *dist_ptr, int nrows, int point1, int point2) {
    if (point1 == point2) {
        return 0.0;
    }

    if (point1 > point2) {
        std::swap(point1, point2);
    }

    int total = (point1 - 1) * (2 * nrows - point1) / 2 + (point2 - point1);
    return dist_ptr[total - 1];
}

struct FitnessWorker : public Worker {
    const RVector<double> dist;
    const RVector<int> medoids;
    const int nrows;
    const int k;
    std::vector<int> &best_k_idx;
    std::vector<double> &min_dists;
    double total;

    FitnessWorker(const NumericVector dist_vec,
                  const IntegerVector medoids_vec,
                  int nrows,
                  int k,
                  std::vector<int> &best_k_idx,
                  std::vector<double> &min_dists)
        : dist(dist_vec),
          medoids(medoids_vec),
          nrows(nrows),
          k(k),
          best_k_idx(best_k_idx),
          min_dists(min_dists),
          total(0.0) {}

    FitnessWorker(const FitnessWorker &other, Split)
        : dist(other.dist),
          medoids(other.medoids),
          nrows(other.nrows),
          k(other.k),
          best_k_idx(other.best_k_idx),
          min_dists(other.min_dists),
          total(0.0) {}

    void operator()(std::size_t begin, std::size_t end) {
        const double *dist_ptr = dist.begin();
        double local_total = 0.0;

        for (std::size_t idx = begin; idx < end; ++idx) {
            int point = static_cast<int>(idx) + 1;
            double min_dist = R_PosInf;
            int best_k = 0;

            for (int m = 0; m < k; ++m) {
                double d = getDistanceHelper(dist_ptr, nrows, point, medoids[m]);
                if (d < min_dist) {
                    min_dist = d;
                    best_k = m;
                }
            }

            best_k_idx[idx] = best_k;
            min_dists[idx] = min_dist;
            local_total += min_dist;
        }

        total += local_total;
    }

    void join(const FitnessWorker &rhs) { total += rhs.total; }
};

//' Fitness function for WOA optimization (C++ implementation)
//'
//' Combines cluster assignment and total distance calculation.
//' Returns Inf if any cluster has fewer than 2 members.
//'
//' @param dist_vec Distance vector (from dist object)
//' @param nrows Number of rows in the original data
//' @param medoids Integer vector of medoid indices (1-indexed)
//' @return total_dist (double), Inf if invalid clustering
//' @keywords internal
// [[Rcpp::export]]
double fitnessFunction_cpp(NumericVector dist_vec, int nrows, IntegerVector medoids) {
    int k = medoids.size();

    std::vector<int> best_k_idx(nrows, 0);
    std::vector<double> min_dists(nrows, 0.0);

    FitnessWorker worker(dist_vec, medoids, nrows, k, best_k_idx, min_dists);
    parallelReduce(0, static_cast<std::size_t>(nrows), worker);

    std::vector<int> cluster_sizes(k, 0);
    for (int i = 0; i < nrows; ++i) {
        cluster_sizes[best_k_idx[i]]++;
    }

    for (int size : cluster_sizes) {
        if (size < 2) return R_PosInf;
    }

    return worker.total;
}

//' Deterministic serial fitness function for final evaluation
//'
//' @param dist_vec Distance vector (from dist object)
//' @param nrows Number of rows in the original data
//' @param medoids Integer vector of medoid indices (1-indexed)
//' @return total_dist (double), Inf if invalid clustering
//' @keywords internal
// [[Rcpp::export]]
double fitnessFunction_serial_cpp(NumericVector dist_vec, int nrows, IntegerVector medoids) {
    int k = medoids.size();
    double total_dist = 0.0;
    std::vector<int> cluster_sizes(k, 0);

    const double *dist_ptr = dist_vec.begin();

    for (int point = 1; point <= nrows; ++point) {
        double min_dist = R_PosInf;
        int best_k = 0;
        for (int m = 0; m < k; ++m) {
            double d = getDistanceHelper(dist_ptr, nrows, point, medoids[m]);
            if (d < min_dist) {
                min_dist = d;
                best_k = m;
            }
        }
        cluster_sizes[best_k]++;
        total_dist += min_dist;
    }

    for (int size : cluster_sizes) {
        if (size < 2) return R_PosInf;
    }

    return total_dist;
}


//' Check if two medoid solutions are identical
//'
//' @param medoids1 First medoid vector (1-indexed)
//' @param medoids2 Second medoid vector (1-indexed)
//' @return TRUE if identical (regardless of order), FALSE otherwise
//' @keywords internal
// [[Rcpp::export]]
bool medoidsEqual_cpp(IntegerVector medoids1, IntegerVector medoids2) {
    if (medoids1.size() != medoids2.size()) {
        return false;
    }

    int k = medoids1.size();

    IntegerVector sorted1 = clone(medoids1);
    IntegerVector sorted2 = clone(medoids2);
    std::sort(sorted1.begin(), sorted1.end());
    std::sort(sorted2.begin(), sorted2.end());

    for (int i = 0; i < k; ++i) {
        if (sorted1[i] != sorted2[i]) {
            return false;
        }
    }

    return true;
}

//' Project MDS coordinates to nearest unique sample indices
//'
//' Finds the nearest sample in the MDS embedded space for each medoid coordinate,
//' while enforcing a one-to-one mapping between medoids and samples (no duplicates).
//'
//' @param medoid_coords NumericMatrix of medoid coordinates in MDS space (ClusNum x mds_dim)
//' @param Z NumericMatrix of sample coordinates in MDS space (nrows x mds_dim)
//' @return IntegerVector of nearest unique sample indices (1-indexed, length = ClusNum)
//' @keywords internal
// [[Rcpp::export]]
IntegerVector projectToIndex_cpp(NumericMatrix medoid_coords, NumericMatrix Z) {
    int n_medoids = medoid_coords.nrow();
    int n_samples = Z.nrow();
    int n_dims = Z.ncol();

    if (medoid_coords.ncol() != n_dims) {
        stop("Dimension mismatch: medoid_coords and Z must have the same number of columns");
    }

    IntegerVector result(n_medoids);
    std::vector<bool> used(n_samples, false);

    for (int i = 0; i < n_medoids; i++) {
        double min_dist_sq = R_PosInf;
        int best_idx = -1;

        for (int j = 0; j < n_samples; j++) {
            if (used[j]) continue;

            double dist_sq = 0.0;
            for (int d = 0; d < n_dims; d++) {
                double diff = medoid_coords(i, d) - Z(j, d);
                dist_sq += diff * diff;
            }

            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_idx = j;
            }
        }

        if (best_idx == -1) {
            stop("Unique index allocation failed in projectToIndex_cpp.");
        }

        result[i] = best_idx + 1;
        used[best_idx] = true;
    }

    return result;
}
