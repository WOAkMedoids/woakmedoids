#' Whale Optimization Algorithm for K-Medoids Clustering
#'
#' This function implements the Whale Optimization Algorithm (WOA) for K-Medoids clustering.
#' Supported distance measures are Dynamic Time Warping (DTW) and Euclidean Distance (ED).
#'
#' @param data Data matrix
#' @param ClusNum Number of clusters
#' @param distance_method Distance calculation method, either "dtw" or "ed"
#' @param learned_w Window size for DTW (only used if distance_method is "dtw")
#' @param Max_iter Maximum number of iterations (default is 20, it can be adjusted according to the size of the dataset)
#' @param n Population size (number of whales, default is 5, itcan be adjusted according to the size of the dataset)
#' @return The `woa_clustering` object containing the clustering result and medoids
#' @references
#' Mirjalili, S., & Lewis, A. (2016). The whale optimization algorithm. Advances in engineering software, 95, 51-67.
#' @author Chenan Huang, Narumasa Tsutsumida
#' @export
#' @examples
#' # NOTE: This example only shows how to implement woa_kmedoids using sample data.
#' # Results do not suggest any meanings.
#' data(Lightning7)
#' Lightning7_data <- Lightning7[, -1]  # Remove the first column of classification data
#'   result <- woa_kmedoids(Lightning7_data, ClusNum = 7, distance_method = "dtw", learned_w = 5)
#'   print(result)


woa_kmedoids <- function(
    data,
    ClusNum,
    distance_method = c("dtw", "ed"),
    learned_w = NULL,  # This parameter is only used for DTW
    Max_iter = 20,
    n = 5) {

  # Calculate distance index based on triangular matrix storage
  calc_distance_index <- function(point1, point2, dist_obj) {
    if (point1 > point2) {
      temp <- point1
      point1 <- point2
      point2 <- temp
    }

    total <- 0
    row <- nrow(dist_obj)
    for (j in 0:(point1 - 1)) {
      total <- total + (ifelse(j == 0, 0, 1) * (row - j))
    }
    total <- total + (point2 - point1)
    return(dist_obj[total])
  }

  # Assign points to the closest medoid for each cluster
  assign_clusters <- function(medoids, dist_obj) {
    n <- attr(dist_obj, "Size")
    cluster_assignments <- numeric(n)

    for (point in 1:n) {
      if (!point %in% medoids) {
        distances <- sapply(medoids, function(medoid) {
          calc_distance_index(point, medoid, dist_obj)
        })
        closest_medoid <- medoids[which.min(distances)]
        cluster_assignments[point] <- which(medoids == closest_medoid)
      } else {
        cluster_assignments[point] <- which(medoids == point)
      }
    }

    return(cluster_assignments)
  }

  # Evaluate fitness of a solution using the totalDistance function implemented in C++
  fitness_function <- function(medoids, dist_mat) {
    clustering <- assign_clusters(medoids, dist_mat)
    cluster_sizes <- table(clustering)

    if (any(cluster_sizes < 2)) {
      return(Inf)
    }

    # Assuming totalDistance is a function accessible from R that was implemented in C++
    fitness_score <- totalDistance(dist_mat, nrow(dist_mat), medoids, clustering)
    return(fitness_score)
  }

  # Update positions by the whale optimization algorithm
  update_position <- function(i, Positions, best_solution, a, a2, b, n) {
    r1 <- runif(1)
    r2 <- runif(1)
    A <- 2 * a * r1 - a
    C <- 2 * r2
    l <- (a2 - 1) * runif(1) + 1

    if (runif(1) < 0.5) {
      if (abs(A) >= 1) {
        rand_leader_index <- sample(1:n, 1)
        X_rand <- Positions[rand_leader_index, ]
        D <- abs(C * X_rand - Positions[i, ])
        return(round(X_rand - A * D))
      } else {
        D <- abs(C * best_solution - Positions[i, ])
        return(round(best_solution - A * D))
      }
    } else {
      D_Leader <- abs(best_solution - Positions[i, ])
      return(round(D_Leader * exp(b * l) * cos(l * 2 * pi) + best_solution))
    }
  }

  # Calculate distance matrix using the selected method
  if (distance_method == "dtw") {
    dist_mat <- proxy::dist(data, method = dtw_basic, norm = "L2", symmetric = TRUE, window.size = learned_w)
  } else if (distance_method == "ed") {
    dist_mat <- proxy::dist(data, method = "Euclidean", norm = "L2", symmetric = TRUE)
  } else {
    stop("Unsupported distance method")
  }

  # Define bounds and parameters
  lb <- 1 # Lower bound
  ub <- nrow(data) # Upper bound
  b <- 1 # Parameter for position update

  # Initialize positions
  Positions <- t(replicate(n, sample(lb:ub, ClusNum)))

  # Main loop for WOA iterations
  for (t in 1:Max_iter) {
    # Calculate fitness values for all positions
    fitness_values <- apply(Positions, 1, function(pos) fitness_function(pos, dist_mat))
    best_solution <- Positions[which.min(fitness_values), ]

    # Update parameters
    a <- 2 - t * ((2) / Max_iter)
    a2 <- -1 + t * ((-1) / Max_iter)

    # Update positions
    results <- lapply(1:n, function(i) update_position(i, Positions, best_solution, a, a2, b, n))
    Positions <- do.call(rbind, results)

    # Ensure best solution remains in positions
    Positions[which.min(fitness_values), ] <- best_solution

    # Ensure positions are within bounds
    Positions[Positions < lb] <- lb
    Positions[Positions > ub] <- ub
    Positions <- t(apply(Positions, 1, function(x) {
      while (length(unique(x)) < length(x)) {
        x[duplicated(x)] <- sample(setdiff(lb:ub, x), sum(duplicated(x)))
      }
      return(x)
    }))
  }

  # Perform K-Medoids clustering using the best solution found by WOA
  woa_clustering <- cluster::pam(dist_mat, k = ClusNum, medoids = best_solution, diss = TRUE, do.swap = FALSE)

  # Return the woa_clustering object
  return(woa_clustering)
}
