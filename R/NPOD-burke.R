burke <- function(Psi) {
  # The burke function solves a convex optimization problem using an
  # interior-point Newton method. It iteratively updates lambda (lam)
  # while ensuring feasibility constraints. A damped step-size approach
  # is used to maintain numerical stability.

  Psi <- abs(Psi)

  # Check Input Specifications
  row <- nrow(Psi)
  col <- ncol(Psi)

  # Check that Psi has only non-negative entries.
  if (min(Psi) < 0) {
    stop("The matrix Psi has a negative element.")
  }

  # Check that Psi * ecol > 0.
  ecol <- matrix(1, nrow = col, ncol = 1) # Ensure ecol is a column vector
  Plam <- Psi %*% ecol
  if (min(Plam) <= 1.e-15) {
    stop("The vector Psi * ecol has a non-positive entry.")
  }

  # Stopping tolerance
  eps <- 1.e-8

  # Initialization
  sig <- 0
  erow <- rep(1, row)
  lam <- ecol
  w <- 1 / Plam
  Ptw <- t(Psi) %*% w
  shrink <- 2 * max(Ptw)
  lam <- lam * shrink
  Plam <- Plam * shrink
  w <- w / shrink
  Ptw <- Ptw / shrink
  y <- ecol - Ptw
  R <- erow - w * Plam
  normR <- max(abs(R))
  gap <- abs((sum(log(w)) + sum(log(Plam)))) / (1 + max(abs(sum(log(Plam)))))
  mu <- (t(lam) %*% y) / col
  mu <- mu[1]

  # Reporting structures
  normNiter <- c(normR)
  muiter <- c(mu)
  iter <- 0

  # Newton Iteration
  while (mu > eps || normR > eps || gap > eps) {
    iter <- iter + 1
    smu <- sig * mu

    # Form the matrix H
    inner <- lam / y
    H <- Psi %*% Diagonal(x = inner, n = col) %*% t(Psi) + Diagonal(x = Plam / w, n = row)
    H <- as.matrix(H)
    # Cholesky factorization of H
    cholH <- try(chol(H), silent = TRUE)
    if ("try-error" %in% class(cholH)) {
      lam <- lam / sum(lam)
      return(list(obj = sum(log(Psi %*% lam)), lam = lam))
    }
    # Right-hand side for dw
    smuyinv <- smu * (ecol / y)
    rhsdw <- erow / w - Psi %*% smuyinv

    # Compute dw
    dw <- solve(cholH, solve(t(cholH), rhsdw))

    # Compute dy
    dy <- -t(Psi) %*% dw

    # Compute dlam
    dlam <- smuyinv - lam - inner * dy

    # Damped Newton Step

    alfpri <- -1 / min(min(dlam / lam), -0.5)
    alfdual <- -1 / min(min(dy / y), -0.5)
    alfdual <- min(alfdual, -1 / min(min(dw / w), -0.5))
    alfpri <- min(1, 0.99995 * alfpri)
    alfdual <- min(1, 0.99995 * alfdual)

    # Newton Update
    if (alfpri == 1) {
      lam <- lam + dlam
    } else {
      lam <- lam + alfpri * dlam
    }


    if (alfdual == 1) {
      w <- w + dw
      y <- y + dy
    } else {
      w <- w + alfdual * dw
      y <- y + alfdual * dy
    }

    # Update mu and intermediary vectors and stopping criteria values
    mu <- (t(lam) %*% y) / col
    mu <- mu[1]
    Plam <- Psi %*% lam
    R <- erow - w * Plam
    Ptw <- Ptw - alfdual * dy
    normR <- max(abs(R))
    gap <- abs((sum(log(w)) + sum(log(Plam)))) / (1 + max(abs(sum(log(Plam)))))

    # Update sig
    if (mu < eps && normR > eps) {
      sig <- 1
    } else {
      c1 <- 2
      c2 <- 1.e+2
      sig <- min(c(0.3, max(c((1 - alfpri)^c1, (1 - alfdual)^c1, (normR - mu) / (normR + c2 * mu)))))
    }

    # Update reporting information
    muiter <- c(muiter, mu)
    normNiter <- c(normNiter, normR)
  }

  lam <- lam / row
  obj <- sum(log(Psi %*% lam))
  lam <- lam / sum(lam)

  list(fobj = obj, lambda = lam)
}


# library(Matrix) # For sparse matrix operations
#
# burke <- function(input_matrix) {
#   input_matrix <- abs(input_matrix)
#
#   # Input Validation
#   num_rows <- nrow(input_matrix)
#   num_cols <- ncol(input_matrix)
#
#   if (any(input_matrix < 0)) stop("The matrix has a negative element.")
#
#   row_sums <- input_matrix %*% matrix(1, nrow = num_cols, ncol = 1)
#   if (any(row_sums <= 1.e-15)) stop("The vector input_matrix * ecol has a non-positive entry.")
#
#   # Initialization
#   tolerance <- 1.e-8
#   damping_factor <- 0
#   row_sum_constraints <- rep(1, num_rows)
#   lambda_vec <- matrix(1, nrow = num_cols, ncol = 1)
#   weight_vec <- 1 / row_sums
#   transposed_weight_product <- t(input_matrix) %*% weight_vec
#   shrink_factor <- 2 * max(transposed_weight_product)
#
#   lambda_vec <- lambda_vec * shrink_factor
#   row_sums <- row_sums * shrink_factor
#   weight_vec <- weight_vec / shrink_factor
#   transposed_weight_product <- transposed_weight_product / shrink_factor
#
#   dual_variable_y <- matrix(1, nrow = num_cols, ncol = 1) - transposed_weight_product
#   residuals <- row_sum_constraints - weight_vec * row_sums
#   constraint_violation <- max(abs(residuals))
#   objective_gap <- abs((sum(log(weight_vec)) + sum(log(row_sums)))) / (1 + abs(sum(log(row_sums))))
#   duality_gap <- (t(lambda_vec) %*% dual_variable_y)[1] / num_cols
#
#   # Iteration Reporting
#   iter <- 0
#   norm_constraints <- c(constraint_violation)
#   duality_iterations <- c(duality_gap)
#
#   while (duality_gap > tolerance || constraint_violation > tolerance || objective_gap > tolerance) {
#     iter <- iter + 1
#     # Compute Hessian
#     scaling_factors <- as.vector(lambda_vec / dual_variable_y)
#     hessian <- input_matrix %*% Diagonal(x = scaling_factors) %*% t(input_matrix) + Diagonal(x = as.vector(row_sums / weight_vec))
#     chol_hessian <- try(chol(hessian), silent = TRUE)
#
#     if (inherits(chol_hessian, "try-error")) {
#       lambda_vec <- lambda_vec / sum(lambda_vec)
#       return(list(obj = sum(log(input_matrix %*% lambda_vec)), lambda = lambda_vec))
#     }
#
#     # Compute Newton Steps
#     scaled_dual_update <- damping_factor * duality_gap * (1 / dual_variable_y)
#     rhs_newton_step <-  damping_factor * duality_gap * (1 / dual_variable_y)
#     delta_weights <- solve(chol_hessian, rhs_newton_step)
#     delta_dual <- -t(input_matrix) %*% delta_weights
#     delta_lambda <- scaled_dual_update - lambda_vec - scaling_factors * delta_dual
#
#     # Damped Updates
#     alpha_primal <- min(1, 0.99995 * (-1 / min(c(delta_lambda / lambda_vec, 1))))
#     alpha_dual <- min(1, 0.99995 * (-1 / min(c(delta_dual / dual_variable_y, delta_weights / weight_vec, 1))))
#
#     lambda_vec <- lambda_vec + alpha_primal * delta_lambda
#     weight_vec <- weight_vec + alpha_dual * delta_weights
#     dual_variable_y <- dual_variable_y + alpha_dual * delta_dual
#
#     # Update Variables
#     duality_gap <- (t(lambda_vec) %*% dual_variable_y)[1] / num_cols
#     row_sums <- input_matrix %*% lambda_vec
#     residuals <- row_sum_constraints - weight_vec * row_sums
#     transposed_weight_product <- transposed_weight_product - alpha_dual * delta_dual
#     constraint_violation <- max(abs(residuals))
#     objective_gap <- abs((sum(log(weight_vec)) + sum(log(row_sums)))) / (1 + abs(sum(log(row_sums))))
#
#     # Update damping factor
#     damping_factor <- ifelse(duality_gap < tolerance && constraint_violation > tolerance, 1,
#                              min(0.3, max((1 - alpha_primal)^2, (1 - alpha_dual)^2, (constraint_violation - duality_gap) / (constraint_violation + 1e2 * duality_gap))))
#
#     # Logging (Optional)
#     norm_constraints <- c(norm_constraints, constraint_violation)
#     duality_iterations <- c(duality_iterations, duality_gap)
#   }
#
#   lambda_vec <- lambda_vec / sum(lambda_vec)
#   obj_value <- sum(log(input_matrix %*% lambda_vec))
#
#   list(fobj = obj_value, lambda = lambda_vec)
# }

