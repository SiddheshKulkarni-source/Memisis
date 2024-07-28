#' Multiview Learner Function
#'
#' This function generates synthetic multi-view data and performs various analyses, including plotting correlations and activity patterns.
#'
#' @param mainDir Character. The main directory where plots and data will be saved.
#' @param subDir Character. The subdirectory within \code{mainDir} for saving plots and data.
#' @param seed.number Numeric. Seed for random number generation.
#' @param M Integer. Number of views.
#' @param p_m Integer vector. Number of features per view.
#' @param n Integer. Train set size.
#' @param nTest Integer. Test set size.
#' @param params List. A list of parameters including:
#' \itemize{
#'   \item \code{fixed_act_pattern} Logical. Whether to use a fixed activity pattern.
#'   \item \code{K0} Integer. Number of factors.
#'   \item \code{nPatternJ} Integer. Number of joint patterns.
#'   \item \code{J0_m} Integer vector. Number of groups in the joint component per view.
#'   \item \code{K0_m} Integer vector. Number of factors in view-specific components.
#'   \item \code{S0_m} Integer vector. Number of groups in view-specific components.
#'   \item \code{nu2_0} Numeric. Variability in the joint loadings.
#'   \item \code{pi2_0} Numeric. Variability in the view-specific loadings.
#'   \item \code{snr_m} Numeric vector. Signal-to-noise ratio for each view.
#'   \item \code{snr_y} Numeric. Signal-to-noise ratio for the response.
#'   \item \code{decaying_loadings} Logical. Whether to enforce loadings decay.
#'   \item \code{decaying_type} Character. Type of decay ('sqrt' or 'log').
#'   \item \code{pJ_zero} Numeric. Probability of group-wise zero in joint component.
#'   \item \code{pJ_sign} Numeric. Probability of group-wise sign-switch in joint component.
#'   \item \code{pJ_zeroEntr} Numeric. Probability of entry-wise zero in joint component.
#'   \item \code{pS_zero} Numeric. Probability of group-wise zero in specific components.
#'   \item \code{pS_sign} Numeric. Probability of group-wise sign-switch in specific components.
#'   \item \code{pS_zeroEntr} Numeric. Probability of entry-wise zero in specific components.
#' }
#' @param Activity_pattern_plot Logical. Whether to save the activity pattern plot.
#' @param Emprical_correlation_plot Logical. Whether to save the empirical correlation plot.
#' @param Crossprd_plot Logical. Whether to save crossproduct plots.
#' @param Data_save Logical. Whether to save the generated data.
#'
#' @return A list containing the generated data and various parameters used in the simulation.
#'
#' @examples
#' mainDir <- tempdir()
#' subDir <- 'plots_simulated_data'
#' seed.number <- 1234
#' M <- 2
#' p_m <- c(500, 700)
#' n <- 100
#' nTest <- 100
#' params <- list(
#'   fixed_act_pattern = TRUE,
#'   K0 = 10,
#'   nPatternJ = 5,
#'   J0_m = c(4, 5),
#'   K0_m = c(5, 6),
#'   S0_m = c(4, 5),
#'   nu2_0 = 0.1,
#'   pi2_0 = 0.1,
#'   snr_m = NULL,
#'   snr_y = 1,
#'   decaying_loadings = FALSE,
#'   decaying_type = 'sqrt',
#'   pJ_zero = 0.6,
#'   pJ_sign = 0.4,
#'   pJ_zeroEntr = 0.5,
#'   pS_zero = 0.6,
#'   pS_sign = 0.4,
#'   pS_zeroEntr = 0.5
#' )
#' Activity_pattern_plot <- TRUE
#' Emprical_correlation_plot <- TRUE
#' Crossprd_plot <- TRUE
#' Data_save <- TRUE
#' Data <- Multiview_learner(mainDir, subDir, seed.number, M, p_m, n, nTest, params, Activity_pattern_plot, Emprical_correlation_plot, Crossprd_plot, Data_save)
#'
#' @export
Multiview_learner <- function(mainDir = '/Users/siddheshkulkarni/Desktop',
                              subDir = 'plots_simulated_data',
                              seed.number = 1234,
                              M = 3, # Number of views
                              p_m = c(500, 700, 1000),  # Number of features per view
                              n = 100,  # Train set size
                              nTest = 100,  # Test set size
                              params = list(
                                fixed_act_pattern = TRUE,  # Set a custom activity pattern in the joint component
                                K0 = 10,  # Number of factors
                                nPatternJ = 5,
                                J0_m = c(5, 6, 7),  # Number of groups in joint component
                                K0_m = c(5, 7, 10),  # Number of factors in view-specific components
                                S0_m = c(6, 7, 8),  # Number of groups in view-specific components
                                nu2_0 = 0.1,  # Input the joint loadings variability
                                pi2_0 = 0.1,  # Input the view-specific loadings variability
                                snr_m = NULL,  # Input the view-specific signal-to-noise ratio (mean)
                                snr_y = 1,  # Input signal-to-noise ratio
                                decaying_loadings = FALSE,  # Enforce loadings decay (re-weighting importance of latent axis of variation)
                                decaying_type = 'sqrt',  # One of 'sqrt' and 'log'
                                pJ_zero = 0.6,  # Probability of group-wise zero in joint component
                                pJ_sign = 0.4,  # Probability of group-wise sign-switch in joint component
                                pJ_zeroEntr = 0.5,  # Probability of entry-wise zero in joint component
                                pS_zero = 0.6,  # Probability of group-wise zero in specific components
                                pS_sign = 0.4,  # Probability of group-wise sign-switch in specific components
                                pS_zeroEntr = 0.5  # Probability of entry-wise zero in specific components
                              ),
                              Activity_pattern_plot = TRUE,  # Save the Activity Pattern Plot
                              Emprical_correlation_plot = TRUE, # Save the Empirical correlation plot
                              Crossprd_plot = TRUE, # Save crossproduct plots
                              Data_save = TRUE  # Save the data in directory
) {
  # Create directory for saving plots
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

  # Set random seed for reproducibility
  set.seed(seed.number)

  # Initialize snr_m if not provided
  if (is.null(params$snr_m)) {
    params$snr_m <- rep(1, M)
  }

  # Generate activity pattern combinations
  if (params$fixed_act_pattern == TRUE) {
    which_comb <- expand.grid(replicate(M + 1, c(0, 1), simplify = FALSE))
    which_comb <- apply(which_comb[rowSums(which_comb) > 1, ], 1, paste, collapse = "")
    which_comb <- sample(which_comb, params$nPatternJ)
    value_comb <- getBlocksDim(params$nPatternJ, params$K0, alpha = 5)
    act_J <- sapply(which_comb, function(v) as.numeric(unlist(strsplit(v, ""))))
    act_J <- t(apply(act_J, 1, function(v) rep(v, value_comb)))
  } else {
    stop("Non-fixed activity patterns are not supported in this implementation.")
  }

  # Generate dimensions for joint and specific components
  dim_J <- lapply(1:M, function(m) getBlocksDim(params$J0_m[m], p_m[m]))
  dim_S <- lapply(1:M, function(m) getBlocksDim(params$S0_m[m], p_m[m]))

  # Initialize group-wise sparsity pattern for joint component
  group_sprs_J <- list()
  for (m in 1:M) {
    group_sprs_J[[m]] <- matrix(rbinom(params$J0_m[m] * params$K0, 1, 1 - params$pJ_zero * (sum(act_J[m,]) / params$K0)), ncol = params$K0) *
      (2 * matrix(rbinom(params$J0_m[m] * params$K0, 1, 1 - params$pJ_sign), ncol = params$K0) - 1)

    # Ensure no active column is zeros only
    col_wrong <- apply(group_sprs_J[[m]], 2, function(v) min(v == 0)) * act_J[m,]
    while (sum(col_wrong) > 0) {
      group_sprs_J[[m]][, as.logical(col_wrong)] <- matrix(rbinom(params$J0_m[m] * sum(col_wrong), 1, 1 - params$pJ_zero * (sum(act_J[m, ]) / params$K0)), ncol = sum(col_wrong)) *
        (2 * matrix(rbinom(params$J0_m[m] * sum(col_wrong), 1, 1 - params$pJ_sign), ncol = sum(col_wrong)) - 1)
      col_wrong <- apply(group_sprs_J[[m]], 2, function(v) min(v == 0)) * act_J[m, ]
    }

    # Ensure no group is zeros only
    grp_wrong <- apply(group_sprs_J[[m]] * (rep(1, params$J0_m[m]) %x% t(act_J[m,])), 1, function(v) min(v == 0))
    while (sum(grp_wrong) > 0) {
      group_sprs_J[[m]][grp_wrong, ] <- act_J[m,] * (2 * matrix(rbinom(sum(grp_wrong) * params$K0, 1, 1 - params$pJ_sign), ncol = params$K0) - 1)
    }
  }

  # Initialize entry-wise sparsity pattern for joint component
  elem_sprs_J <- lapply(1:M, function(m) matrix(rbinom(p_m[m] * params$K0, 1, 1 - params$pJ_zeroEntr), ncol = params$K0))
  for (m in 1:M) {
    sprs_tot <- elem_sprs_J[[m]] * apply(group_sprs_J[[m]], 2, function(v) rep(v, dim_J[[m]]))
    if (max(apply(sprs_tot, 1, function(v) min(v == 0))) > 0) {
      elem_sprs_J[[m]][which(apply(sprs_tot, 1, function(v) min(v == 0)) > 0), act_J[m, ]] <- 1
    }
  }

  # Initialize group-wise sparsity pattern for specific components
  group_sprs_S <- list()
  for (m in 1:M) {
    group_sprs_S[[m]] <- matrix(rbinom(params$S0_m[m] * params$K0_m[m], 1, 1 - params$pS_zero), ncol = params$K0_m[m]) *
      (2 * matrix(rbinom(params$S0_m[m] * params$K0_m[m], 1, 1 - params$pS_sign), ncol = params$K0_m[m]) - 1)

    # Ensure no group is zeros only
    while (max(apply(group_sprs_S[[m]], 1, function(v) min(v == 0))) > 0) {
      group_sprs_S[[m]] <- matrix(rbinom(params$S0_m[m] * params$K0_m[m], 1, 1 - params$pS_zero), ncol = params$K0_m[m]) *
        (2 * matrix(rbinom(params$S0_m[m] * params$K0_m[m], 1, 1 - params$pS_sign), ncol = params$K0_m[m]) - 1)
    }
  }

  # Initialize entry-wise sparsity pattern for specific components
  elem_sprs_S <- lapply(1:M, function(m) matrix(rbinom(p_m[m] * params$K0_m[m], 1, 1 - params$pS_zeroEntr), ncol = params$K0_m[m]))
  for (m in 1:M) {
    sprs_tot <- elem_sprs_S[[m]] * apply(group_sprs_S[[m]], 2, function(v) rep(v, dim_S[[m]]))
    if (max(apply(sprs_tot, 1, function(v) min(v == 0))) > 0) {
      elem_sprs_S[[m]][which(apply(sprs_tot, 1, function(v) min(v == 0)) > 0), ] <- 1
    }
  }

  # Generate the true loading matrices
  mean_J <- lapply(1:M, function(m) (-1) ^ c(1:params$J0_m[m]) * rbeta(params$J0_m[m], shape1 = 5, shape2 = 3) / sqrt(params$K0))
  mean_S <- lapply(1:M, function(m) (-1) ^ c(1:params$S0_m[m]) * rbeta(params$S0_m[m], shape1 = 5, shape2 = 3) / sqrt(params$K0_m[m]))
  nu2_m <- params$nu2_0 / rep(params$K0, M)
  pi2_m <- params$nu2_0 / params$K0_m

  # Assign features to a group
  col_J <- sapply(1:M, function(m) unlist(sapply(1:params$J0_m[m], function(l) rep(l, dim_J[[m]][l]))))
  col_S <- sapply(1:M, function(m) unlist(sapply(1:params$S0_m[m], function(l) rep(l, dim_S[[m]][l]))))

  # Sample the full loadings matrices
  load_J <- lapply(1:M, function(m)
    matrix(mean_J[[m]][col_J[[m]]], nrow = p_m[m], ncol = params$K0) +
      sqrt(nu2_m[m]) * matrix(rnorm(p_m[m] * params$K0), nrow = p_m[m], ncol = params$K0))
  load_S <- lapply(1:M, function(m)
    matrix(mean_S[[m]][col_S[[m]]], nrow = p_m[m], ncol = params$K0_m[m]) +
      sqrt(pi2_m[m]) * matrix(rnorm(p_m[m] * params$K0_m[m]), nrow = p_m[m], ncol = params$K0_m[m])
  )

  # Enforce desired activity in the joint loadings
  for (m in 1:M) { load_J[[m]][, !act_J[m,]] <- 0 }
  # Enforce group-wise sparsity pattern
  for (m in 1:M) { load_J[[m]] <- load_J[[m]] * apply(group_sprs_J[[m]], 2, function(v) rep(v, dim_J[[m]])) }
  for (m in 1:M) { load_S[[m]] <- load_S[[m]] * apply(group_sprs_S[[m]], 2, function(v) rep(v, dim_S[[m]])) }
  # Enforce column-wise sparsity pattern
  for (m in 1:M) { load_J[[m]] <- load_J[[m]] * elem_sprs_J[[m]] }
  for (m in 1:M) { load_S[[m]] <- load_S[[m]] * elem_sprs_S[[m]] }

  # Superimpose unstructured loadings
  for (m in 1:M) {
    load_J[[m]] <- load_J[[m]] + sqrt(nu2_m[m] / 10) * matrix(rnorm(p_m[m] * params$K0), ncol = params$K0)
    load_S[[m]] <- load_S[[m]] + sqrt(pi2_m[m] / 10) * matrix(rnorm(p_m[m] * params$K0_m[m]), ncol = params$K0_m[m])
  }

  # Re-weight loadings columns by index to enforce least eigenvalue problem
  if (params$decaying_loadings) {
    if (params$decaying_type == 'sqrt') {
      for (m in 1:M) {
        load_J[[m]] <- load_J[[m]] * matrix(sqrt(params$K0 / sum(1 / sqrt(1:params$K0))) / sqrt(1:params$K0), nrow = p_m[m], ncol = params$K0, byrow = TRUE)
        load_S[[m]] <- load_S[[m]] * matrix(sqrt(params$K0_m[m] / sum(1 / sqrt(1:params$K0_m[m]))) / sqrt(1:params$K0_m[m]), nrow = p_m[m], ncol = params$K0_m[m], byrow = TRUE)
      }
    } else if (params$decaying_type == 'log') {
      for (m in 1:M) {
        load_J[[m]] <- load_J[[m]] * matrix(sqrt(params$K0 / sum(1 / log(1 + 1:params$K0))) / log(1 + 1:params$K0), nrow = p_m[m], ncol = params$K0, byrow = TRUE)
        load_S[[m]] <- load_S[[m]] * matrix(sqrt(params$K0_m[m] / sum(1 / log(1 + 1:params$K0_m[m]))) / log(1 + 1:params$K0_m[m]), nrow = p_m[m], ncol = params$K0_m[m], byrow = TRUE)
      }
    } else {
      stop("decaying_type must be one of 'sqrt' and 'log'")
    }
  }

  # Sample idiosyncratic noises
  s2m <- lapply(1:M, function(m)
    (1 / rgamma(p_m[m], shape = 1 / 0.1, rate = 1 / 0.1)) * (rowSums(load_J[[m]]^2) + rowSums(load_S[[m]]^2)) / params$snr_m[m]
  )
  # Ensure no noise variance is zero
  for (m in 1:M) {
    if (min(s2m[[m]]) == 0) { s2m[[m]][s2m[[m]] == 0] <- 1 }
  }

  # Close all graphics devices
  graphics.off()

  # Initialize list for storing sorted correlation matrices
  sorted_correlation_matrix <- list()

  # Loop over each view to calculate and plot the correlation matrix
  for (m in 1:M) {
    cor_m <- cov2cor(tcrossprod(load_J[[m]]) + tcrossprod(load_S[[m]]) + diag(s2m[[m]]))
    idxHC_m <- hclust(as.dist(1 - cor_m), method = 'average')
    hHC_m <- idxHC_m$height[p_m[m] - 10]
    cutree_m <- cutree(tree = idxHC_m, h = hHC_m)
    idx_sort_HC_m <- sort(cutree_m, index.return = TRUE)$ix
    plotmat(cor_m[idx_sort_HC_m, idx_sort_HC_m], filename = paste0('cor_', m), is_cor = TRUE, folder = file.path(mainDir, subDir))
    load_J[[m]] <- load_J[[m]][idx_sort_HC_m, ]
    load_S[[m]] <- load_S[[m]][idx_sort_HC_m, ]
    s2m[[m]] <- s2m[[m]][idx_sort_HC_m]
    plotmat(load_J[[m]], filename = paste0('Lambda_', m), folder = file.path(mainDir, subDir))
    plotmat(load_S[[m]], filename = paste0('Gamma_', m), folder = file.path(mainDir, subDir))
    pdf(file = file.path(mainDir, subDir, paste('snr_', m, '.pdf', sep = "")), height = 3, width = 5)
    par(mar = c(5, 5, 1, 1))
    hist((rowSums(load_J[[m]]^2) + rowSums(load_S[[m]]^2)) / s2m[[m]], main = NULL, xlab = 'Signal to Noise')
    dev.off()
  }

  # Optionally plot crossproduct matrices
  if (Crossprd_plot == TRUE) {
    graphics.off()
    for (m in 1:M) {
      plotmat(tcrossprod(load_J[[m]]), filename = paste0('LLt_', m), folder = file.path(mainDir, subDir))
      plotmat(tcrossprod(load_S[[m]]), filename = paste0('GGt_', m), folder = file.path(mainDir, subDir))
    }
  }

  # Optionally save generated data
  if (Data_save == TRUE) {
    eta <- matrix(rnorm(n * params$K0), nrow = n)
    eta_test <- matrix(rnorm(nTest * params$K0), nrow = nTest)
    phi_m <- lapply(params$K0_m, function(kk) matrix(rnorm(n * kk), nrow = n))
    phi_m_test <- lapply(params$K0_m, function(kk) matrix(rnorm(nTest * kk), nrow = nTest))
    Theta <- (2 * rbinom(params$K0, 1, 0.5) - 1) * rbeta(params$K0, shape1 = 5, shape2 = 3) / sqrt(params$K0)
    Theta <- unname(Theta * act_J[nrow(act_J), ]) # Use last row of act_J for Theta
    s2y <- sum(Theta^2) / params$snr_y
    X_m <- lapply(1:M, function(m)
      t(tcrossprod(load_J[[m]], eta) + tcrossprod(load_S[[m]], phi_m[[m]]) +
          sqrt(s2m[[m]]) * matrix(rnorm(p_m[m] * n), ncol = n)))
    X_m_test <- lapply(1:M, function(m)
      t(tcrossprod(load_J[[m]], eta_test) + tcrossprod(load_S[[m]], phi_m_test[[m]]) +
          sqrt(s2m[[m]]) * matrix(rnorm(p_m[m] * nTest), ncol = nTest)))
    yTrain <- eta %*% Theta + sqrt(s2y) * rnorm(n)
    yTest <- eta_test %*% Theta + sqrt(s2y) * rnorm(nTest)
    colnames(yTest) <- c("response")
    colnames(yTrain) <- c("response")
    preprocess_y <- preProcess(yTrain, method = c("center", "scale"))
    yTrain <- as.matrix(predict(preprocess_y, yTrain))
    yTest <- as.matrix(predict(preprocess_y, yTest))
    preprocess_X_m <- list()
    for (m in 1:M) {
      colnames(X_m[[m]]) <- seq(ncol(X_m[[m]]))
      colnames(X_m_test[[m]]) <- seq(ncol(X_m_test[[m]]))
      preprocess_X_m[[m]] <- preProcess(X_m[[m]], method = c("center", "scale"))
      X_m[[m]] <- as.matrix(predict(preprocess_X_m[[m]], X_m[[m]]))
      X_m_test[[m]] <- as.matrix(predict(preprocess_X_m[[m]], X_m_test[[m]]))
    }
    Data <- list(M = M, n = n, nTest = nTest, p_m = p_m, yTrain = yTrain, yTest = yTest,
                 X_m = X_m, X_m_test = X_m_test, preprocess_X_m = preprocess_X_m,
                 preprocess_y = preprocess_y, s2y = s2y, s2m = s2m, Theta = Theta,
                 Lambda_m = load_J, Gamma_m = load_S, eta = eta, eta_test = eta_test,
                 phi_m = phi_m, phi_m_test = phi_m_test)
    saveRDS(Data, file = file.path(mainDir, "Simulated_multiview.rds"))
  }

  # Optionally plot empirical correlations
  if (Emprical_correlation_plot == TRUE) {
    graphics.off()
    for (m in 1:M) {
      plotmat(cor(X_m[[m]]), filename = paste0('cor_emp_', m), is_cor = TRUE, folder = file.path(mainDir, subDir))
    }
  }

  # Optionally plot the activity pattern
  if (Activity_pattern_plot == TRUE) {
    pdf(file = file.path(mainDir, subDir, 'venn_sim_true.pdf'), height = 4, width = 8)
    grid::grid.newpage()

    # Calculate areas for each set
    areas <- sapply(1:(M + 1), function(i) sum(act_J[i, ]))

    # Calculate intersections
    n_intersections <- matrix(0, nrow = M + 1, ncol = choose(M + 1, 2))
    for (i in 1:(M + 1)) {
      combs <- combn(1:(M + 1), i)
      for (j in 1:ncol(combs)) {
        if (length(combs[, j]) == 1) {
          n_intersections[i, j] <- sum(act_J[combs[, j], ])
        } else {
          n_intersections[i, j] <- sum(apply(act_J[combs[, j], , drop = FALSE], 2, prod))
        }
      }
    }

    # Create Venn diagram based on the number of views
    category <- c(expression(Lambda[1]), expression(Lambda[2]), expression(Lambda[3]), expression(theta))[1:(M + 1)]
    create_venn_diagram(M, areas, n_intersections, category, mainDir, subDir)

    dev.off()
  }

  return(Data)
}

