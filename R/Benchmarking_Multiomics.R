#' Simulated Multi-Omics Data Generation
#'
#' This function generates simulated multi-omics data for benchmarking and testing.
#'
#' @param mainDir Character. The main directory for saving plots. Default is '/Users/siddheshkulkarni/Desktop'.
#' @param subDir Character. The subdirectory for saving plots. Default is 'plots_simulated_data'.
#' @param seed.number Integer. The seed number for random number generation. Default is 1234.
#' @param M Integer. Number of views (omics layers). Default is 3.
#' @param p_m Integer vector. Number of features per view. Default is c(500, 700, 1000).
#' @param n Integer. Train set size. Default is 100.
#' @param nTest Integer. Test set size. Default is 100.
#' @param fixed_act_patter Logical. Set a custom activity pattern in the joint component. Default is TRUE.
#' @param J0_m Integer vector. Number of groups in the joint component. Default is c(5, 6, 7).
#' @param K0_m Integer vector. Number of factors in view-specific components. Default is c(5, 7, 10).
#' @param S0_m Integer vector. Number of groups in view-specific components. Default is c(6, 7, 8).
#' @param nu2_0 Numeric. Input the joint loadings variability. Default is 0.1.
#' @param pi2_0 Numeric. Input the view-specific loadings variability. Default is 0.1.
#' @param decaying_loadings Logical. Enforce loadings decay (re-weighting importance of latent axis of variation). Default is FALSE.
#' @param decaying_type Character. Type of decay ('sqrt' or 'log'). Default is 'sqrt'.
#' @param pJ_zero Numeric. Probability of group-wise zero in joint component. Default is 0.6.
#' @param pJ_sign Numeric. Probability of group-wise sign-switch in joint component. Default is 0.4.
#' @param pJ_zeroEntr Numeric. Probability of entry-wise zero in joint component. Default is 0.5.
#' @param pS_zero Numeric. Probability of group-wise zero in specific components. Default is 0.6.
#' @param pS_sign Numeric. Probability of group-wise sign-switch in specific components. Default is 0.4.
#' @param pS_zeroEntr Numeric. Probability of entry-wise zero in specific components. Default is 0.5.
#' @param Activity_pattern_plot Logical. Save the Activity Pattern Plot. Default is TRUE.
#' @param Emprical_correlation_plot Logical. Save the Empirical correlation plot. Default is TRUE.
#' @param Crossprd_plot Logical. Save the cross-product matrices plot. Default is TRUE.
#' @param Data_save Logical. Save the data. Default is TRUE.
#'
#' @return List containing simulated multi-omics data and parameters.
#' @export
### load required packages

# Clear the workspace and close all graphics devices

# rm(list = ls())
#
# graphics.off()
#
# # Load necessary libraries
# library(VennDiagram)  # For drawing Venn diagrams
# library(grid)  # For graphical functions
# library(caret)  # For data pre-processing and model training
# library(fields)  # For image plotting and data manipulation
# library(ggplot2)  # For data visualization

# Custom Functions Definition --------------------------------------------------

# Function to generate block dimensions based on Dirichlet distribution

# getBlocksDim <- function(K, p, alpha = 2) {
#   nK = rep(0, K)
#   while (min(nK) == 0) {
#     pK <- LaplacesDemon::rdirichlet(1, rep(alpha, K))
#     nK = pmax(1, floor(p * pK))
#     nK = c(nK[1:K - 1], p - sum(nK[1:K - 1]))
#   }
#   return(sort(nK, decreasing = TRUE))
# }

# Function to plot a matrix with specified settings
# plotmat <- function(mat_m, n_colors = 256, lab_col = NULL, is_cor = FALSE,
#                     folder = '~/Desktop/', filename = 'prova') {
#
#   if (is_cor) {
#     zL = -1
#     zR = 1
#   } else {
#     zM = max(abs(mat_m))
#     zL = -zM
#     zR = zM
#     mat_m = t(mat_m)
#   }
#
#   file_path <- file.path(folder, paste(filename, '.png', sep = ""))
#   png(file = file_path, height = 5, width = 5.6, res = 300, pointsize = 3.5, unit = 'cm')
#
#   colors <- colorRampPalette(c("#AA4499", "white", "#117733"))(n_colors)
#   ticks <- sapply(c(0.7 * zL, 0, 0.7 * zR), round, 2)
#
#   if (is_cor) {
#     par(pty = "s", mar = c(1, 0, 1, 4))
#     image(mat_m, useRaster = TRUE, asp = 1, axes = FALSE, col = colors,
#           xlim = c(0, 1), ylim = c(0, 1), zlim = c(zL, zR))
#     rect(0, 0, 1, 1, border = "black")
#   } else {
#     par(pty = "s", mar = c(0, 0, 0, 4))
#     image(mat_m, useRaster = TRUE, asp = 1, axes = FALSE, col = colors,
#           xlim = c(0 - 0.1 * 15 / nrow(mat_m), 1 + 0.1 * 15 / nrow(mat_m)), ylim = c(0, 1), zlim = c(zL, zR))
#     rect(0 - 0.04 * 15 / nrow(mat_m), 0, 1 + 0.04 * 15 / nrow(mat_m), 1, border = "black")
#   }
#
#   image.plot(zlim = c(zL, zR), col = colors, legend.only = TRUE, side = 4,
#              axis.args = list(at = ticks, labels = TRUE, cex.axis = 0.9), legend.shrink = 0.5,
#              legend.width = 0.9, legend.mar = 4.5)
#
#   if (!is.null(lab_col)) {
#     mtext(lab_col, side = 4, line = 1, cex = 1.3, las = 1, at = 0.8)
#   }
#
#   dev.off()
# }


Multiview_learner <- function( mainDir = '/Users/siddheshkulkarni/Desktop',
                                     subDir = 'plots_simulated_data',
                                     seed.number = 1234,
                                     M = 3, # Number of views (NB: the inputs are currently set up only for this specific case!!!)
                                     p_m = c(500, 700, 1000),  # Number of features per view
                                     n = 100,  # Train set size
                                     nTest = 100,  # Test set size
                                     fixed_act_pattern = TRUE,  # Set a custom activity pattern in the joint component
                                     K0 = 15,  #number of factors
                                     nPatternJ = 7,
                                     J0_m = c(5, 6, 7),  # Number of groups in joint component
                                     K0_m = c(5, 7, 10),  # Number of factors in view-specific components
                                     S0_m = c(6, 7, 8),  # Number of groups in view-specific components

                                     nu2_0 = 0.1,  # Input the joint loadings variability
                                     pi2_0 = 0.1,  # Input the view-specific loadings variability

                                     snr_y = 1,  # Input signal-to-noise ratio

                                     decaying_loadings = FALSE,  # Enforce loadings decay (re-weighting importance of latent axis of variation)

                                     decaying_type = 'sqrt',  # One of 'sqrt' and 'log'

                                     pJ_zero = 0.6,  # Probability of group-wise zero in joint component
                                     pJ_sign = 0.4,  # Probability of group-wise sign-switch in joint component
                                     pJ_zeroEntr = 0.5,  # Probability of entry-wise zero in joint component

                                     pS_zero = 0.6,  # Probability of group-wise zero in specific components
                                     pS_sign = 0.4,  # Probability of group-wise sign-switch in specific components
                                     pS_zeroEntr = 0.5, # Probability of entry-wise zero in specific components

                                     Activity_pattern_plot = TRUE,  #Save the Activity Pattern Plot

                                     Emprical_correlation_plot = TRUE, #Save the Emprical correlation plot

                                     Crossprd_plot=TRUE, #Save crossporduct plots

                                     Data_save=TRUE  ## save the data in directory
){
# Inputs & Dataset Dimensions --------------------------------------------------

# Create directory for saving plots

#mainDir = '~/Desktop'
#subDir = 'plots_simulated_data'

# mainDir = '/Users/siddheshkulkarni/Desktop'
# subDir = 'plots_simulated_data'
#
#

dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

set.seed(seed.number)  # Set random seed for reproducibility

# M = 3  # Number of views (NB: the inputs are currently set up only for this specific case!!!)
# p_m = c(500, 700, 1000)  # Number of features per view
# n = 100  # Train set size
# nTest = 100  # Test set size

# fixed_act_patter = TRUE  # Set a custom activity pattern in the joint component

if (fixed_act_pattern==TRUE & M == 3) {

  # Input desired non-zero number of active columns in selected configurations

  which_comb <- c('1110', '1010', '1001', '0111', '1101', '1011')

  value_comb <- c(2, 3, 1, 2, 3, 4)

  K0 <- sum(value_comb)  # Number of factors in joint component
} else {
  # Generate encodings for all possible activity configurations

  which_comb <- expand.grid(replicate(M + 1, c(0, 1), simplify = FALSE))

  which_comb <- apply(which_comb[rowSums(which_comb) > 1, ], 1, paste, collapse = "")

  # Select number of factors and number of appearing
  #K0 = 15
  #nPatternJ = 7
  #nPatternJ = min(nPatternJ, length(which_comb))

  # Sample appearing configurations and dedicated number of factors

  which_comb <- sample(which_comb, nPatternJ)

  value_comb <- getBlocksDim(nPatternJ, K0, alpha = 5)
}

# #### Make changes in the code from here ##
#
# J0_m <- c(5, 6, 7)  # Number of groups in joint component
# K0_m <- c(5, 7, 10)  # Number of factors in view-specific components
# S0_m <- c(6, 7, 8)  # Number of groups in view-specific components
#
# nu2_0 <- 0.1  # Input the joint loadings variability
# pi2_0 <- 0.1  # Input the view-specific loadings variability
#

snr_m <- rep(1, M)  # Input the view-specific signal-to-noise ratio (mean)

#
# decaying_loadings = FALSE  # Enforce loadings decay (re-weighting importance of latent axis of variation)
# decaying_type = 'sqrt'  # One of 'sqrt' and 'log'
#
# pJ_zero <- 0.6  # Probability of group-wise zero in joint component
# pJ_sign <- 0.4  # Probability of group-wise sign-switch in joint component
# pJ_zeroEntr <- 0.5  # Probability of entry-wise zero in joint component
#
# pS_zero <- 0.6  # Probability of group-wise zero in specific components
# pS_sign <- 0.4  # Probability of group-wise sign-switch in specific components
# pS_zeroEntr <- 0.5  # Probability of entry-wise zero in specific components

# Generate Structured Loadings -----

## Set Sparsity Pattern in Joint Component -------------------------------------

# Sample the group sizes. ## what is meant by group sizes?

dim_J <- lapply(1:M, function(m) getBlocksDim(J0_m[m], p_m[m]))

dim_S <- lapply(1:M, function(m) getBlocksDim(S0_m[m], p_m[m]))

# Structured activity pattern in joint component. ## what is meant by group sizes?

act_J <- sapply(which_comb, function(v) as.numeric(unlist(strsplit(v, ""))))

act_J <- t(apply(act_J, 1, function(v) rep(v, value_comb)))

# Group-wise sparsity pattern & sign switching

group_sprs_J <- list()

m<-1

for (m in 1:M){

  group_sprs_J[[m]] <-    matrix(rbinom(J0_m[m] * K0, 1, 1 - pJ_zero * (sum(act_J[m,]) / K0)), ncol = K0) *
                          (2 * matrix(rbinom(J0_m[m] * K0, 1, 1 - pJ_sign), ncol = K0) - 1)

  # Check that no active column is zeros only

  col_wrong <-  apply(group_sprs_J[[m]], 2, function(v)
                min(v == 0)) * act_J[m,]

  while (sum(col_wrong) > 0) {

    group_sprs_J[[m]][, as.logical(col_wrong)] <- matrix(rbinom(J0_m[m] * sum(col_wrong), 1, 1 - pJ_zero * (sum(act_J[m, ]) / K0)), ncol = sum(col_wrong)) *
                                                  (2 * matrix(rbinom(J0_m[m] * sum(col_wrong), 1, 1 - pJ_sign), ncol = sum(col_wrong)) - 1)

    col_wrong <- apply(group_sprs_J[[m]], 2, function(v) min(v == 0)) * act_J[m, ]
  }

  # Check that no group is zeros only
  grp_wrong <- apply(group_sprs_J[[m]] * (rep(1, J0_m[m]) %x% t(act_J[m,])), 1, function(v) min(v == 0))

  while (sum(grp_wrong) > 0) {
    # Replace the zero-only groups with new values
    group_sprs_J[[m]][grp_wrong, ] <- act_J[m,] * (2 * matrix(rbinom(sum(grp_wrong) * K0, 1, 1 - pJ_sign), ncol = K0) - 1)
  }

}
  ## Column-wise sparsity pattern (entry-wise)

  elem_sprs_J <- lapply(1:M, function(m) matrix(rbinom(p_m[m] * K0, 1, 1 - pJ_zeroEntr), ncol = K0))

  # Check that no row is zeros only
  for (m in 1:M) {
    sprs_tot <- elem_sprs_J[[m]] * apply(group_sprs_J[[m]], 2, function(v) rep(v, dim_J[[m]]))
    if (max(apply(sprs_tot, 1, function(v) min(v == 0))) > 0) {
      # Replace the zero-only rows with new values
      elem_sprs_J[[m]][which(apply(sprs_tot, 1, function(v) min(v == 0)) > 0), act_J[m, ]] <- 1
    }
  }


  ## Set group-wise sparsity pattern in view-specific components

  ## Group-wise sparsity pattern & sign switching
  group_sprs_S <- list()

  for (m in 1:M) {
    group_sprs_S[[m]] <- matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_zero), ncol = K0_m[m]) *
                        (2 * matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_sign), ncol = K0_m[m]) - 1)

    # Check that no group is zeros only

    while (max(apply(group_sprs_S[[m]], 1, function(v) min(v == 0))) > 0) {
      # Replace the zero-only groups with new values
      group_sprs_S[[m]] <- matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_zero), ncol = K0_m[m]) *
        (2 * matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_sign), ncol = K0_m[m]) - 1)
    }
  }

  ## Column-wise sparsity pattern (entry-wise)
  elem_sprs_S <- lapply(1:M, function(m) matrix(rbinom(p_m[m] * K0_m[m], 1, 1 - pS_zeroEntr), ncol = K0_m[m]))
  # Check that no row is zeros only
  for (m in 1:M) {
    sprs_tot <- elem_sprs_S[[m]] * apply(group_sprs_S[[m]], 2, function(v) rep(v, dim_S[[m]]))
    if (max(apply(sprs_tot, 1, function(v) min(v == 0))) > 0) {
      # Replace the zero-only rows with new values
      elem_sprs_S[[m]][which(apply(sprs_tot, 1, function(v) min(v == 0)) > 0), ] <- 1
    }
  }

  ## Generate the true loading matrices

  # Set the means for each loading group

  mean_J <- lapply(1:M, function(m) (-1) ^ c(1:J0_m[m]) * rbeta(J0_m[m], shape1 = 5, shape2 = 3) / sqrt(K0))

  mean_S <- lapply(1:M, function(m) (-1) ^ c(1:S0_m[m]) * rbeta(S0_m[m], shape1 = 5, shape2 = 3) / sqrt(K0_m[m]))

  # Set the view-specific variability

  nu2_m <- nu2_0 / rep(K0, M)

  pi2_m <- nu2_0 / K0_m

  # Assign the feature to a group

  col_J <- sapply(1:M, function(m) unlist(sapply(1:J0_m[m], function(l) rep(l, dim_J[[m]][l]))))

  col_S <- sapply(1:M, function(m) unlist(sapply(1:S0_m[m], function(l) rep(l, dim_S[[m]][l]))))

  # Sample the full loadings matrices

  load_J <- lapply(1:M, function(m)
    matrix(mean_J[[m]][col_J[[m]]], nrow = p_m[m], ncol = K0) +
      sqrt(nu2_m[m]) * matrix(rnorm(p_m[m] * K0), nrow = p_m[m], ncol = K0))

  load_S <- lapply(1:M, function(m)
    matrix(mean_S[[m]][col_S[[m]]], nrow = p_m[m], ncol = K0_m[m]) +
      sqrt(pi2_m[m]) * matrix(rnorm(p_m[m] * K0_m[m]), nrow = p_m[m], ncol = K0_m[m])
  )

  # Enforce the desired activity in the joint loadings
  for (m in 1:M) { load_J[[m]][, !act_J[m,]] <- 0 }

  # Enforce group-wise sparsity pattern
  for (m in 1:M) { load_J[[m]] <- load_J[[m]] * apply(group_sprs_J[[m]], 2, function(v) rep(v, dim_J[[m]])) }
  for (m in 1:M) { load_S[[m]] <- load_S[[m]] * apply(group_sprs_S[[m]], 2, function(v) rep(v, dim_S[[m]])) }

  # Enforce column-wise sparsity pattern
  for (m in 1:M) { load_J[[m]] <- load_J[[m]] * elem_sprs_J[[m]] }
  for (m in 1:M) { load_S[[m]] <- load_S[[m]] * elem_sprs_S[[m]] }

  # Superimpose unstructured loadings
  for (m in 1:M) {
    load_J[[m]] <- load_J[[m]] + sqrt(nu2_m[m] / 10) * matrix(rnorm(p_m[m] * K0), ncol = K0)
    load_S[[m]] <- load_S[[m]] + sqrt(pi2_m[m] / 10) * matrix(rnorm(p_m[m] * K0_m[m]), ncol = K0_m[m])
  }

  # Re-weight loadings columns by index (to enforce least eigenvalue problem)

  if (decaying_loadings==TRUE) {
    if (decaying_type == 'sqrt') {
      for (m in 1:M) {
        load_J[[m]] <- load_J[[m]] * matrix(sqrt(K0 / sum(1 / sqrt(1:K0))) / sqrt(1:K0), nrow = p_m[m], ncol = K0, byrow = TRUE)
        load_S[[m]] <- load_S[[m]] * matrix(sqrt(K0_m[m] / sum(1 / sqrt(1:K0_m[m]))) / sqrt(1:K0_m[m]), nrow = p_m[m], ncol = K0_m[m], byrow = TRUE)
      }
    } else if (decaying_type == 'log') {
      for (m in 1:M) {
        load_J[[m]] <- load_J[[m]] * matrix(sqrt(K0 / sum(1 / log(1 + 1:K0))) / log(1 + 1:K0), nrow = p_m[m], ncol = K0, byrow = TRUE)
        load_S[[m]] <- load_S[[m]] * matrix(sqrt(K0_m[m] / sum(1 / log(1 + 1:K0_m[m]))) / log(1 + 1:K0_m[m]), nrow = p_m[m], ncol = K0_m[m], byrow = TRUE)
      }
    } else {
      stop("decaying_type must be one of 'sqrt' and 'log'")
    }
  }

  # Sample idiosyncratic noises

  s2m <- lapply(1:M, function(m)
    (1 / rgamma(p_m[m], shape = 1 / 0.1, rate = 1 / 0.1)) * (rowSums(load_J[[m]]^2) + rowSums(load_S[[m]]^2)) / snr_m[m]
  )
  # Check that no one is zero
  for (m in 1:M) {
    if (min(s2m[[m]]) == 0) {
      s2m[[m]][s2m[[m]] == 0] <- 1
    }
  }

  # Inspect Results

  # Close all graphics devices
  graphics.off()

  sorted_correlation_matrix <- list()


  # Loop over each view (omics layer)
  for (m in 1:M) {

    # Calculate the correlation matrix for the view
    cor_m <- cov2cor(tcrossprod(load_J[[m]]) + tcrossprod(load_S[[m]]) + diag(s2m[[m]]))

    # Perform hierarchical clustering on the distance matrix (1 - correlation)

    idxHC_m  <- hclust(as.dist(1 - cor_m), method = 'average')

    hHC_m    <- idxHC_m$height[p_m[m] - 10]

    cutree_m <- cutree(tree = idxHC_m, h = hHC_m)

    # Sort indices based on hierarchical clustering
    idx_sort_HC_m <- sort(cutree_m, index.return = TRUE)$ix

    # Plot the sorted correlation matrix

    plotmat(cor_m[idx_sort_HC_m, idx_sort_HC_m], filename = paste0('cor_', m), is_cor = TRUE, folder = file.path(mainDir, subDir))

    # Sort loadings and variances based on hierarchical clustering

    load_J[[m]] <- load_J[[m]][idx_sort_HC_m, ]

    load_S[[m]] <- load_S[[m]][idx_sort_HC_m, ]

    s2m[[m]] <- s2m[[m]][idx_sort_HC_m]

    # Plot the sorted joint and specific loadings matrices

    plotmat(load_J[[m]], filename = paste0('Lambda_', m), folder = file.path(mainDir, subDir))

    plotmat(load_S[[m]], filename = paste0('Gamma_', m), folder = file.path(mainDir, subDir))

    # Plot the signal-to-noise ratio histogram
    file_path <- pdf(file = file.path(mainDir, subDir, paste('snr_', m, '.pdf', sep = "")), height = 3, width = 5)

    par(mar = c(5, 5, 1, 1))

    hist((rowSums(load_J[[m]]^2) + rowSums(load_S[[m]]^2)) / s2m[[m]], main = NULL, xlab = 'Signal to Noise')

    dev.off()
  }


  # Optionally, plot additional matrices if needed

  if (Crossprd_plot==TRUE) {

    graphics.off()

    for (m in 1:M) {

      plotmat(tcrossprod(load_J[[m]]), filename = paste0('LLt_', m), folder = file.path(mainDir, subDir))

      plotmat(tcrossprod(load_S[[m]]), filename = paste0('GGt_', m), folder = file.path(mainDir, subDir))
    }
  }

  if (Data_save==TRUE) {

    # Generate Data --------------------------------------------------------------

    ## Generate true latent factors

    # Joint factors

    eta <- matrix(rnorm(n * K0), nrow = n)

    eta_test <- matrix(rnorm(nTest * K0), nrow = nTest)

    # View-specific factors

    phi_m <- lapply(K0_m, function(kk) matrix(rnorm(n * kk), nrow = n))

    phi_m_test <- lapply(K0_m, function(kk) matrix(rnorm(nTest * kk), nrow = nTest))

    ## Generate response coefficients and noise

    # Response coefficients

    Theta <- (2 * rbinom(K0, 1, 0.5) - 1) * rbeta(K0, shape1 = 5, shape2 = 3) / sqrt(K0)

    Theta <- unname(Theta * act_J[M + 1, ])

    # Response noise variance

    s2y <- sum(Theta^2) / snr_y

    ## Generate predictors

    # Train predictors

    X_m <- lapply(1:M, function(m)
                 t(tcrossprod(load_J[[m]], eta) + tcrossprod(load_S[[m]], phi_m[[m]]) +
                 sqrt(s2m[[m]]) * matrix(rnorm(p_m[m] * n), ncol = n)))

    # Test predictors

    X_m_test <- lapply(1:M, function(m)
                 t(tcrossprod(load_J[[m]], eta_test) + tcrossprod(load_S[[m]], phi_m_test[[m]]) +
                sqrt(s2m[[m]]) * matrix(rnorm(p_m[m] * nTest), ncol = nTest)))

    ## Generate responses

    yTrain <- eta %*% Theta + sqrt(s2y) * rnorm(n)

    yTest  <- eta_test %*% Theta + sqrt(s2y) * rnorm(nTest)

    # Rescale Data ---------------------------------------------------------------

    colnames(yTest)  <- c("response")
    colnames(yTrain) <- c("response")

    # Center and scale responses

    preprocess_y <- preProcess(yTrain, method = c("center", "scale"))

    yTrain <- as.matrix(predict(preprocess_y, yTrain))

    yTest <- as.matrix(predict(preprocess_y, yTest))

    preprocess_X_m <- list()

    for (m in 1:M) {

      colnames(X_m[[m]])  <-  seq(ncol(X_m[[m]]))
      colnames(X_m_test[[m]]) <-  seq(ncol(X_m_test[[m]]))

      # Center and scale omics data

      preprocess_X_m[[m]] <-  preProcess(X_m[[m]], method = c("center", "scale"))
      X_m[[m]] <-  as.matrix(predict(preprocess_X_m[[m]], X_m[[m]]))
      X_m_test[[m]] <- as.matrix(predict(preprocess_X_m[[m]], X_m_test[[m]]))

    }
  }

## Return Data Objects ###

  if (Data_save==TRUE) {
    # Save Data ------------------------------------------------------------------

    Data <- list(M = M, n = n, nTest = nTest, p_m = p_m, yTrain = yTrain, yTest = yTest,
                 X_m = X_m, X_m_test = X_m_test, preprocess_X_m = preprocess_X_m,
                 preprocess_y = preprocess_y, s2y = s2y, s2m = s2m, Theta = Theta,
                 Lambda_m = load_J, Gamma_m = load_S, eta = eta, eta_test = eta_test,
                 phi_m = phi_m, phi_m_test = phi_m_test)

    saveRDS(Data, file = file.path(mainDir, "Simulated_multiview.rds"))
  }

  # Plot empirical correlations --------------------------------------------------

  if (Emprical_correlation_plot==TRUE) {
    graphics.off()
    for (m in 1:M) {
      plotmat(cor(X_m[[m]]), filename = paste0('cor_emp_', m), is_cor = TRUE, folder = file.path(mainDir, subDir))
    }
  }

  # Visualizing enforced activity pattern ----------------------------------------

  if (Activity_pattern_plot==TRUE & M == 3) {
    # Visualization of resulting activity pattern

    pdf(file = file.path(mainDir, subDir, 'venn_sim_true.pdf'), height = 4, width = 8)

    grid::grid.newpage()

    prova <- draw.quad.venn(area1 = sum(act_J[1, ]),
                            area2 = sum(act_J[2, ]),
                            area3 = sum(act_J[3, ]),
                            area4 = sum(act_J[4, ]),
                            # Pairwise
                            n12 = sum(apply(act_J[c(1, 2), ], 2, prod)),
                            n23 = sum(apply(act_J[c(2, 3), ], 2, prod)),
                            n13 = sum(apply(act_J[c(1, 3), ], 2, prod)),
                            n14 = sum(apply(act_J[c(1, 4), ], 2, prod)),
                            n24 = sum(apply(act_J[c(2, 4), ], 2, prod)),
                            n34 = sum(apply(act_J[c(3, 4), ], 2, prod)),
                            # Triple
                            n123 = sum(apply(act_J[c(1, 2, 3), ], 2, prod)),
                            n124 = sum(apply(act_J[c(1, 2, 4), ], 2, prod)),
                            n234 = sum(apply(act_J[c(2, 3, 4), ], 2, prod)),
                            n134 = sum(apply(act_J[c(1, 3, 4), ], 2, prod)),
                            # Quadruple
                            n1234 = sum(apply(act_J, 2, prod)),
                            category = c(expression(Lambda[1]),
                                         expression(Lambda[2]),
                                         expression(Lambda[3]),
                                         expression(theta)),
                            col = "black", fill = c('#332288', '#882255', '#CC6677', '#117733'),
                            alpha = c(0.4, 0.4, 0.35, 0.35), scaled = TRUE,
                            lty = 'blank', cex = 1.5, cat.cex = 1.6)

    dev.off()
  }


return(Data)

}


#   ### Input for Functions
#
# Benchmarking_output <- Benchmarking_multiomics(
#
#   mainDir = '/Users/siddheshkulkarni/Desktop',
#
#   subDir = 'plots_simulated_data',
#
#   seed.number = 1234,
#
#   M = 3, # Number of views (NB: the inputs are currently set up only for this specific case!!!)
#   p_m = c(500, 700, 1000),  # Number of features per view
#   n = 100,  # Train set size
#   nTest = 100,  # Test set size
#   fixed_act_patter = TRUE,  # Set a custom activity pattern in the joint component
#
#   J0_m = c(5, 6, 7),  # Number of groups in joint component
#   K0_m = c(5, 7, 10),  # Number of factors in view-specific components
#   S0_m = c(6, 7, 8),  # Number of groups in view-specific components
#
#   nu2_0 = 0.1,  # Input the joint loadings variability
#   pi2_0 = 0.1,  # Input the view-specific loadings variability
#
#   decaying_loadings = FALSE,  # Enforce loadings decay (re-weighting importance of latent axis of variation)
#   decaying_type = 'sqrt',  # One of 'sqrt' and 'log'
#
#   pJ_zero = 0.6,  # Probability of group-wise zero in joint component
#   pJ_sign = 0.4,  # Probability of group-wise sign-switch in joint component
#   pJ_zeroEntr = 0.5,  # Probability of entry-wise zero in joint component
#
#   pS_zero = 0.6,  # Probability of group-wise zero in specific components
#   pS_sign = 0.4,  # Probability of group-wise sign-switch in specific components
#   pS_zeroEntr = 0.5, # Probability of entry-wise zero in specific components
#
#   ### Additional Matrix Plots ##
#
#   Activity_pattern_plot = TRUE,  #Save the Activity Pattern Plot
#
#   Emprical_correlation_plot = TRUE, #Save the Empirical correlation plot
#
#   Crossprd_plot=TRUE, ## crossporduct matrices
#
#   Data_save=TRUE  ## save the data in directory
#
# )
#


