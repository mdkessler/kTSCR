C <- 100  # represents samples
R <- 200 # represents features
y <- rnorm(C) # represents outcome variable
X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
Is <- get_pairwise_rank_matrices(X)
pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
View(sorted_corrs)
feature_sum_corr_scores <- calc_feature_sum_corr_scores(sorted_corrs)
View(feature_sum_corr_scores)
names(feature_sum_corr_scores)

pairs <- rbind(c(1,2),c(4,6), c(10,5))
pairs
