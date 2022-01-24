TDC = function(decoy_wins, target_wins, BC1 = 1, c = 0.5, lambda = 0.5){
  nTD <- cumsum(target_wins)
  nDD <- cumsum(decoy_wins)
  fdps <- pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))
  qvals <- rev(cummin(rev(fdps)))
  return(qvals)
}

winning_scores <- c(1:1000)
labels_vec <- matrix(0, nrow = 10, ncol = 1000)
decoy_wins_vec <- matrix(0, nrow = 10, ncol = 1000)
target_wins_vec <- matrix(0, nrow = 10, ncol = 1000)
all_group_ids <- rep(1, 1000)
for (i in 1:10){
  labels_vec[i, ] <- sample(c(-1, 1), size = 1000, replace = T)
  decoy_wins_vec[i, ] <- (labels_vec[i, ] != 1)
  target_wins_vec[i, ] <- (labels_vec[i, ] == 1)
}

test_that("Group_walk agrees with TDC", {
  for (i in 1:10){
    expect_equal(group_walk(winning_scores, labels_vec[i, ], all_group_ids), rev(TDC(rev(decoy_wins_vec[i, ]), rev(target_wins_vec[i, ]))))
  }
})
