#' Implements group-walk algorithm
#'
#' This function returns a list of q-values corresponding to hypotheses that have been partitioned into groups.
#' For FDR control, users should report the target-hypotheses with q-values less than or equal to their choice of threshold, alpha.
#' For further details about how group-walk works, see: https://www.biorxiv.org/content/10.1101/2022.01.30.478144v1
#'
#' @param winning_scores A numerical vector of winning scores generated from the target-decoy competitions for each hypothesis.
#' @param labels A vector of winning labels indicating whether it was a target (= 1) or a decoy (!= 1) for each hypothesis.
#' @param all_group_ids A vector of group IDs associated to each hypothesis (can be recorded as integers, factors, characters).
#' @param K A window size parameter (integer).
#' @param return_frontier A boolean indicating whether the function should return the complete sequence of frontiers.
#' @param correction A correction factor used to in the numerator of the estimated false discovery rate (FDR) (Use 1 for FDR control).
#' @export
#' @examples
#' create_uncalibrated_hypotheses <- function(m_vec, pi_0_vec, mus, sds) {
#'   total <- sum(m_vec)
#'   g_total <- length(m_vec)
#'   data <- matrix(0, ncol = 4, nrow = total)
#'   for (g in 1:length(m_vec)){
#'     m <- m_vec[g]
#'     pi_0 <- pi_0_vec[g]
#'     mu <- mus[g]
#'     sd <- sds[g]
#'     if (g == 1) {
#'       start <- 0
#'     } else {
#'       start <- sum(m_vec[1:(g - 1)])
#'     }
#'     targets_nonnull <- rnorm(floor(m*pi_0), mean = mu, sd = sd)
#'     targets_null <- rnorm(m - floor(m*pi_0), mean = 0, sd = 1)
#'     decoys <- rnorm(m, mean = 0, sd = 1)
#'     targets <- c(targets_nonnull, targets_null)
#'     W <- pmax(targets, decoys)
#'     data[(start + 1):(start + m), 1] <- W
#'     data[(start + 1):(start + m), 2] <- g
#'     decoy_inds <- which(decoys > targets)
#'     inc_native_inds <- (which(targets_null > decoys[(floor(m*pi_0) + 1):m])) + floor(m*pi_0)
#'     X <- rep(0, m)
#'     X[decoy_inds] <- -1
#'     X[inc_native_inds] <- 1
#'     Y <- X
#'     X[X == 0] <- 1
#'     data[(start + 1):(start + m), 3] <- Y
#'     data[(start + 1):(start + m), 4] <- X
#'   }
#'   return(data)
#' }
#'
#' data <- create_uncalibrated_hypotheses(m_vec = rep(1000, 3),
#'            pi_0_vec = rep(0.6, 3), mus = c(2.5, 3, 3.5), sds = rep(1, 3))
#'
#' winning_scores <- data[, 1]
#' all_group_ids <- data[, 2]
#' labels <- data[, 4]
#' q_vals <- group_walk(winning_scores, labels, all_group_ids)
#'
#' @return A sequence of q-values for each hypothesis. If return_frontier = T, additionally the sequence of frontiers will be returned.
group_walk <- function(winning_scores, labels, all_group_ids, K = 40, return_frontier = FALSE, correction = 1) {
  lengths_all <- lengths(list(winning_scores, labels, all_group_ids))
  if (length(unique(lengths_all)) > 1) {
    stop("winning_scores, labels, all_group_ids are not of the same length")
  }

  ordered_inds <- order(winning_scores) # record the indices in order of the scores from smallest to largest
  groups_and_labels <- data.frame(all_group_ids[ordered_inds], labels[ordered_inds]) # reorder the group_ids and winning labels according to the ordered score
  colnames(groups_and_labels) <- c("ordered_groups", "ordered_labels")

  q_vals <- rep(1, length(ordered_inds)) # initialize the vector of q-values
  q_val <- 1
  group_ids <- unique(sort(all_group_ids)) # distinct groups
  g_total <- length(group_ids) # no. of distinct groups
  totals <- tapply(winning_scores, all_group_ids, length) # the number of hypotheses in each group

  labels_sorted <- vector(mode = "list", length = g_total) # initialize the labels broken up into groups (as booleans where TRUE corresponds to a target win)
  for (g in 1:g_total) {
    labels_sorted[[g]] <- groups_and_labels$ordered_labels[groups_and_labels$ordered_groups == group_ids[g]] == 1
  }

  starts <- rep(1, g_total) # frontier starting from the worst score and proceeding to higher scores
  weights <- rep(0, g_total) # weights for each group to determine which group we should step into next

  decoys_plus_one <- sum(!unlist(labels_sorted)) + correction # initialize the number of decoys (plus one) yet to have been seen
  rejections <- sum(unlist(labels_sorted)) # initialize the number of targets yet to have been seen

  switch <- TRUE

  frontiers <- matrix(1, nrow = sum(totals), ncol = 2)
  counter <- 1

  while (any(starts <= totals)) {
    if (return_frontier == TRUE) {
      frontiers[counter, ] <- starts
      counter <- counter + 1
    }

    FDR_e <- decoys_plus_one / max(rejections, 1) # calculates the estimated FDR for scores above the frontier

    q_val <- min(q_val, FDR_e)

    if (any(starts[starts <= totals] <= K)) { # frontier in each group need to exceed K, otherwise we randomly select a group to step forward in.
      inds <- which(starts[starts <= totals] == min(starts[starts <= totals]))
    } else { # calculates the number of decoys in the last K hypotheses within each group

      if (switch) {
        for (g in 1:g_total) {
          if (starts[g] <= totals[g]) {
            weights[g] <- sum(!labels_sorted[[g]][(starts[g] - K):(starts[g] - 1)])
          }
        }
        switch <- FALSE
      } else {
        # weight of the current group drops out the (K + 1)th hypothesis and adds the most recent hypothesis
        weights[index_update] <- weights[index_update] - (!labels_sorted[[index_update]][(starts[index_update] - K - 1)]) +
          (!labels_sorted[[index_update]][(starts[index_update] - 1)])
      }

      inds <- which(weights[starts <= totals] == max(weights[starts <= totals])) # which groups have the most decoys
    }
    if (length(inds) == 1) randind <- inds else randind <- sample(inds, size = 1) # selects the index at random of non-exhausted groups that have the most decoys
    index_update <- which(starts <= totals)[randind] # the index of the corresponding group that needs to be updated
    q_vals[groups_and_labels$ordered_groups == group_ids[index_update]][starts[index_update]] <- q_val # stores the q_val in the corresponding position in order of the winning scores
    label_at_update <- labels_sorted[[index_update]][starts[index_update]] # label associated to the group that needs to be updated

    decoys_plus_one <- decoys_plus_one - !label_at_update # updating the number of decoys
    rejections <- rejections - label_at_update # updating the number of targets

    starts[index_update] <- starts[index_update] + 1
  }
  q_vals[ordered_inds] <- q_vals # reorder the q_vals according to the original order of the winning scores

  if (return_frontier) {
    results <- list(q_vals, frontiers)
    names(results) <- c("q_vals", "frontiers")
    return(results)
  } else {
    return(as.numeric(q_vals))
  }
}
