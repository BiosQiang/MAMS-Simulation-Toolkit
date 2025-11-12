#########################################################################################################
# Function: MAMS_SIM
# Purpose : Perform a single simulation run for a multi-arm multi-stage (MAMS) trial
#
# Author  : Qiang Zhang (Optimised and documented version)
# Date    : 12 Sep 2023
# Version : 1.2
#
# Description:
#   This function simulates one realisation of a MAMS design trial with multiple
#   experimental arms versus a shared control. At each interim analysis, each arm
#   can stop for futility or efficacy according to given boundaries. If eff_stop = TRUE,
#   declaring efficacy in any arm terminates the whole trial.
#
# Inputs:
#   arm         : Number of experimental arms (excluding control)
#   stage_num   : Number of stages (including final)
#   cum_sam     : Vector of cumulative sample sizes per arm at each stage
#   ratio       : Control-to-treatment sample size ratio
#   delta1      : Vector of true means for each arm (length = arm + 1, last is control)
#   sd          : Common within-group standard deviation
#   futi_bound  : Futility boundaries (Z-scale) per stage
#   eff_bound   : Efficacy boundaries (Z-scale) per stage
#   eff_stop    : Logical; if TRUE, efficacy in any arm stops the trial
#   seed        : Optional random seed for reproducibility
#
# Outputs (list):
#   global                : “Rej global null” or “Acc global null”
#   rej_detail            : Detailed rejection statement (e.g., Rej H01+H02)
#   trial_overview         : Short summary of stopping stage and global decision
#   trial_overview_det     : Detailed version with individual rejections
#   final_act              : Data frame of each arm’s decision and stage
#   sta_matrix             : Z-statistics by stage and arm
#   judge_matrix           : Decisions (1=futility,2=efficacy,0=continue)
#   stop_matrix            : Matrix showing at which stage each arm stopped
#   enrolled_sample_size   : Total expected enrolled sample size in this simulation
#########################################################################################################

MAMS_SIM <- function(arm,
                     stage_num,
                     cum_sam,
                     ratio,
                     delta1,
                     sd,
                     futi_bound,
                     eff_bound,
                     eff_stop = TRUE,
                     seed = NULL) {
  
  # ==================== 0. Input validation ==================== #
  stopifnot(length(cum_sam) == stage_num)
  stopifnot(all(diff(cum_sam) > 0))
  stopifnot(length(delta1) == arm + 1)
  stopifnot(length(futi_bound) == stage_num,
            length(eff_bound)  == stage_num)
  if (!is.null(seed)) set.seed(seed)
  
  # ==================== 1. Data generation ==================== #
  n_per_arm_max <- max(cum_sam)
  n_ctrl_max    <- n_per_arm_max * ratio
  
  # Each column = one treatment arm
  trt_mat <- matrix(NA_real_, nrow = n_per_arm_max, ncol = arm)
  for (k in seq_len(arm)) {
    trt_mat[, k] <- rnorm(n_per_arm_max, mean = delta1[k], sd = sd)
  }
  
  # Control arm data
  ctrl_vec <- rnorm(n_ctrl_max, mean = delta1[arm + 1], sd = sd)
  
  # ==================== 2. Compute stage-wise Z statistics ==================== #
  # Z = (mean_treat - mean_ctrl) / sqrt(sd^2*(1/nT + 1/nC))
  Z_mat <- matrix(NA_real_, nrow = stage_num, ncol = arm)
  for (i in seq_len(stage_num)) {
    nT <- cum_sam[i]
    nC <- ratio * nT
    trt_means <- colMeans(trt_mat[seq_len(nT), , drop = FALSE])
    ctrl_mean <- mean(ctrl_vec[seq_len(nC)])
    se        <- sd * sqrt(1 / nT + 1 / nC)
    Z_mat[i, ] <- (trt_means - ctrl_mean) / se
  }
  
  # ==================== 3. Decision classification ==================== #
  # judge_matrix: 1=futility, 2=efficacy, 0=continue
  judge_mat <- matrix(0L, nrow = stage_num, ncol = arm)
  for (i in seq_len(stage_num)) {
    judge_mat[i, ] <- ifelse(Z_mat[i, ] <= futi_bound[i], 1L,
                             ifelse(Z_mat[i, ] > eff_bound[i], 2L, 0L))
  }
  
  # ==================== 4. Identify first action per arm ==================== #
  first_hit <- function(v) {
    idx <- which(v %in% c(1L, 2L))
    if (length(idx) == 0) Inf else idx[1]
  }
  fir_vec <- apply(judge_mat, 2, first_hit)
  
  # Build summary table for each arm
  final_act <- vector("list", arm)
  for (k in seq_len(arm)) {
    if (is.finite(fir_vec[k]) && fir_vec[k] < stage_num) {
      # Early stop (interim)
      if (judge_mat[fir_vec[k], k] == 1L) {
        dec <- sprintf("Drop arm %d at IA%d", k, fir_vec[k])
      } else {
        dec <- sprintf("Rej H0%d at IA%d", k, fir_vec[k])
      }
      if (fir_vec[k] + 1 <= stage_num)
        judge_mat[(fir_vec[k] + 1):stage_num, k] <- 0L
    } else if (is.finite(fir_vec[k]) && fir_vec[k] == stage_num) {
      # Final analysis stop
      if (judge_mat[fir_vec[k], k] == 1L) {
        dec <- sprintf("Acc H0%d at FA", k)
      } else {
        dec <- sprintf("Rej H0%d at FA", k)
      }
    } else {
      # Never triggered any boundary -> accept at final
      fir_vec[k] <- stage_num
      dec <- sprintf("Acc H0%d at FA", k)
    }
    final_act[[k]] <- data.frame(arm = k,
                                 dec = dec,
                                 fir_act = fir_vec[k],
                                 stringsAsFactors = FALSE)
  }
  final_act <- do.call(rbind, final_act)
  
  # ==================== 5. Trial-level stopping logic ==================== #
  stop_matrix <- matrix(0L, nrow = stage_num, ncol = arm + 1)
  
  if (eff_stop) {
    # If any efficacy found, stop entire trial at first efficacy stage
    eff_stage <- suppressWarnings(min(final_act$fir_act[grepl("Rej", final_act$dec)], na.rm = TRUE))
    if (is.finite(eff_stage)) {
      for (i in seq_len(nrow(final_act))) {
        if (final_act$fir_act[i] > eff_stage &&
            !grepl("Rej", final_act$dec[i])) {
          final_act$fir_act[i] <- eff_stage
          final_act$dec[i] <- "Stop"
        }
      }
    }
  }
  
  # Construct stop matrix: 1 = stopped arm, ratio = control arm indicator
  stop_vec <- c(final_act$fir_act, max(final_act$fir_act))
  for (i in seq_along(stop_vec)) {
    stop_matrix[stop_vec[i], i] <- if (i <= arm) 1L else ratio
  }
  
  # ==================== 6. Global null hypothesis decision ==================== #
  if (any(grepl("Rej", final_act$dec))) {
    global <- "Rej global null"
    rej_arms <- paste0("H0", final_act$arm[grepl("Rej", final_act$dec)], collapse = "+")
    rej_detail <- paste("Rej", rej_arms)
  } else {
    global <- "Acc global null"
    rej_detail <- NA
  }
  
  # ==================== 7. Trial overview text ==================== #
  trial_overview <- paste("Stop at stage", max(stop_vec), ",", global)
  trial_overview_det <- if (!is.na(rej_detail)) {
    paste("Stop at stage", max(stop_vec), ",", rej_detail)
  } else {
    NA
  }
  
  # ==================== 8. Enrolled sample size ==================== #
  enrolled_sample_size <- sum(cum_sam %*% stop_matrix)
  
  # ==================== 9. Combine results ==================== #
  output_info <- list(
    global = global,
    rej_detail = rej_detail,
    trial_overview = trial_overview,
    trial_overview_det = trial_overview_det,
    final_act = final_act,
    sta_matrix = Z_mat,
    judge_matrix = judge_mat,
    stop_matrix = stop_matrix,
    enrolled_sample_size = enrolled_sample_size
  )
  
  return(output_info)
}


#########################################################################################################
# Function: MAMS_SIM_RES
# Purpose : Perform multiple replications of MAMS_SIM and summarise simulation results
#
# Description:
#   This function repeatedly calls MAMS_SIM() for nsim times, summarises the
#   frequency of each adaptation/decision, and computes expected sample size.
#
# Outputs:
#   A list containing frequency tables, expected sample size, and all combined results.
#########################################################################################################

MAMS_SIM_RES <- function(arm,
                         stage_num,
                         cum_sam,
                         ratio,
                         delta1,
                         sd,
                         futi_bound,
                         eff_bound,
                         eff_stop = TRUE,
                         nsim = 1000,
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  sample_res <- numeric(nsim)
  arm_res <- data.frame()
  trial_stop_res <- data.frame()
  trial_stop_det_res <- data.frame()
  trial_hy_res <- data.frame()
  trial_hy_rej_del_res <- data.frame()
  
  pb <- txtProgressBar(style = 3)
  start_time <- Sys.time()
  
  for (sim in seq_len(nsim)) {
    sim_out <- MAMS_SIM(arm, stage_num, cum_sam, ratio,
                        delta1, sd, futi_bound, eff_bound, eff_stop)
    sample_res[sim] <- sim_out$enrolled_sample_size
    arm_res <- rbind(arm_res, sim_out$final_act)
    trial_stop_res <- rbind(trial_stop_res, data.frame(act = sim_out$trial_overview))
    trial_stop_det_res <- rbind(trial_stop_det_res, data.frame(act = sim_out$trial_overview_det))
    trial_hy_res <- rbind(trial_hy_res, data.frame(act = sim_out$global))
    trial_hy_rej_del_res <- rbind(trial_hy_rej_del_res, data.frame(act = sim_out$rej_detail))
    setTxtProgressBar(pb, sim / nsim)
  }
  
  close(pb)
  end_time <- Sys.time()
  
  exp_sample <- mean(sample_res)
  exp_text <- paste("Expected sample size:", round(exp_sample, 2))
  
  # Frequency of actions
  freq <- as.data.frame(table(arm_res$dec))
  names(freq) <- c("act", "Freq")
  freq$prop <- sprintf("%.1f%%", 100 * freq$Freq / nsim)
  
  cat("\nSimulation completed in", round(difftime(end_time, start_time, units = "secs"), 1), "seconds\n")
  cat(exp_text, "\n\n")
  
  final_res <- freq[order(freq$act), ]
  result_list <- list(
    nsim = nsim,
    final_res = final_res,
    expected_sample_size = exp_sample,
    sample_distribution = sample_res
  )
  
  return(result_list)
}


#########################################################################################################
# Example of use
#########################################################################################################

# library(MAMS)
# example boundaries (2-stage, 1 experimental arm)
# a <- mams(K=1, J=2, alpha=0.025, power=0.85, r=1:2, r0=1:2,
#           p=pnorm(0.65/(1.2*sqrt(2))), p0=pnorm(0.25/(1.2*sqrt(2))),
#           ushape='fixed', ufix=10, lshape='fixed', lfix=0.93)
#
# ss <- MAMS_SIM_RES(
#   arm = 1,
#   stage_num = 2,
#   cum_sam = a$n * a$rMat[1,],
#   ratio = 1,
#   delta1 = c(0.45, 0),
#   sd = 1.2,
#   futi_bound = a$l,
#   eff_bound = a$u,
#   nsim = 10000
# )
# print(ss$final_res)
