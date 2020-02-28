#' Synthetic Control Method
#'
#' Applies the synthetic control method to a variable within a data.frame. The user must specify
#' which variable receives the treatment, a dummy variable taking a value of '0' before treatment
#' and '1' after the treatment, a variable specifying groups, and the co-variates used to
#' determine weights.
#'
#' @param data data.frame
#' @param y character, symbol, or number denoting the column containing the dependent variable
#' @param treat character, symbol, or number denoting the column containing the treatment dummy
#' @param time character, symbol, or number denoting the column containing the time variable
#' @param group character, symbol, or number denoting the column containing the grouping variable
#' @param ... characters, symbols, or numbers denoting the columns containing co-variates
#' @param op operation used on X0 for determining weights. Default is 'mean'
#' @param method optimization method accepted by \code{optim}. Default is "Nelder-Mead"
#'
#' @return list of tables containing \code{y}, the estimated counterfactual, weights, plots, the regression results
#' @importFrom kernlab ipop primal
#' @importFrom dplyr group_by_at select_at summarise_all ungroup pull
#' @importFrom tidyr gather
#' @import magrittr
#' @import ggplot2
#' @export
#'
#' @examples
#' data(ecuador)
#' results <- synth(
#'   pwt_mod,                                         # data
#'   GDPpc, Treat, year, country,                     # y, treat, time, group
#'   pop, PR, KL, xr, delta, ex, im, s, Pex, Pim, Pk  # co-variates
#' )
#'
#' # Plot Results
#' results$plot_path
#' results$plot_gap
#' results$plot_jitter
#'
#' # Show Weights
#' results$w_group
#' results$w_pred
#'
#' # Show Differences before and after
#' results$difference
#'
#' # Regression summary of Treatment on the Gap
#' summary(results$lm)
#'
#' # Data.frame with Treated and Synthetic Series
#' results$results
#'
synth <- function(
  data, y, treat, time, group, ..., op = mean, method = "Nelder-Mead"
){
  # Check if Data.Frame
  if (!is.data.frame(data)){stop("A 'data.frame' object must be supplied.")}

  # Obtain Variables
  y     <- as.character(substitute(y))
  treat <- as.character(substitute(treat))
  time  <- as.character(substitute(time))
  group <- as.character(substitute(group))
  predictors <- sapply(eval(substitute(alist(...))), function(i){as.character(i)})

  # Group receiving Treatment
  treat_group   <- unique(data[[group]][data[[treat]] == 1])
  if (length(treat_group) != 1){stop("Only one group can receive the treatment.")}
  treat_data    <- data[data[[group]] == treat_group, ]
  control_data  <- data[data[[group]] != treat_group, ]

  # Pre- & Post- Times
  pre_treat  <- treat_data[treat_data[[treat]] == 0, ]
  post_treat <- treat_data[treat_data[[treat]] == 1, ]
  pre_time   <- pre_treat[[time]]
  post_time  <- post_treat[[time]]
  treat_date <- min(post_time)

  # X1, Z1 (Treatment Group, Prior; Independent, Dependent)
  X1 <- as.matrix(apply(pre_treat[, predictors], 2, op, na.rm = TRUE))
  Z1 <- as.matrix(pre_treat[, y])

  # X0, Z0 (Control Groups, Prior; Independent, Dependent)
  pre_control <- control_data[control_data[[time]] %in% pre_time, ]
  pre_control %>%
    group_by_at(group) %>%
    select_at(predictors) %>%
    summarise_all(op, na.rm=TRUE) %>%
    ungroup() -> X0
  control_group <- pull(X0, 1)
  J  <- length(control_group)
  X0 <- t(as.matrix(X0[,-1]))
  Z0 <- matrix(pull(pre_control[, y],1), byrow = FALSE,
               nrow=length(pre_time), ncol = J)

  # Scale X's
  X <- cbind(X0, X1)
  scaled.X <- t(t(X) %*% (1/(sqrt(apply(X, 1, var))) * diag(dim(X)[1])))
  X0.scaled <- scaled.X[,-dim(scaled.X)[2]]
  X1.scaled <- scaled.X[, dim(scaled.X)[2]]
  if (is.vector(X0.scaled)){X0.scaled <- t(X0.scaled)}

  # FUNCTION TO OPTIMIZE
  cal_weights <- function(par, loss = TRUE){
    abs_par <- abs(par)
    V  <- diag(abs_par/sum(abs_par))
    H  <- t(X0.scaled) %*% V %*% (X0.scaled)
    c  <- -1*c(t(X1.scaled) %*% V %*% X0.scaled)
    nc <- length(c)
    A  <- t(rep(1, nc))
    l  <- rep(0, nc)
    u  <- rep(1, nc)
    results <- ipop(c=c, H=H, A=A, b=1, l=l, u=u, r=0)
    w    <- as.matrix(primal(results))
    if (loss){
      dW   <- Z1 - as.matrix(Z0) %*% as.matrix(w)
      return(as.numeric(t(dW) %*% dW))
    } else {return(w)}
  }

  # Optimization
  # ------------
  # First Version
  n_pred <- length(predictors)
  par1   <- rep(1/n_pred, n_pred)
  opt1   <- optim(par1, cal_weights, method = method)

  # Second Version
  X <- cbind(1, t(cbind(X1.scaled, X0.scaled)))
  Z <- cbind(Z1, Z0)
  beta  <- solve(t(X) %*% X) %*% t(X) %*% t(Z)
  beta0 <- beta[-1, ]
  V     <- diag(beta0 %*% t(beta0))
  par2  <- V/sum(V)
  opt2  <- optim(par2, cal_weights, method = method)

  # Compare Optimizations & Get Weights
  if (opt1$value < opt2$value){
    opt <- opt1$par
  } else {
    opt <- opt2$par
  }
  w_pred <- abs(opt)/sum(abs(opt))
  w_raw  <- cal_weights(w_pred, loss = FALSE)
  w_grp  <- as.numeric(abs(w_raw)/sum(abs(w_raw)))

  # Weight Table
  pred_data <- data.frame(Variable = predictors, Weights = w_pred)
  pred_data <- cbind(predictors, as.data.frame(round(cbind(
    w_pred, X1, X0 %*% as.matrix(w_grp),apply(X0, 1, mean)
  ), 4)))
  rownames(pred_data) <- NULL
  colnames(pred_data) <- c("Variable", "Weights", "Treated", "Synthetic", "Sample Mean")

  # Predictions
  results <- data.frame(
    Time      = c(pre_time, post_time),
    Treated   = pull(treat_data[, y], 1),
    Synthetic = as.numeric(matrix(pull(control_data[, y],1), byrow=FALSE, ncol = J) %*% w_grp)
  )
  results["Gap"] <- results[["Treated"]] - results[["Synthetic"]]
  results["D"]   <- c(rep(0, length(pre_time)), rep(1, length(post_time)))
  colnames(results)[1] <- time

  reg <- lm(Gap ~ D, data = results)


  # Difference Table
  pre_diff  <- mean(results[["Gap"]][results[[time]] %in% pre_time])
  post_diff <- mean(results[["Gap"]][results[[time]] %in% post_time])
  did       <- post_diff - pre_diff
  diff_tab  <- data.frame(
    Time = c("Before", "After", "Difference"),
    Difference = c(pre_diff, post_diff, did)
  )

  # plotting
  results %>%
    gather(key = "Series", value = "Value", Treated, Synthetic) ->
    results2

  results2 %>%
    ggplot(aes_string(x = time, y= "Value", color = "Series", linetype="Series")) +
    geom_line() + xlab(time)+ylab(y)+
    scale_linetype_manual(values=c("dashed", "solid")) +
    geom_vline(xintercept = treat_date, color="black", size = 1.5) ->
    path_plot

  results %>%
    ggplot(aes_string(x = time, y = "Gap")) + geom_line()+
    geom_vline(xintercept = treat_date, color="black", size = 1.5) ->
    gap_plot

  results2 %>%
    ggplot(aes_string(x = "D", y = "Gap")) +
    geom_jitter(aes(color = Series)) +
    geom_smooth(aes(color = Series), method = "lm", fill = NA) ->
    jitter_plot

  return(structure(list(
    "results"    = results,
    "w_group"    = data.frame(Group = control_group, Weights = round(w_grp, 6)),
    "w_pred"     = pred_data,
    "treat_date" = treat_date,
    "difference" = diff_tab,
    "lm"         = reg,
    "plot_path"  = path_plot,
    "plot_gap"   = gap_plot,
    "plot_jitter"= jitter_plot
  ), class = c("SynthModel", "list")))
}
#' Apply the Synthetic Control Method to Multiple Treated Groups
#'
#' Applies the synthetic control method to multiple treatment groups within a data.frame. See
#' \code{synth} for the case with only one group receiving treatment.
#'
#' @param data data.frame
#' @param y character, symbol, or number denoting the column containing the dependent variable
#' @param treat character, symbol, or number denoting the column containing the treatment dummy
#' @param time character, symbol, or number denoting the column containing the time variable
#' @param group character, symbol, or number denoting the column containing the grouping variable
#' @param ... characters, symbols, or numbers denoting the columns containing co-variates
#' @param op operation used on X0 for determining weights. Default is 'mean'
#' @param method optimization method accepted by \code{optim}. Default is "Nelder-Mead"
#'
#' @return list of tables containing \code{y}, the estimated counterfactual, weights, plots, the regression results
#' @importFrom tidyr gather
#' @import ggplot2
#' @import magrittr
#' @export
#'
#' @examples
#'
synth_many <- function(
  data, y, treat, time, group, ..., op = mean, method = "Nelder-Mead"
){
  treat <- substitute(treat)
  if (is.numeric(treat)){treat <- colnames(data)[treat]
  } else {treat <- as.character(treat)}
  group <- substitute(group)
  if (is.numeric(group)){group <- colnames(data)[group]
  } else {group <- as.character(group)}


  # Groups receiving Treatment
  treat_groups   <- unique(data[[group]][data[[treat]] == 1])
  control_groups <- unique(data[[group]][data[[treat]] == 0])
  control_data   <- data[data[[group]] != treat_group, ]

  # Loop across Treatment Groups, Run Synthetic
  Loop <- lapply(treat_groups, function(i){
    df <- rbind(data[data[[group]] == i, ], control_data)
    return(synth(df, y, treat, time, group, ..., op=op, method=method))
  })
  models <- list()
  for (i in 1:length(Loop)){
    models[[treat_groups[i]]] <- Loop[[i]]
  }

  # Compile Synthetic Results from each group into data.frame, Regression
  dat <- Reduce(rbind, lapply(models, function(i){i$results}), accumulate = FALSE)
  reg <- lm(Gap ~ D, data = dat)

  # Jitter Plot
  dat %>%
    gather(key = "Series", value = "Value", Treated, Synthetic) ->
    dat2

  dat2 %>%
    ggplot(aes_string(x = "D", y = "Gap")) +
    geom_jitter(aes(color = Series)) +
    geom_smooth(aes(color = Series), method = "lm", fill = NA) ->
    jitter_plot


  return(structure(list(
    "models"      = models,
    "lm"          = reg,
    "plot_jitter" = jitter_plot
  ), class = c("SynthMany", "list")))
}
