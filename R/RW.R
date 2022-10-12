RW<- function(formula, data, alpha = 0.05, verbose = TRUE, p_shape) {

  data <- model.frame(formula, data)
  dp <- as.character(formula)




  a <- param_estimates(formula, data,p_shape)$a
  n <- param_estimates(formula,data,p_shape)$n
  m <- param_estimates(formula,data,p_shape)$m
  muhat <- param_estimates(formula,data,p_shape)$muhat
  sigmahat <- param_estimates(formula,data,p_shape)$sigmahat

  p <- c(rep(p_shape,a))
  k <-2*p-3
  nu <- 2*p-1


  bigm <- (2*p*m)/k

  w_MML <- bigm/(sigmahat^2)

  TW <-sum(muhat^2*w_MML)-((sum((muhat*w_MML)))^2/sum(w_MML));
  W_star <- (TW/(a-1))/(1+(2*(a-2)/((a^2)-1))*(sum((1/(n-1))*(1-(w_MML/sum(w_MML)))^2)))

  df_denom <- ((3/(a^2-1))*(sum((1/(n-1))*(1-(w_MML/sum(w_MML)))^2)))^(-1)


  pvalue <- 1-pf(W_star,a-1,df_denom);

  if (verbose) {
    cat("\n", "","Robust Welch Test based on MML Estimators", paste("(alpha = ",alpha,")",sep = ""), "\n",
        sep = " ")
    cat({"-------------------------------------------------------------"},
        "\n", sep = " ")
    cat("  data :", paste(dp[[2L]], "and", dp[[3L]]), "\n\n", sep = " ")
    cat("  statistic  :", W_star, "\n", sep = " ")
    cat("  num df     :", a-1, "\n", sep = " ")
    cat("  denom df   :", df_denom, "\n", sep = " ")
    cat("  p.value    :", pvalue, "\n\n", sep = " ")
    cat(if (pvalue > alpha) {
      "  Result     : Difference is not statistically significant."
    }
    else {
      "  Result     : Difference is statistically significant."
    }, "\n")
    cat({"-------------------------------------------------------------"},
        "\n", sep = " ")
  }

  result <- list()
  result$statistic <- W_star
  result$dfs <- c(a-1,df_denom)
  result$p.value <- pvalue
  result$alpha <- alpha
  result$method <- "Robust Welch Test based on MML Estimators"
  result$data <- data
  result$formula <- formula


  attr(result, "class") <- "htest"
  invisible(result)


}




























