RPB <- function(formula, data, alpha = 0.05, verbose = TRUE, p_shape, repn=5000){


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
  var_MML <- (sigmahat^2)/bigm
  w_MML <- bigm/(sigmahat^2)


  gmean <- sum((w_MML*muhat))/sum(w_MML)
  T0 <- sum(((muhat-gmean)/(sqrt(var_MML)))^2)


  pv=0;
  for(i in 1:repn){
    zi <- rnorm(a);
    ui <- rchisq(a,n-1);
    PB <- (sum(zi^2*(n-1)/ui))-(((sum((sqrt(bigm)*zi*(n-1))/(sigmahat*ui)))^2)/(sum(bigm*(n-1)/((sigmahat^2)*ui))));
    if(PB>T0){pv=pv+1}
  }
  pvalue <- pv/repn;

  if (verbose) {
    cat("\n", "","Robust Parametric Bootstrap Test based on MML Estimators", paste("(alpha = ",alpha,")",sep = ""), "\n",
        sep = " ")
    cat({"-------------------------------------------------------------"},
        "\n", sep = " ")
    cat("  data :", paste(dp[[2L]], "and", dp[[3L]]), "\n\n", sep = " ")
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
  result$p.value <- pvalue
  result$alpha <- alpha
  result$method <- "Robust Parametric Bootstrap Test based on MML Estimators"
  result$data <- data
  result$formula <- formula


  attr(result, "class") <- "htest"
  invisible(result)


}
RPB <- function(formula, data, alpha = 0.05, verbose = TRUE, p_shape, repn=5000){


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
  var_MML <- (sigmahat^2)/bigm
  w_MML <- bigm/(sigmahat^2)


  gmean <- sum((w_MML*muhat))/sum(w_MML)
  T0 <- sum(((muhat-gmean)/(sqrt(var_MML)))^2)


  pv=0;
  for(i in 1:repn){
    zi <- rnorm(a);
    ui <- rchisq(a,n-1);
    PB <- (sum(zi^2*(n-1)/ui))-(((sum((sqrt(bigm)*zi*(n-1))/(sigmahat*ui)))^2)/(sum(bigm*(n-1)/((sigmahat^2)*ui))));
    if(PB>T0){pv=pv+1}
  }
  pvalue <- pv/repn;

  if (verbose) {
    cat("\n", "","Robust Parametric Bootstrap Test based on MML Estimators", paste("(alpha = ",alpha,")",sep = ""), "\n",
        sep = " ")
    cat({"-------------------------------------------------------------"},
        "\n", sep = " ")
    cat("  data :", paste(dp[[2L]], "and", dp[[3L]]), "\n\n", sep = " ")
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
  result$p.value <- pvalue
  result$alpha <- alpha
  result$method <- "Robust Parametric Bootstrap Test based on MML Estimators"
  result$data <- data
  result$formula <- formula


  attr(result, "class") <- "htest"
  invisible(result)


}
