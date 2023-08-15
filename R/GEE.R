## ref: Lipsitz, Stuart R., et al. “Generalized Estimating Equations for Correlated Binary Data: Using the Odds Ratio as a Measure of Association.” Biometrika, vol. 78, no. 1, 1991, pp. 153–60. JSTOR, https://doi.org/10.2307/2336905.

#' Estimate of the familial association for disease
#'
#' \code{GEE} is a function for analyzing dichotomous outcomes in unselected 
#' samples of twin pairs using generalized estimating equations. 
#'
#' @param data a dataframe, raw data
#' @param famid character, family Number.
#' @param outcome character, outcome variable.
#' @param outcome_present character or numeric, present in outcome variable. Default is 1.
#' @param outcome_absent character or numeric, absent in outcome variable. Default is 0.
#' @param zyg character, indicates the type of twin, in the example data 1 is monochorionic and 0 is dichorionic.
#' @param Monochorionic_value character or numeric, the value indicating monochorionic twins. Default is 1.
#' @param Dichorionic_value character or numeric, the value indicating dichorionic twins. Default is 0.
#' @param seed numeric, the methodology will use this seed to randomly designating one member of each twin
#' pair as an “index” twin and the other member as the “co-twin.”
#' @param print.out a boolean, Whether or not to print results.
#' 
#' @return
#' Return a list. 
#' result.twin.index is the result of function \code{twin.index}.
#' chisq_G and  p_value_G are the Woolf-Haldane test statistics for the difference in the zygosity-specific odds ratios.
#'
#' @examples
#' data <- load_example_data();
#' GEE(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
#' 
#' @export
GEE <- function(data, famid, 
                outcome, 
                outcome_present = 1,
                outcome_absent = 0,
                zyg,
                Monochorionic_value = 1,
                Dichorionic_value = 0, 
                seed = 42, print.out=TRUE){
  ## Check parameters
  all_input_para <- c('data','famid','outcome','zyg')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('print.out',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  
  chisq_res <- twin.chisqTest(data=data,famid=famid,outcome=outcome,zyg=zyg,print.out=FALSE)
  stat_mat <- chisq_res$stat_mat
  n11 <- stat_mat['MZ','Both'];n12 <- stat_mat['DZ','Both'];n1 <- stat_mat['Total','Both']
  n21 <- stat_mat['MZ','One'];n22 <- stat_mat['DZ','One'];n2 <- stat_mat['Total','One']
  n31 <- stat_mat['MZ','Neither'];n32 <- stat_mat['DZ','Neither'];n3 <- stat_mat['Total','Neither']
  
  result.twin.index <- twin.index(data=data,famid=famid,outcome=outcome,zyg=zyg,seed = seed)
  twin.index.MZ <- result.twin.index$stat_MZ
  twin.index.DZ <- result.twin.index$stat_DZ
  a1 <- twin.index.MZ$a1;b1 <- twin.index.MZ$b1;c1 <- twin.index.MZ$c1
  d1 <- twin.index.MZ$d1;a2 <- twin.index.DZ$a2;b2 <- twin.index.DZ$b2
  c2 <- twin.index.DZ$c2;d2 <- twin.index.DZ$d2
  OR_MZ <- n11*n31/(n21/2)^2
  OR_DZ <- n12*n32/(n22/2)^2
  Var_OR_MZ <- exp(1/n11+4/n21+1/n31)
  Var_OR_DZ <- exp(1/n12+4/n22+1/n32)
  WMZ <- (1/a1+1/b1+1/c1+1/d1)^(-1);WDZ <- (1/a2+1/b2+1/c2+1/d2)^(-1)
  # Woolf-Haldane test
  # H0: OR(MZ) = OR(DZ) 
  chisq_G <- (log(OR_MZ)-log(OR_DZ))^2/(1/WMZ+1/WDZ)
  p_value_G <- 1-stats::pchisq(chisq_G,df=1)
  
  if(print.out==TRUE){
    cat('\n   *Woolf-Haldane test     \n\n')
    cat(sprintf('OR(MZ) = %s, OR(DZ) = %s\n',round(OR_MZ,3),round(OR_DZ,3)))
    cat(sprintf('Test of H\u2080: OR(MZ) = OR(DZ)\n  X-squared = %s, df = 1, p-value = %s\n',round(chisq_G,3),format(p_value_G, scientific = TRUE, digits = 3)))
    cat('p-value < 0.05 and OR(MZ) > OR(DZ) indicates that a specific underlying genetic influence on disease development.
p-value > 0.05 indicates that the absence of a genetic influence on disease development.\n')
  }
  # Return
  OR_df <- dplyr::tribble(
    ~Method, ~OddsRatio, ~Variance, 
    "Monozygotic", OR_MZ, Var_OR_MZ,
    "Dizygotic", OR_DZ, Var_OR_DZ
  ) %>% 
    dplyr::mutate(OR_lower=OddsRatio-1.96*sqrt(Variance),
                  OR_upper=OddsRatio+1.96*sqrt(Variance))
  
  Test_df <- dplyr::tribble(
    ~Method, ~NullHypothesis, ~Chisq, ~P_value,
    "Generalized Estimating Equations", "OR(MZ) = OR(DZ)", chisq_G, p_value_G
  )
  result <- list()
  result$result.twin.index <- result.twin.index
  result$OR <- OR_df
  result$test <- Test_df
  return(invisible(result))
}


