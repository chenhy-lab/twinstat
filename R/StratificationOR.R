## ref: Ramakrishnan V, Goldberg J, Henderson WG, Eisen SA, True W, Lyons MJ, Tsuang MT. Elementary methods for the analysis of dichotomous outcomes in unselected samples of twins. Genet Epidemiol. 1992;9(4):273-87. doi: 10.1002/gepi.1370090406IF: 2.1 Q3 . PMID: 1398046IF: 2.1 Q3 .

#' Estimate of the familial association for disease
#'
#' \code{compareMD} is a function for analyzing dichotomous outcomes in unselected 
#' samples of twin pairs using stratified estimators of the odds ratio. 
#'
#' @param result.twin.index a list, result of function twin.index.
#' @param print.out a boolean, Whether or not to print results.
#' 
#' @return
#' Return two dataframe. 
#' OR is different effects odds ratios.
#' test is different test statistics for the common odds ratio.
#'
#' @examples
#' data <- load_example_data();
#' twin.index_res <- twin.index(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
#' compareMD(twin.index_res)
#' 
#' @export
compareMD <- function(result.twin.index, print.out=TRUE){
  ## Check parameters
  all_input_para <- c('result.twin.index')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('print.out',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  
  twin.index.MZ <- result.twin.index$stat_MZ
  twin.index.DZ <- result.twin.index$stat_DZ
  a1 <- twin.index.MZ$a1;b1 <- twin.index.MZ$b1;c1 <- twin.index.MZ$c1
  d1 <- twin.index.MZ$d1;a2 <- twin.index.DZ$a2;b2 <- twin.index.DZ$b2
  c2 <- twin.index.DZ$c2;d2 <- twin.index.DZ$d2
  n11 <- a1+c1;n01 <- b1+d1;m11 <- a1+b1;m01 <- c1+d1;N1 <- n11+n01+m11+m01
  n12 <- a2+c2;n02 <- b2+d2;m12 <- a2+b2;m02 <- c2+d2;N2 <- n12+n02+m12+m02
  # zygosity-specific odds ratios
  OR_MZ <- (a1*d1)/(b1*c1);OR_DZ <- (a2*d2)/(b2*c2)
  # common environmental effects
  OR_c <- OR_DZ^2/OR_MZ;OR_a <- OR_DZ/OR_c
  WMZ <- (1/a1+1/b1+1/c1+1/d1)^(-1);WDZ <- (1/a2+1/b2+1/c2+1/d2)^(-1)
  # common odds ratio
  OR_common <- exp((WMZ*log(OR_MZ)+WDZ*log(OR_DZ))/(WMZ+WDZ))
  Var_OR_common <- 1/(WMZ+WDZ)
  Var_OR_MZ <- exp(1/WMZ);Var_OR_DZ <- exp(1/WDZ)
  Var_OR_c <- OR_DZ^2/Var_OR_MZ;Var_OR_a <- Var_OR_DZ/Var_OR_c
  # Mantel-Haenszel test
  # H0: OR_common = 1
  chisq_MH <- (((a1+a2)-(n11*m11/N1+n12*m12/N2))^2)/sum((n11*n01*m11*m01)/(N1^2*(N1-1)),(n12*n02*m12*m02)/(N2^2*(N2-1)))
  p_value_MH <- 1-stats::pchisq(chisq_MH,df=1)
  
  # Woolf-Haldane test
  # H0: OR(MZ) = OR(DZ) 
  chisq_G <- (log(OR_MZ)-log(OR_DZ))^2/(1/WMZ+1/WDZ)
  p_value_G <- 1-stats::pchisq(chisq_G,df=1)
  
  # Common environmental effect OR(c)
  # H0: OR(MZ) = OR(DZ)^2
  # A significant value of chisq_c (when OR(MZ) < OR(DZ)^2, OR_c > 1) would suggest that the disease
  # association between index twins and co-twins is due to both the common environment
  # and additive genes. Alternatively, if there is no evidence to reject the null hypothesis
  # we would conclude that the disease association in index twins and co-twins is due to
  # purely additive genetic effects.
  # A significant value of (when chisq_c (OR(MZ) > OR(DZ)^2, OR_c < 1) suggest that the 
  # disease association between index twins and co-twins is due to additive genes as 
  # well as some other type of genetic effects.
  chisq_c <- (log(OR_MZ)-2*log(OR_DZ))^2/(1/WMZ+4/WDZ)
  p_value_c <- 1-stats::pchisq(chisq_c,df=1)
  
  if(print.out==TRUE){
    cat('\n   *Mantel-Haenszel test     \n\n')
    cat(sprintf('Common OR = %s\n',round(OR_common,3)))
    cat(sprintf('Test of H\u2080: OR = 1\n  X-squared = %s, df = 1, p-value = %s\n',round(chisq_MH,3),format(p_value_MH, scientific = TRUE, digits = 3))) #There is no evidence for a familial influence on disease
    cat('p-value < 0.05 indicates that the familial aggregation of disease, there is a significantly higher risk of disease in the index twin if the co-twin has the disease.\n')
    cat('\n   *Woolf-Haldane test     \n\n')
    cat(sprintf('OR(MZ) = %s, OR(DZ) = %s\n',round(OR_MZ,3),round(OR_DZ,3)))
    cat(sprintf('Test of H\u2080: OR(MZ) = OR(DZ)\n  X-squared = %s, df = 1, p-value = %s\n',round(chisq_G,3),format(p_value_G, scientific = TRUE, digits = 3)))
    cat('p-value < 0.05 and OR(MZ) > OR(DZ) indicates that a specific underlying genetic influence on disease development.
p-value > 0.05 indicates that the absence of a genetic influence on disease development.\n')
    cat('\n   *Viswanathan-Ramakrishnan test     \n\n')
    cat(sprintf('OR(a) = %s, OR(c) = %s\n',round(OR_c,3),round(OR_a,3)))
    cat(sprintf('Test of H\u2080: OR(MZ) = OR(DZ)^2\n  X-squared = %s, df = 1, p-value = %s\n',round(chisq_c,3),format(p_value_c, scientific = TRUE, digits = 3)))
    cat('p-value < 0.05 and OR(MZ) > OR(DZ) indicates that disease association between index twins and co-twins is due to additive genes and some other type of genetic effects.
p-value < 0.05 and OR(MZ) < OR(DZ) indicates that disease association between index twins and co-twins is due to both the common environment and additive genes.
p-value > 0.05 indicates that the zygosity specific odds ratios are due solely to additive genetic effects and not to common environment.')
  }
  # Return
  OR_df <- dplyr::tribble(
    ~Components, ~OddsRatio, ~Variance, 
    "Monozygotic", OR_MZ, Var_OR_MZ,
    "Dizygotic", OR_DZ, Var_OR_DZ,
    "Common odds ratio", OR_common, Var_OR_common, 
    "Additive genetic", OR_a, Var_OR_a,
    "Common environmental", OR_c, Var_OR_c
  ) %>% 
    dplyr::mutate(OR_lower=OddsRatio-1.96*sqrt(Variance),
                  OR_upper=OddsRatio+1.96*sqrt(Variance))
  
  Test_df <- dplyr::tribble(
    ~Method, ~NullHypothesis, ~Chisq, ~P_value,
    "Mantel-Haenszel test", "OR(Common environmental effects)=1", chisq_MH, p_value_MH,
    "Woolf-Haldane test", "OR(MZ) = OR(DZ)", chisq_G, p_value_G,
    "Viswanathan-Ramakrishnan test", "OR(MZ) = OR(DZ)^2", chisq_c, p_value_c
  )
  result <- list()
  result$OR <- OR_df
  result$test <- Test_df
  return(invisible(result))
}







