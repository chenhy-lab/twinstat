## ref: Ramakrishnan V, Goldberg J, Henderson WG, Eisen SA, True W, Lyons MJ, Tsuang MT. Elementary methods for the analysis of dichotomous outcomes in unselected samples of twins. Genet Epidemiol. 1992;9(4):273-87. doi: 10.1002/gepi.1370090406IF: 2.1 Q3 . PMID: 1398046IF: 2.1 Q3 .

#' Estimate of the familial association for disease
#'
#' \code{logisticMD} is a function for analyzing  dichotomous outcomes in unselected 
#' samples of twins using logistic regression based on unconditional
#' maximum likelihood estimation. 
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
#' Return a dataframe 
#' The result of logistic regression method.
#'
#' @examples
#' data <- load_example_data();
#' logisticMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
#' 
#' @export
logisticMD <- function(data, famid, 
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
  
  tmp_data <- data %>% dplyr::select(dplyr::all_of(c(famid,outcome,zyg))) %>% 
    dplyr::rename(FAMID=famid,Outcome=outcome,Zygo=zyg) %>% 
    dplyr::mutate(Outcome=factor(Outcome,levels = c(outcome_present,outcome_absent),labels = c('Present','Absent')),
                  Zygo=factor(Zygo,levels = c(Monochorionic_value,Dichorionic_value),labels = c('MZ','DZ'))) %>% 
    stats::na.omit()
  count_table <- tmp_data %>% dplyr::group_by(FAMID) %>% dplyr::summarise(n=dplyr::n())
  
  if(all(count_table$n == 2)){
    tmp_data$id <- rownames(tmp_data)
    ## Cross-tabulation table
    set.seed(seed)
    index_twin <- tmp_data %>% dplyr::group_by(FAMID) %>% dplyr::arrange(Outcome) %>% 
      dplyr::slice_sample(n=1) %>% dplyr::mutate(index_outcome=0)
    co_twin <- tmp_data %>% dplyr::anti_join(index_twin) %>% 
      dplyr::mutate(index_outcome=1)
    data <- merge(dplyr::select(index_twin,FAMID,indextwin=Outcome,Zygo),dplyr::select(co_twin,FAMID,cotwin=Outcome),by='FAMID') %>% 
      dplyr::mutate(indextwin=ifelse(indextwin=='Present',1,0),
                    cotwin=ifelse(cotwin=='Present',1,0),
                    Zygo=ifelse(Zygo=='MZ',1,0.5))
    # fit model
    fit1 <- stats::glm(indextwin ~ cotwin + Zygo, data = data)
    unisum <- summary(fit1)
    fit1_Beta1 <- coef(fit1)[2]
    fit1_Beta1_p <- unisum$coefficients[2, "Pr(>|t|)"]
    
    fit2 <- stats::glm(indextwin ~ cotwin + Zygo + cotwin*Zygo, data = data)
    unisum <- summary(fit2)
    fit2_Beta1 <- coef(fit2)[2]
    Beta3 <- coef(fit2)[4]
    # Beta_MZ <- fit2_Beta1+Beta3
    # Beta_DZ <- fit2_Beta1+0.5*Beta3
    fit2_Beta1_p <- unisum$coefficients[2, "Pr(>|t|)"]
    Beta3_p <- unisum$coefficients[4, "Pr(>|t|)"]
    # The independent variables used in the regression models are: X1 is the indicator for the co-twin effect, X2 is the indicator for zygosity (1 if MZ and 0.5 if DZ), and X3 is the indicator for Vietnam. 
    res_df <- dplyr::tribble(
      ~Components, ~RegressionModel, ~H0, ~Estimate, ~P_value, 
      "Familial", "Y = beta0 + beta1*X1 + beta2*X2", "beta1 = 0", sprintf("beta1 = %s",fit1_Beta1), fit1_Beta1_p,
      "Additive genetic", "Y = beta0 + beta1*X1 + beta2*X2 + beta3*X1*X2", "beta3 = 0", sprintf("beta3 = %s",Beta3), fit2_Beta1_p,
      "Common environmental", "Y = beta0 + beta1*X1 + beta2*X2 + beta3*X1*X2", "beta1 = 0", sprintf("beta1 = %s",fit2_Beta1), Beta3_p
    )
  }
  if(print.out==TRUE){
    print(res_df)
    cat("\n*X1 is the indicator for the co-twin effect, X2 is the indicator for zygosity.")
  }
  # Return
  return(res_df)
}


