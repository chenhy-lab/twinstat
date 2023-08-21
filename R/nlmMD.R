## ref: Rabe-Hesketh S, Skrondal A, Gjessing HK. Biometrical modeling of twin and family data using standard mixed model software. Biometrics. 2008 Mar;64(1):280-8. doi: 10.1111/j.1541-0420.2007.00803.xIF: 1.9 Q2 . Epub 2007 May 2. PMID: 17484777.

#' Estimate of the familial association for disease
#'
#' \code{nlmMD} is a function for analyzing  dichotomous outcomes in unselected 
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
nlmMD <- function(data, famid, 
                  outcome, 
                  outcome_present = 1,
                  outcome_absent = 0,
                  zyg,
                  Monochorionic_value = 1,
                  Dichorionic_value = 0, 
                  cov_var = NULL,
                  print.out=TRUE){
  ## Check parameters
  all_input_para <- c('data','famid','outcome','zyg')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('print.out',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('print.out can be only TRUE or FALSE!');return(FALSE)}
  if(!all(cov_var %in% colnames(data))){message('Some of the cov_var is not in the data!');return(FALSE)}
  
  tmp_data <- data %>% dplyr::select(dplyr::all_of(c(famid,outcome,zyg,cov_var))) %>% 
    dplyr::rename(FAMID=famid,Outcome=outcome,Zygo=zyg) %>% 
    dplyr::mutate(Outcome=factor(Outcome,levels = c(outcome_present,outcome_absent),labels = c('Present','Absent')),
                  Zygo=factor(Zygo,levels = c(Monochorionic_value,Dichorionic_value),labels = c('MZ','DZ')),
                  id=rownames(.)) %>% 
    stats::na.omit()
  
  dt_list <- split(tmp_data,tmp_data[,famid])
  dt_list <- lapply(dt_list,function(df){
    df <- df %>% dplyr::mutate(M=ifelse(.$Zygo=='MZ', rep(.$id[1],2), .$id))
    return(df)
  })
  dt_model <- dplyr::bind_rows(dt_list)
  dt_model$M <- as.factor(dt_model$M)
  if(is.null(cov_var)){
    use_formula <- as.formula('Outcome ~ 1 + (1|FAMID) + (1|M)')
  }else{
    use_formula <- as.formula(sprintf('Outcome ~ %s + (1|FAMID) + (1|M)', paste0(cov_var,collapse = '+')))
  }
  model_glmer <- lme4::glmer(use_formula, 
                 data = dt_model, family = binomial,
                 control=glmerControl(optimizer=c("Nelder_Mead","bobyqa"), 
                                      optCtrl=list(maxfun=100000)), 
                 na.action = "na.exclude")
  
  
}




# 输出模型结果
print(summary(model_glmer));
res <- summary(model_glmer)
sig2E.es <- res$sigma**2
sig2A.es <- 2*res$varcor[1]$M[1]
sig2C.es <- res$varcor[2]$FAMID[1]-res$varcor[1]$M[1]
h2.es <- sig2A.es/(sig2A.es+sig2C.es+sig2E.es)

sig2E_ADE.es <- res$sigma**2;
sig2A_ADE.es <- 3*res$varcor[2]$FAMID[1]-res$varcor[1]$M[1];
sig2D_ADE.es <- -2*res$varcor[2]$FAMID[1]+2*res$varcor[1]$M[1]
h2_ADE.es <- (sig2A_ADE.es+sig2D_ADE.es)/(sig2A_ADE.es+sig2D_ADE.es+sig2E_ADE.es)



sig2A.es <- 2*res$varcor["M","(Intercept)"] 
sig2C.es <- res$varcor["Pair","(Intercept)"] - res$varcor["M","(Intercept)"]
(h2.es <- sig2A.es/(sig2A.es + sig2C.es))


