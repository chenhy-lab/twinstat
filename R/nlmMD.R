## ref: Ramakrishnan V, Goldberg J, Henderson WG, Eisen SA, True W, Lyons MJ, Tsuang MT. Elementary methods for the analysis of dichotomous outcomes in unselected samples of twins. Genet Epidemiol. 1992;9(4):273-87. doi: 10.1002/gepi.1370090406IF: 2.1 Q3 . PMID: 1398046IF: 2.1 Q3 .

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
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  
  tmp_data <- data %>% dplyr::select(dplyr::all_of(c(famid,outcome,zyg))) %>% 
    dplyr::rename(FAMID=famid,Outcome=outcome,Zygo=zyg) %>% 
    dplyr::mutate(Outcome=factor(Outcome,levels = c(outcome_present,outcome_absent),labels = c('Present','Absent')),
                  Zygo=factor(Zygo,levels = c(Monochorionic_value,Dichorionic_value),labels = c('MZ','DZ'))) %>% 
    stats::na.omit()
  
  dt_list <- split(data,data[,famid])
  dt_list <- lapply(dt_list,function(df){
    df <- df %>% dplyr::mutate(M=ifelse(.[[]]==1, rep(.$id[1],2), .$id))
    return(df)
  })
  two <- dplyr::bind_rows(two_list)
  two$M <- as.factor(two$M)
  model <- glmer(BPD ~ RDS + (1|FAMID) + (1|M), 
                 data = two, family = binomial,
                 control=glmerControl(optimizer=c("Nelder_Mead","bobyqa"), 
                                      optCtrl=list(maxfun=100000)), 
                 na.action = "na.exclude");
  
}




# 输出模型结果
print(summary(model));
res <- summary(model)
sig2E.es <- res$sigma**2
sig2A.es <- 2*res$varcor[1]$M[1]
sig2C.es <- res$varcor[2]$FAMID[1]-res$varcor[1]$M[1]
h2.es <- sig2A.es/(sig2A.es+sig2C.es+sig2E.es)

sig2E_ADE.es <- res$sigma**2;
sig2A_ADE.es <- 3*res$varcor[2]$FAMID[1]-res$varcor[1]$M[1];
sig2D_ADE.es <- -2*res$varcor[2]$FAMID[1]+2*res$varcor[1]$M[1]
h2_ADE.es <- (sig2A_ADE.es+sig2D_ADE.es)/(sig2A_ADE.es+sig2D_ADE.es+sig2E_ADE.es)

