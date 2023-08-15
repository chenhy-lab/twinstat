## ref: Donner A, Klar N, Eliasziw M. Statistical methodology for estimating twin similarity with respect to a dichotomous trait. Genet Epidemiol. 1995;12(3):267-77. doi: 10.1002/gepi.1370120304. PMID: 7557348.
#' Intraclass correlation
#'
#' \code{get_ICC} is a function measures twin similarity with the intraclass correlation (ICC) 
#' of the phenotype separately within monochorionic and dichorionic pairs.
#'
#' @param n_neither a numeric, number of samples in which neither of the twins showed the relevant phenotype.
#' @param n_one a numeric, number of samples in which only one of the twins showed the relevant phenotype.
#' @param n_both a numeric, number of samples in which both twins showed the relevant phenotype.
#' @param method character, methods for parameter estimation (e.g. "auto","normal","gof"). Default is "auto".
#' Recommended method is 'normal' when total number of sample N >100 and 
#' the prevalence parameter pi (0.3<pi<0.7) and the level of twin similarity
#' with respect to the dichotomous trait in question Rho is no greater than 0.60.
#' Recommended method is 'gof' when total number of sample N is small or moderate in size (<=100), 
#' pi is fairly extreme (<0.3 or >0.7), or p is anticipated to be of substantial magnitude (>=0.6)
#' @param print.out a boolean, Wwether or not to print results.
#' 
#' @return
#' Return a dataframe. 
#' pi is the estimated prevalence of trait.
#' Vh is the standard error square of ICC.
#' ICC means the the level of twin similarity with respect to the dichotomous trait.
#' Z0 and p_value is the statistic using the normal theory method.
#' CI_lower and CI_upper is the 95% confidence interval of ICC.
#'
#' @examples
#' data <- load_example_data();
#' chisq_res <- twin.chisqTest(data = data,famid = 'FAMID',outcome = 'BPD',zyg = 'zygo')
#' stat_mat <- chisq_res$stat_mat
#' get_ICC(n_neither = stat_mat['Monochorionic','Neither'],
#'         n_one     = stat_mat['Monochorionic','One'],
#'         n_both    = stat_mat['Monochorionic','Both'])
#' 
#' @export
get_ICC <- function(n_neither, n_one, n_both, method='auto', print.out=TRUE){
  ## Check parameters
  all_input_para <- c('n_neither','n_one','n_both','method')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('method',c('auto','normal','gof'),envir=environment()),
                 check_option('print.out',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  
  ## Statistical analysis
  N <- n_both+n_one+n_neither
  pi_est <- (2*n_both+n_one)/(2*N)
  p <- 1-n_one/(2*N*pi_est*(1-pi_est))
  Z0 <- p*sqrt(N)
  p_value <- 2*stats::pnorm(abs(Z0),lower.tail = FALSE)
  Vh <- ((1-p)/N)*((1-p)*(1-2*p)+(p*(2-p))/(2*pi_est*(1-pi_est)))
  if(method=='auto'){
    if(N >100 & pi_est>0.3 & pi_est<0.7 & p<0.6){
      method <- 'normal'
    }
    if(N<=100 | pi_est<0.3 | pi_est>0.7 | p>=0.6){
      method <- 'gof'
    }
  }
  if(method=='normal'){
    se_p <- sqrt(Vh)
    z_alpha <- 1.96
    p_low <- round(p-z_alpha*se_p,3)
    p_up <- round(p+z_alpha*se_p,3)
  }
  if(method=='gof'){
    # p1 <- n_both/N;p2 <- n_one/N;p3 <- n_neither/N
    # chisq_g <- function(nx,px){(nx-N*px)^2/(N*px)}
    # chisq_G <- sum(chisq_g(n_both,p1),chisq_g(n_one,p2),chisq_g(n_neither,p3))
    # p_value <- 1-stats::pchisq(chisq_G,df=1)
    chisq_G <- 3.84
    y1 <- ((n_one-2*N*pi_est*(1-pi_est))^2+4*N^2*pi_est^2*(1-pi_est)^2)/(4*N*pi_est^2*(1-pi_est)^2*(chisq_G+N))-1
    y2 <- (n_one^2-4*N*pi_est*(1-pi_est)*(1-4*pi_est*(1-pi_est))*chisq_G)/(4*N*pi_est^2*(1-pi_est)^2*(chisq_G+N))-1
    y3 <- (n_one+(1-2*pi_est*(1-pi_est))*chisq_G)/(pi_est*(1-pi_est)*(chisq_G+N))-1
    V <- 1/27*y3^3-1/6*(y2*y3-3*y1)
    W <- sqrt((1/9*y3^2-1/3*y2)^3)
    theta <- acos(V/W)
    pi_const <- 3.1415927
    p_low <- sqrt(1/9*y3^2-1/3*y2)*(cos((theta+2*pi_const)/3)+sqrt(3)*sin((theta+2*pi_const)/3))-1/3*y3
    p_up <- 2*sqrt(1/9*y3^2-1/3*y2)*cos((theta+5*pi_const)/3)-1/3*y3
  }
  
  ## Output results
  method_name <- ifelse(method=='normal','Normal theory method','Goodness-of-fit approach')
  if(print.out==TRUE){
    cat(sprintf('\n   Using method: %s     \n\n',method_name))
    cat(sprintf('Estimated prevalence of trait \u03C0: %s\n',round(pi_est,3)))
    cat(sprintf('Estimated prevalence of trait \u03C1(ICC): %s\n',round(p,3)))
    cat(sprintf('Test of H\u2080: \u03C1 = 0\n  Z\u2080 = %s, (p-value = %s)\n',round(Z0,3),format(p_value, scientific = TRUE, digits = 3)))
    cat(sprintf('95%% Confidence interval for \u03C1: Lower:%s,Upper:%s\n\n',round(p_low,3),round(p_up,3)))
  }
  result <- list()
  result$pi <- pi_est
  result$Vh <- Vh
  result$ICC <- p
  result$Z0 <- Z0
  result$p_value <- p_value
  result$CI_lower <- p_low
  result$CI_upper <- p_up
  return(invisible(as.data.frame(result)))
}


#' Compare intraclass correlation between two groups
#'
#' \code{compareICC} is a function to report CIs for the ICCs using a goodness-of fit
#' approach and tested whether the 2 ICCs differ.
#'
#' @param chisq_res a dataframe, the result of the twin.chisqTest function.
#' @param method character, methods for parameter estimation (e.g. "auto","normal","gof"). Default is "auto".
#' Recommended method is 'normal' when total number of monochorionic sample >=100 or dichorionic sample >=100.
#' #' Recommended method is 'gof' when total number of monochorionic sample <100 or dichorionic sample <100.
#' @param print.out a boolean, Wwether or not to print results.
#'
#' @return
#' Return a data frame. Rows are genes/drivers, columns are "ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Z-statistics", "Ave.G1" and "Ave.G0".
#' Names of the columns may vary from different group names. Sorted by P-values.
#'
#'
#' @examples
#' data <- load_example_data();
#' chisq_res <- twin.chisqTest(data = data,famid = 'FAMID',outcome = 'BPD',zyg = 'zygo')
#' compareICC(chisq_res$stat_mat)
#' 
#' @export
compareICC <- function(chisq_res,method = 'auto',print.out=TRUE){
  ## Check parameters
  all_input_para <- c('chisq_res')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('method',c('auto','normal','gof'),envir=environment()),
                 check_option('print.out',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  
  M_Neither=chisq_res['MZ','Neither']
  M_One=chisq_res['MZ','One']
  M_Both=chisq_res['MZ','Both']
  D_Neither=chisq_res['DZ','Neither']
  D_One=chisq_res['DZ','One']
  D_Both=chisq_res['DZ','Both']
  
  ## Statistical analysis
  Z <- chisq_G <- NULL
  N1 <- sum(M_Both,M_One,M_Neither);N2 <- sum(D_Both,D_One,D_Neither)
  p1 <- get_ICC(M_Neither, M_One, M_Both, print.out=FALSE)$ICC
  p2 <- get_ICC(D_Neither, D_One, D_Both, print.out=FALSE)$ICC
  pi1 <- get_ICC(M_Neither, M_One, M_Both, print.out=FALSE)$pi
  pi2 <- get_ICC(D_Neither, D_One, D_Both, print.out=FALSE)$pi
  
  if(method=='auto'){
    if(N1>100 | N2>100){
      method <- 'normal'
    }
    if(N1<100 & N2<100){
      method <- 'gof'
    }
  }
  
  if(method=='normal'){
    V1 <- get_ICC(M_Neither,M_One,M_Both, print.out=FALSE)$Vh
    V2 <- get_ICC(D_Neither,D_One,D_Both, print.out=FALSE)$Vh
    Z <- (p1-p2)/sqrt(V1+V2)
    p_value <- 2*stats::pnorm(abs(Z),lower.tail = FALSE)
  }
  if(method=='gof'){
    get_molecule <- function(Nh,pih,ph){Nh*pih*(1-pih)*ph}
    get_denominator <- function(Nh,pih,ph){Nh*pih*(1-pih)}
    p_mean <- sum(get_molecule(N1,pi1,p1),get_molecule(N2,pi2,p2))/sum(get_denominator(N1,pi1,p1),get_denominator(N2,pi2,p2))
    p1 <- p2 <- p_mean
    p11 <- pi1^2+pi1*(1-pi1)*p1;p21 <- 2*pi1*(1-pi1)*(1-p1);p31 <- (1-pi1)^2+pi1*(1-pi1)*p1
    p12 <- pi2^2+pi2*(1-pi2)*p2;p22 <- 2*pi2*(1-pi2)*(1-p2);p32 <- (1-pi2)^2+pi2*(1-pi2)*p2
    chisq_g <- function(nh,ph,Nh){(nh-Nh*ph)^2/(Nh*ph)}
    chisq_G <- sum(chisq_g(M_Both,p11,N1),chisq_g(M_One,p21,N1),chisq_g(M_Neither,p31,N1),
                   chisq_g(D_Both,p12,N2),chisq_g(D_One,p22,N2),chisq_g(D_Neither,p32,N2))
    p_value <- 1-stats::pchisq(chisq_G,df=1)
  }
  
  ## Output results
  method_name <- ifelse(method=='normal','Normal theory method','Goodness-of-fit approach')
  if(print.out==TRUE){
    cat(sprintf('\n   Using method: %s     \n\n',method_name))
    cat('Test of H\u2080: \u03C1\u2081 = \u03C1\u2082\n')
    if(method=='normal'){
      cat(sprintf('  Z = %s, (p-value = %s)\n',round(Z,3),format(p_value, scientific = TRUE, digits = 3)))
    }else{
      cat(sprintf('  X-squared = %s, (p-value = %s)\n',round(chisq_G,3),format(p_value, scientific = TRUE, digits = 3)))
    }
  }
  result <- list()
  if(method=='normal'){
    result$Z <- Z
  }else{
    result$`X-squared` <- chisq_G
  }
  result$p_value <- p_value
  return(invisible(result))
}


