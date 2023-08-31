## ref: Rabe-Hesketh S, Skrondal A, Gjessing HK. Biometrical modeling of twin and family data using standard mixed model software. Biometrics. 2008 Mar;64(1):280-8. doi: 10.1111/j.1541-0420.2007.00803.xIF: 1.9 Q2 . Epub 2007 May 2. PMID: 17484777.

#' Analysis of bivariate binomial data: Twin analysis
#'
#' \code{nlmMD} is a function for analyzing bivariate binomial data with correct for some covariates
#' based on R package mets (https://kkholst.github.io/mets/index.html).
#'
#' @param data a dataframe, raw data
#' @param famid character, family Number.
#' @param outcome character, outcome variable.
#' @param outcome_present character or numeric, present in outcome variable. Default is 1.
#' @param outcome_absent character or numeric, absent in outcome variable. Default is 0.
#' @param zyg character, indicates the type of twin, in the example data 1 is monochorionic and 0 is dichorionic.
#' @param Monochorionic_value character or numeric, the value indicating monochorionic twins. Default is 1.
#' @param Dichorionic_value character or numeric, the value indicating dichorionic twins. Default is 0.
#' @param cov_var vector, the covariant variables.
#' @param use_model character, fit the model: Pairwise odds ratio model (or), Bivariate Probit model (bp), additive gamma random effects model (gamma). Default is "bp".
#' @param type character, With random effects (un) and special functionality for polygenic random effects modelling such as ACE (ace), ADE (ade), AE (ae) and so forth. Default is "ace".
#' @param print.out a boolean, Whether or not to print results.
#' 
#' @return
#' Return a class of summary.mets.twostage
#'
#' @examples
#' data <- load_example_data();
#' nlmMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo',cov_var='RDS',use_model='bp')
#' nlmMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo',cov_var='RDS',use_model='bp',type='ade')
#' nlmMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo',cov_var='RDS',use_model='gamma',type='ace')
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
                  use_model = "bp",
                  type="ace",
                  print.out=TRUE){
  ## Check parameters
  all_input_para <- c('data','famid','outcome','zyg')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('use_model',c('or','bp','gamma'),envir=environment()),
                 check_option('type',c('un','ace','ade','ae'),envir=environment()),
                 check_option('print.out',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  if(!all(cov_var %in% colnames(data))){message('Some of the cov_var is not in the data!');return(FALSE)}
  
  tmp_data <- data %>% dplyr::select(dplyr::all_of(c(famid,outcome,zyg,cov_var))) %>% 
    dplyr::rename(FAMID=famid,Outcome=outcome,Zygo=zyg) %>% 
    dplyr::mutate(Outcome=factor(Outcome,levels = c(outcome_present,outcome_absent),labels = c('Present','Absent')),
                  Zygo=factor(Zygo,levels = c(Monochorionic_value,Dichorionic_value),labels = c('MZ','DZ')),
                  id=rownames(.)) %>% 
    dplyr::mutate(Outcome_num=as.numeric(Outcome)-1) %>% 
    stats::na.omit()
    
  if(use_model == "or"){
    use_formula <- as.formula(sprintf('Outcome ~ %s', paste0(cov_var,collapse = '+')))
    out.or <- mets::easy.binomial.twostage(use_formula,data=tmp_data,
                                  response="Outcome_num",id="FAMID",var.link=1,
                                  theta.formula=~-1+Zygo1)
    res <- summary(out.or)
  }
  
  if(use_model == "bp"){
    use_formula <- as.formula(sprintf('Outcome_num ~ %s', paste0(cov_var,collapse = '+')))
    out.bp <- mets::bptwin(use_formula,data=tmp_data,id="FAMID",zyg="Zygo",DZ="DZ",type=type)
    res <- summary(out.bp)
  }
  
  if(use_model == "gamma"){
    use_formula <- as.formula(sprintf('Outcome_num ~ %s', paste0(cov_var,collapse = '+')))
    theta.des <- model.matrix( ~-1+Zygo,data=tmp_data)
    margbin <- stats::glm(use_formula,data=tmp_data,family=binomial())
    if(type == 'un'){
      out.gamma <- mets::binomial.twostage(margbin,data=tmp_data,model="gamma",
                                           clusters=tmp_data$FAMID,detail=0,theta=c(0.1)/1,var.link=1,
                                           theta.des=theta.des)
      res <- summary(out.gamma)
    }else{
      out.tmp <- mets::twin.polygen.design(tmp_data,id="FAMID",zygname="Zygo",zyg="DZ",type=type)
      out.polygen <- mets::binomial.twostage(margbin,data=tmp_data,
                                             clusters=tmp_data$FAMID,detail=0,theta=c(0.1)/1,var.link=0,
                                             random.design=out.tmp$des.rv,theta.des=out.tmp$pardes)
      res <- summary(out.polygen)
      # mets::concordanceTwinACE(out.polygen,type="ace")
    }
  }
  if(print.out==TRUE){
    print(res)
  }
  return(res)
}





