#' @import dplyr
#' @importFrom tibble column_to_rownames
#' @importFrom lme4 glmer

library(tibble)
library(dplyr)
library(lme4)

## Auxiliary function
check_para <- function(para_name,envir){
  if(base::exists(para_name,envir=envir)==FALSE){message(sprintf('%s missing !',para_name));return(0)}
  if(is.null(base::get(para_name,envir=envir))==TRUE){message(sprintf('%s is NULL !',para_name));return(0)}
  return(1)
}
check_option <- function(para_name,option_list,envir){
  if(!base::get(para_name,envir=envir) %in% option_list){
    message(sprintf('Only accept %s set at: %s !',para_name,base::paste(option_list,collapse=';')));return(0)
  }
  return(1)
}

#' Load example data
#'
#' This function loads example data for testing and learning purposes.
#'
#' @return The loaded example data
#' @export
#'
#' @examples
#' data <- load_example_data()
load_example_data <- function() {
  data_path <- system.file("extdata","exampleData.RData",package = "twinstat",mustWork = TRUE)
  # Load example data here
  load(data_path)
  # Return the loaded data
  return(data)
}



#' Sample size for counting monochorionic and dichorionic twins
#'
#' \code{twin.chisqTest} is a function to count the difference in the number of monochorionic and dichorionic twins
#'
#' @param data a dataframe, raw data
#' @param famid character, family id.
#' @param outcome character, outcome variable.
#' @param outcome_present character or numeric, present in outcome variable. Default is 1.
#' @param outcome_absent character or numeric, absent in outcome variable. Default is 0.
#' @param zyg character, indicates the type of twin, in the example data 1 is monochorionic and 0 is dichorionic.
#' @param Monochorionic_value character or numeric, the value indicating monochorionic twins. Default is 1.
#' @param Dichorionic_value character or numeric, the value indicating dichorionic twins. Default is 0.
#' 
#' 
#' @return
#' Return a list. Chisq_value and P_value are the results of the chi-square test and 
#' stat_mat is the cross-tabulation table.
#'
#'
#' @examples
#' data <- load_example_data();
#' twin.chisqTest(data = data,famid = 'FAMID',outcome = 'BPD',zyg = 'zygo')
#' 
#' @export
twin.chisqTest <- function(data, famid, 
                           outcome, 
                           outcome_present = 1,
                           outcome_absent = 0,
                           zyg,
                           Monochorionic_value = 1,
                           Dichorionic_value = 0,
                           print.out=TRUE,
                           ...){
  ## Check parameters
  all_input_para <- c('data','famid','outcome','zyg')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  ## Cross-tabulation table
  tmp_data <- data %>% dplyr::select(dplyr::all_of(c(famid,outcome,zyg))) %>% 
    dplyr::rename(FAMID=famid,Outcome=outcome,Zygo=zyg) %>% 
    dplyr::mutate(Outcome=factor(Outcome,levels = c(outcome_present,outcome_absent),labels = c('Present','Absent')),
                  Zygo=factor(Zygo,levels = c(Monochorionic_value,Dichorionic_value),labels = c('MZ','DZ'))) %>% 
    stats::na.omit()
  count_table <- tmp_data %>% dplyr::group_by(FAMID) %>% dplyr::summarise(n=dplyr::n())
  if(all(count_table$n == 2)){
    tmp_data <- tmp_data %>% 
      dplyr::group_by(FAMID) %>% 
      dplyr::mutate(time=rank(FAMID,ties.method = 'first')) %>% 
      dplyr::ungroup() %>% as.data.frame() %>% 
      stats::reshape(direction = 'wide',
                     v.names = 'Outcome',      
                     idvar = 'FAMID',
                     timevar = 'time',
                     sep = '') %>% 
      dplyr::mutate(Outcome_comb=dplyr::case_when(Outcome1=='Present' & Outcome2=='Present' ~ 'Both',
                                           Outcome1=='Absent' & Outcome2=='Absent' ~ 'Neither',
                                           (Outcome1=='Present' & Outcome2=='Absent') | (Outcome1=='Absent' & Outcome2=='Present') ~ 'One',
                                           .default = 'Outlier')) %>% 
      dplyr::filter(Outcome_comb != 'Outlier')
    stat_mat <- as.data.frame.matrix(table(tmp_data$Zygo,tmp_data$Outcome_comb)) %>% dplyr::select('Neither','One','Both')
    if(any(stat_mat==0)){
      cat('The count in the cell has zero!')
      return(NULL)
    }else{
      ## chi-square test
      chisq_stat <- stats::chisq.test(stat_mat,...)
      if(print.out==TRUE){
        print(chisq_stat)
      }
      stat_mat <- stat_mat %>% 
        dplyr::mutate(Total=base::rowSums(dplyr::across(dplyr::where(is.numeric)))) %>% base::t() %>% as.data.frame() %>% 
        dplyr::mutate(Total=base::rowSums(dplyr::across(dplyr::where(is.numeric)))) %>% base::t() %>% as.data.frame()
      
      ## Output results
      result <- list()
      result$Chisq_value <- chisq_stat$statistic
      result$P_value <- chisq_stat$p.value
      result$stat_mat <- stat_mat
      return(invisible(result))
    }
  }else{
    cat('There must be two children per family!')
    return(NULL)
  }
    
}



#' Organization of Twin Data for the Genetic Analysis of Dichotomous Outcomes
#'
#' \code{twin.index} is a function to count the the number of index twins and Co-twins
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
#'
#' @return
#' Return a list. stat_MZ_mat and stat_DZ_mat are cross-tabulation tables.
#'
#'
#' @examples
#' data <- load_example_data();
#' result.twin.index <- twin.index(data = data,famid = 'FAMID',outcome = 'BPD',zyg = 'zygo')
#' result.twin.index
#' 
#' @export
twin.index <- function(data, famid, 
                       outcome, 
                       outcome_present = 1,
                       outcome_absent = 0,
                       zyg,
                       Monochorionic_value = 1,
                       Dichorionic_value = 0, 
                       seed = 22){
  ## Check parameters
  all_input_para <- c('data','famid','outcome','zyg')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
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
    data <- merge(dplyr::select(index_twin,FAMID,indextwin=Outcome,Zygo),dplyr::select(co_twin,FAMID,cotwin=Outcome),by='FAMID')
    
    t1 <- table(Co_twin_outcome=data$cotwin,Index_twin_outcome=data$indextwin,Zygo=data$Zygo) %>% as.data.frame()
    if(nrow(t1) == 8){
      t1_MZ <- t1 %>% dplyr::filter(Zygo=='MZ') %>% 
        dplyr::select(-Zygo) %>% 
        stats::reshape(idvar = "Co_twin_outcome", timevar = c("Index_twin_outcome"), direction = "wide") %>% 
        tibble::remove_rownames() %>% 
        tibble::column_to_rownames('Co_twin_outcome') %>% 
        dplyr::mutate(Total=rowSums(dplyr::across(dplyr::where(is.numeric)))) %>% base::t() %>% as.data.frame() %>% 
        dplyr::mutate(Total=rowSums(dplyr::across(dplyr::where(is.numeric)))) %>% base::t() %>% as.data.frame() %>% 
        dplyr::rename(Index.Absent=Freq.Absent,Index.Present=Freq.Present)
      t1_MZ <- t1_MZ[c('Present','Absent','Total'),c('Index.Present','Index.Absent','Total')]
      
      a1 <- t1_MZ['Present','Index.Present']
      c1 <- t1_MZ['Absent','Index.Present']
      b1 <- t1_MZ['Present','Index.Absent']
      d1 <- t1_MZ['Absent','Index.Absent']
      n11 <- sum(t1_MZ$Index.Present)
      n01 <- sum(t1_MZ$Index.Absent)
      
      t1_DZ <- t1 %>% dplyr::filter(Zygo=='DZ') %>% 
        dplyr::select(-Zygo) %>% 
        stats::reshape(idvar = "Co_twin_outcome", timevar = c("Index_twin_outcome"), direction = "wide") %>% 
        tibble::remove_rownames() %>% 
        tibble::column_to_rownames('Co_twin_outcome') %>% 
        dplyr::mutate(Total=rowSums(dplyr::across(dplyr::where(is.numeric)))) %>% base::t() %>% as.data.frame() %>% 
        dplyr::mutate(Total=rowSums(dplyr::across(dplyr::where(is.numeric)))) %>% base::t() %>% as.data.frame() %>% 
        dplyr::rename(Index.Absent=Freq.Absent,Index.Present=Freq.Present)
      t1_DZ <- t1_DZ[c('Present','Absent','Total'),c('Index.Present','Index.Absent','Total')]
      
      a2 <- t1_DZ['Present','Index.Present']
      c2 <- t1_DZ['Absent','Index.Present']
      b2 <- t1_DZ['Present','Index.Absent']
      d2 <- t1_DZ['Absent','Index.Absent']
      n12 <- sum(t1_DZ$Index.Present)
      n02 <- sum(t1_DZ$Index.Absent)
      
      ## Output results
      result <- list()
      result$stat_MZ_mat <- t1_MZ
      result$stat_DZ_mat <- t1_DZ
      result$stat_DZ <- result$stat_MZ <- list()
      result$stat_MZ$a1 <- a1
      result$stat_MZ$c1 <- c1
      result$stat_MZ$b1 <- b1
      result$stat_MZ$d1 <- d1
      result$stat_MZ$n11 <- n11
      result$stat_MZ$n01 <- n01
      result$stat_DZ$a2 <- a2
      result$stat_DZ$c2 <- c2
      result$stat_DZ$b2 <- b2
      result$stat_DZ$d2 <- d2
      result$stat_DZ$n12 <- n12
      result$stat_DZ$n02 <- n02
      return(invisible(result))
    }else{
      cat('The count in the cell has zero!')
      return(NULL)
    }
  }else{
    cat('There must be two children per family!')
    return(NULL)
  }
}



