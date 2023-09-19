#' Knit report using working environment
#'
#' Warning: It is quite likely that this will be called within an Rmd file
#' implying a recursive call to \code{knit()}. This will generate "duplicate label"
#' errors for unlabelled chunks. To avoid this, all code chunks
#' in your Rmd file should be named.
#' Supposedly this error can also be avoided by setting the following option:
#' \code{options(knitr.duplicate.label = 'allow')}.
#' I tried this but it didn't seem to help.
#'
#' @param input_filename Rmd file.
#' @param output_filename Markdown or HTML output file.  An HTML file
#' is specified using the .htm, .html, .HTM or .HTML file extension.
#' When html is specified, a similarly named markdown file is also
#' generated.
#' All output files including cache and figures will appear in the
#' same folder as \code{output_filename}.
#' 
#' @return NULL
#' @keywords internal
knit_report <- function(input_filename, output_filename)
{
    output_filename <- normalizePath(output_filename)

    output_dir <- dirname(output_filename)
    if (!file.exists(output_dir))
        dir.create(output_dir)

    current_dir <- getwd()
    on.exit(setwd(current_dir))
    setwd(output_dir)

    name <- gsub("\\.[^.]+$", "", basename(output_filename))
    suffix <- gsub(".*\\.([^.]+)$", "\\1", output_filename)

    is.html <- tolower(suffix) %in% c("htm","html")
    is.pdf <- tolower(suffix) == "pdf"
    is.docx <- tolower(suffix) %in% c("doc", "docx", "word")
    is.md <- tolower(suffix) %in% c("md", "markdown")

    if (is.html)
        return(knitr::knit2html(input_filename, output=paste0(name, ".html"), envir=parent.frame()))
    else if (is.md)
        return(knitr::knit(input_filename, output=paste0(name, ".md"), envir=parent.frame()))
    else if (is.pdf)
    {        
        return(rmarkdown::render(input_filename, rmarkdown::pdf_document(), intermediates_dir=getwd(), output_dir=getwd(), output_file=paste0(name, ".pdf"), clean = TRUE, envir=parent.frame()))
    }
    else if (is.docx)
    {        
        return(rmarkdown::render(input_filename, rmarkdown::word_document(), intermediates_dir=getwd(), output_dir=getwd(), output_file=paste0(name, ".docx"), clean = TRUE, envir=parent.frame()))
    }
    else
        stop("Please choose a filename with pdf, html, docx or md suffix")
}


#' Generate twinstat report
#'
#' \code{twinstat_report} this report will generate a report containing tables and graphs summarising the results.
#'
#' @param data a dataframe, raw data
#' @param famid character, family Number.
#' @param outcome character, outcome variable.
#' @param outcome_present character or numeric, present in outcome variable. Default is 1.
#' @param outcome_absent character or numeric, absent in outcome variable. Default is 0.
#' @param zyg character, indicates the type of twin, in the example data 1 is monochorionic and 0 is dichorionic.
#' @param Monochorionic_value character or numeric, the value indicating monochorionic twins. Default is 1.
#' @param Dichorionic_value character or numeric, the value indicating dichorionic twins. Default is 0.
#' @param seed numeric, see function [logisticMD()] and [GEE()].
#' @param cov_var vector, see function [metsMD()].
#' @param filename character, Markdown or HTML output file.  An HTML file
#' is specified using the .htm, .html, .HTM or .HTML file extension.
#' When html is specified, a similarly named markdown file is also
#' generated.
#' @examples
#' \dontrun{
#' data <- load_example_data()
#' data <- data %>% mutate_at(c("Sex","RDS","INST"),function(x)as.factor(x))
#' cov_var <- c("Sex","GA","BW","RDS","INST")
#' twinstat_report(data=data,famid='FAMID',outcome='BPD',zyg='zygo',cov_var=cov_var)
#' }
#' @export
#' @return NULL
twinstat_report <- function(data, famid, 
                            outcome, 
                            outcome_present = 1,
                            outcome_absent = 0,
                            zyg,
                            Monochorionic_value = 1,
                            Dichorionic_value = 0, 
                            seed = 42,
                            cov_var = NULL,
                            filename = "./report.html") {
  #
  all_input_para <- c('data','famid','outcome','zyg')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  message("Performing analysis")
  # ICC
  {
    chisq_res <- twin.chisqTest(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,print.out=F)
    stat_mat <- chisq_res$stat_mat
    ICC_MZ <- get_ICC(n_neither=stat_mat['MZ','Neither'],n_one=stat_mat['MZ','One'],n_both=stat_mat['MZ','Both'],print.out = F)
    ICC_DZ <- get_ICC(n_neither=stat_mat['DZ','Neither'],n_one=stat_mat['DZ','One'],n_both=stat_mat['DZ','Both'],print.out = F)
    MZvsDZ <- compareICC(stat_mat,print.out =F)
    stat_mat$Chisq <- NA;stat_mat[1,'Chisq'] <- base::round(chisq_res$Chisq_value,2)
    stat_mat$P <- NA;stat_mat[1,'P'] <- base::round(chisq_res$P_value,2)
    stat_mat$ICC <- NA;
    stat_mat['MZ','ICC'] <- base::paste0(base::round(ICC_MZ$ICC*100,1),'% (',base::round(ICC_MZ$CI_lower*100,1),'%, ',base::round(ICC_MZ$CI_upper*100,1),'%)')
    stat_mat['DZ','ICC'] <- base::paste0(base::round(ICC_DZ$ICC*100,1),'% (',base::round(ICC_DZ$CI_lower*100,1),'%, ',base::round(ICC_DZ$CI_upper*100,1),'%)')
    stat_mat$`P(MZvsDZ)` <- NA;stat_mat[1,'P(MZvsDZ)'] <- base::round(MZvsDZ$p_value,2)
    ICC_res <- stat_mat
  }
  # strataOR
  {
    twin.index_res <- twin.index(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,seed = seed)
    strataOR_res <- compareMD(twin.index_res,print.out = F)
  }
  # GEE
  {
    GEE_res <- GEE(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,seed = seed,print.out = F)
  }
  # logisticMD
  {
    glm_res <- logisticMD(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,seed = seed,print.out = F)
  }
  # mets
  {
    mets_list <- list()
    mets_list[['ace']] <- metsMD(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,cov_var=cov_var,use_model='bp',type = 'ace',print.out=F)
    mets_list[['ade']] <- metsMD(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,cov_var=cov_var,use_model='bp',type = 'ade',print.out=F)
    mets_list[['ae']] <- metsMD(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,cov_var=cov_var,use_model='bp',type = 'ae',print.out=F)
    mets_list[['ce']] <- metsMD(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,cov_var=cov_var,use_model='bp',type = 'ce',print.out=F)
    mets_list[['de']] <- metsMD(data=data,famid=famid,outcome=outcome,outcome_present=outcome_present,outcome_absent=outcome_absent,zyg=zyg,Monochorionic_value=Monochorionic_value,Dichorionic_value=Dichorionic_value,cov_var=cov_var,use_model='bp',type = 'de',print.out=F)
    coef_list <- list()
    for (model_type in c('ace','ade','ae','ce','de')) {
      model_res <- mets_list[[model_type]]
      if(model_type %in% c('ace')){
        coef_list[[model_type]] <- c(paste0(round(model_res$coef[,1],3),' (',round(model_res$coef[,3],3),' ~ ',round(model_res$coef[,4],3),')')[1:3],round(-2*as.numeric(model_res$logLik),2),round(model_res$AIC,2))
        names(coef_list[[model_type]]) <- c('A','C','E','-2LL','AIC')
      }
      if(model_type %in% c('ade')){
        coef_list[[model_type]] <- c(paste0(round(model_res$coef[,1],3),' (',round(model_res$coef[,3],3),' ~ ',round(model_res$coef[,4],3),')')[1:3],round(-2*as.numeric(model_res$logLik),2),round(model_res$AIC,2))
        names(coef_list[[model_type]]) <- c('A','D','E','-2LL','AIC')
      }
      if(model_type %in% c('ae')){
        coef_list[[model_type]] <- c(paste0(round(model_res$coef[,1],3),' (',round(model_res$coef[,3],3),' ~ ',round(model_res$coef[,4],3),')')[1:2],round(-2*as.numeric(model_res$logLik),2),round(model_res$AIC,2))
        names(coef_list[[model_type]]) <- c('A','E','-2LL','AIC')
      }
      if(model_type %in% c('ce')){
        coef_list[[model_type]] <- c(paste0(round(model_res$coef[,1],3),' (',round(model_res$coef[,3],3),' ~ ',round(model_res$coef[,4],3),')')[1:2],round(-2*as.numeric(model_res$logLik),2),round(model_res$AIC,2))
        names(coef_list[[model_type]]) <- c('C','E','-2LL','AIC')
      }
      if(model_type %in% c('de')){
        coef_list[[model_type]] <- paste0(round(model_res$coef[,1],3),' (',round(model_res$coef[,3],3),' ~ ',round(model_res$coef[,4],3),')')[1:2]
        names(coef_list[[model_type]]) <- c('D','E')
      }
    }
    mets_res <- list()
    mets_res$adj_var <- paste0(cov_var,collapse = ',')
    mets_res$ace_cov <- round(mets_list[['ace']]$par,4) %>% as.data.frame()
    mets_res$ade_cov <- round(mets_list[['ade']]$par,4) %>% as.data.frame()
    mets_res$ace <- do.call('bind_rows',coef_list[c('ace','ae','ce')]) %>% as.data.frame() %>% {rownames(.) <- c('ACE','AE','CE');.}
    mets_res$ade <- do.call('bind_rows',coef_list[c('ade','ae','de')]) %>% as.data.frame() %>% {rownames(.) <- c('ADE','AE','DE');.}
  }
  
  tablist <- list(
    ICC_res=ICC_res,
    use_seed=seed,
    twin.index_res=twin.index_res,
    strataOR_res=strataOR_res,
    GEE_res=GEE_res,
    glm_res=glm_res,
    mets_res=mets_res
  )
  
  knit_report(system.file('Rmd/report.Rmd',package = "twinstat"), filename)
}

