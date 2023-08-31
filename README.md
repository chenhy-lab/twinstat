# twinstat:Analysis of Twins in R <img src="inst/figures/imgfile.png" align="right" height="138"/>

## :writing_hand: Author

Huiyao Chen

------------------------------------------------------------------------

## :arrow_double_down: Installation

Install the twinstat package via the Github repository:

``` r
# install.package("remotes")   #In case you have not installed it.
remotes::install_github("chenhy-lab/twinstat")
```

------------------------------------------------------------------------

## Examples

> `load_example_data` function will loads example data for testing and learning purposes
> printing according to <https://ysph.yale.edu/c2s2/software/twin-analysis/real-data-set/>

### Demo script

``` r
library(twinstat)
data <- load_example_data()
# 1 is homozygous, 0 is dizygous
twin.chisqTest(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
chisq_res <- twin.chisqTest(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
stat_mat <- chisq_res$stat_mat

## Public correlation models
# Estimate the intraclass correlation
get_ICC(n_neither=stat_mat['MZ','Neither'],
        n_one=stat_mat['MZ','One'],
        n_both=stat_mat['MZ','Both'])
get_ICC(n_neither=stat_mat['DZ','Neither'],
        n_one=stat_mat['DZ','One'],
        n_both=stat_mat['DZ','Both'])
ICC_res <- get_ICC(n_neither=stat_mat['MZ','Neither'],
                   n_one=stat_mat['MZ','One'],
                   n_both=stat_mat['MZ','Both'])
# Compare intraclass correlation between two groups
compareICC(stat_mat)


## Stratified ratio methods
# StratificationOR
twin.index_res <- twin.index(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
compareMD(twin.index_res)

# Generalized estimating equations (GEE)
GEE_res <- GEE(data=data,famid='FAMID',outcome='BPD',zyg='zygo')

## logistic regression
glm_res <- logisticMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo')

## mets package
nlmMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo',cov_var='RDS',use_model='bp')
```


## Related Tools

- [TwinAnalysis](https://github.com/IvanVoronin/TwinAnalysis/): A package for 
  structural equation modeling and twin analysis using OpenMx package for R.
- [mets](https://github.com/kkholst/mets): Implementation of 
	various statistical models for multivariate event history data.
- [TwinAnalysisR](https://github.com/SherryDong/TwinAnalysisR): 
  A function for ACE model for twin analysis.
