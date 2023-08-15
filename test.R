library(dplyr)
library(twinstat)
# 1 is homozygous, 0 is dizygous
data <- load_example_data()
twin.chisqTest(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
chisq_res <- twin.chisqTest(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
stat_mat <- chisq_res$stat_mat

get_ICC(n_neither=stat_mat['MZ','Neither'],
        n_one=stat_mat['MZ','One'],
        n_both=stat_mat['MZ','Both'])
get_ICC(n_neither=stat_mat['DZ','Neither'],
        n_one=stat_mat['DZ','One'],
        n_both=stat_mat['DZ','Both'])
ICC_res <- get_ICC(n_neither=stat_mat['MZ','Neither'],
                   n_one=stat_mat['MZ','One'],
                   n_both=stat_mat['MZ','Both'])

compareICC(stat_mat)


# StratificationOR
twin.index_res <- twin.index(data=data,famid='FAMID',outcome='BPD',zyg='zygo')
compareMD(twin.index_res)
compareMD_res <- compareMD(twin.index_res)

# Generalized estimating equations (GEE)
GEE_res <- GEE(data=data,famid='FAMID',outcome='BPD',zyg='zygo')

# logistic regression
glm_res <- logisticMD(data=data,famid='FAMID',outcome='BPD',zyg='zygo')





