
Args<-commandArgs(TRUE)
proteinname=Args[1]
print(proteinname)
modality=Args[2]
threshold=as.numeric(Args[3])
threshold_folder=paste0('MR_',Args[4])
matching_files <- list.files('UKB_pGWAS_annoProtein/', full.names = FALSE, recursive = FALSE)
proteinfile <- matching_files[grepl(paste0("^", paste(proteinname,'_',sep=''),sep=''), matching_files)]
allAssocaitions <- readRDS("allAssociations_allModality.rds")
if(length(proteinfile)!=1){
stop('protein GWAS missing')
}else{
bfile='/Project/UKB/GWAS/LD_ref_1kgenome/EUR'
library(TwoSampleMR)
exposure_dat_0 <- read_exposure_data(
    filename = paste0('/Project/UKB/UKB_pGWAS_annoProtein/', proteinfile),
    sep = " ",
    snp_col = "rsid",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "A1FREQ",
    pval_col = "LOG10P",
    samplesize_col = "N",
    log_pval = TRUE
)
exposure_dat_0$exposure <- rep(proteinname,nrow(exposure_dat_0))
# https://mrcieu.github.io/ieugwasr/articles/local_ld.html
library('ieugwasr')
exposure_dat_0$rsid <- exposure_dat_0$SNP
exposure_dat_0$trait_id <- rep(proteinname,nrow(exposure_dat_0))
exposure_clumped <- ld_clump(
    exposure_dat_0,
    plink_bin = '/Project/UKB/GWAS/LD_ref_1kgenome/plink',
    bfile = bfile,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p = threshold,
    pop = "EUR"
)
if(modality=='volume'){
    sampleSize_outcome=25576
    regionNames <- read.table('/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_regions_uniform.csv',sep=',',header=T)
    regionNames <- regionNames[,1]
    outputpath=paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/GWAS/crossModality/volume/',threshold_folder,'/')
    dir.create(outputpath, recursive=TRUE)
    outcomepath <- '/Project/UKB/GWAS/volume/GWAS/GWAS_summary/'
    combines <- allAssocaitions[['volume']]
    library(tidyr)
    combines <- data.frame(combines) %>% separate('combines',into=c('protein', 'metric'),sep='\\+')
    outs <- sapply(combines[combines$protein==proteinname, 'metric'],function(x) paste0('GWAS_volume.',x,'.glm.linear'))
    names(outs) <- NULL
}else if(modality=='thick'){
    sampleSize_outcome=25576
    regionNames <- read.table('/Project/data/imputated_median_scaled_bl_intersection_LRthick_aparc_original_regions_uniform.csv',sep=',',header=T)
    regionNames <- regionNames[,1]
    outputpath=paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/GWAS/crossModality/thickness/',threshold_folder,'/')
    dir.create(outputpath, recursive=TRUE)
    outcomepath <- '/Project/UKB/GWAS/thickness/GWAS/GWAS_summary/'
    combines <- allAssocaitions[['thick']]
    library(tidyr)
    combines <- data.frame(combines) %>% separate('combines',into=c('protein', 'metric'),sep='\\+')
    outs <- sapply(combines[combines$protein==proteinname, 'metric'],function(x) paste0('GWAS_thickness.',x,'.glm.linear'))
    names(outs) <- NULL
}else if (modality=='area'){
    sampleSize_outcome=25576
    regionNames <- read.table('/Project/data/imputated_median_scaled_bl_intersection_LRarea_aparc_original_regions_uniform.csv',sep=',',header=T)
    regionNames <- regionNames[,1]
    outputpath=paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/GWAS/crossModality/area/',threshold_folder,'/')
    dir.create(outputpath, recursive=TRUE)
    outcomepath <- '/Project/UKB/GWAS/Area/GWAS/GWAS_summary/'
    combines <- allAssocaitions[['area']]
    library(tidyr)
    combines <- data.frame(combines) %>% separate('combines',into=c('protein', 'metric'),sep='\\+')
    outs <- sapply(combines[combines$protein==proteinname, 'metric'],function(x) paste0('GWAS_area.',x,'.glm.linear'))
    names(outs) <- NULL
}else if(modality=='DTIFA'){
    sampleSize_outcome=24544
    regionNames <- sapply(1:27,function(x) paste0('FA', as.character(x)))
    outputpath=paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/GWAS/crossModality/FA/',threshold_folder,'/')
    dir.create(outputpath, recursive=TRUE)
    outcomepath <- '/Project/UKB/GWAS/DTI/GWAS/FA/GWAS_summary/'
    combines <- allAssocaitions[['DTIFA']]
    library(tidyr)
    combines <- data.frame(combines) %>% separate('combines',into=c('protein', 'metric'),sep='\\+')
    outs <- sapply(combines[combines$protein==proteinname, 'metric'],function(x) paste0('GWAS_LRFA.',x,'.glm.linear'))
    names(outs) <- NULL
}else if(modality=='DTIMD'){
    sampleSize_outcome=24544
    regionNames <- sapply(1:27,function(x) paste0('MD', as.character(x)))
    outputpath=paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/GWAS/crossModality/MD/',threshold_folder,'/')
    dir.create(outputpath, recursive=TRUE)
    outcomepath <- '/Project/UKB/GWAS/DTI/GWAS/MD/GWAS_summary/'
    combines <- allAssocaitions[['DTIMD']]
    library(tidyr)
    combines <- data.frame(combines) %>% separate('combines',into=c('protein', 'metric'),sep='\\+')
    outs <- sapply(combines[combines$protein==proteinname, 'metric'],function(x) paste0('GWAS_LRMD.',x,'.glm.linear'))
    names(outs) <- NULL
}

for(out in outs){
    outcomefilename=out
    outcome <- strsplit(out, '\\.')[[1]][2]
    outcome <- regionNames[as.numeric(gsub('PHENO','',outcome))]
    #########################################  2.MR Analysis preparation #########################################
    outcomefile <- paste0(outcomepath, outcomefilename)

    snps=exposure_clumped$SNP
    outcome_dat_0 <- read_outcome_data(
        snps = snps,
        filename = outcomefile,
        sep = "\t",
        snp_col = "ID",
        beta_col = "BETA",
        se_col = "SE",
        effect_allele_col = "A1",
        other_allele_col = "OMITTED",
        eaf_col = "A1_FREQ",
        pval_col = "P"
    )
    outcome_dat_0$outcome <- rep(outcome,nrow(outcome_dat_0))
    outcome_dat_0$samplesize.outcome <- rep(sampleSize_outcome,nrow(outcome_dat_0))

    # d) Quality control of instrument variants
    output1<-paste0(outputpath, outcomefilename,'_',proteinname,'_QC_IVW') #output data IVW
    output2<-paste0(outputpath, outcomefilename,'_',proteinname,'_QC_MR-egger') #output data MR-egger
    output3<-paste0(outputpath, outcomefilename,'_',proteinname,'_QC_outlierSNPs') #output data outlier snp
    
    library(RadialMR)
    dat <- harmonise_data(exposure_dat = exposure_clumped,outcome_dat = outcome_dat_0,action=2)
    data3<-format_radial(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome,dat$SNP)
    res1<-ivw_radial(data3,0.05,1,0.0001)
    res2<-egger_radial(data3,0.05,1)
    write.table(as.data.frame(res1$data),output1,quote=F,row.names=F,col.names=T,sep="\t")
    write.table(as.data.frame(res2$data),output2,quote=F,row.names=F,col.names=T,sep="\t")
    if (class(res1$outliers)=='data.frame'){
        a1<-res1$outliers$SNP
    }else{
        a1<-c('None')
    }
    if (class(res2$outliers)=='data.frame'){
        a2<-res2$outliers$SNP
    }else{
        a2<-c('None')
    }
    a3<-unique(c(a1,a2))
    write.table(data.frame(a3),output3,quote=F,row.names=F,col.names=F,sep="\t")


    ######################################### 3. MR analysis #########################################
    output1 <- paste0(outputpath, outcome,'_',proteinname,'_MR_WM') #weight median
    output2 <- paste0(outputpath, outcome,'_',proteinname,'_MR_IVW') #ivw
    output3 <- paste0(outputpath, outcome,'_',proteinname,'_MR_egger') #egger slope
    output4 <- paste0(outputpath, outcome,'_',proteinname,'_MR_W') #weight mode
    output5 <- paste0(outputpath, outcome,'_',proteinname,'_MR_robustAdjProfile') #Robust adjusted profile score
    output6 <- paste0(outputpath, outcome,'_',proteinname,'_MR_waldRatio') #Wald Ratio

    library(TwoSampleMR)
    # remove outlier snps after quality control
    dat <- dat[!(dat$SNP %in% a3),]
    tsmr1<-mr(dat, method_list=c("mr_weighted_median"))
    tsmr2<-mr(dat, method_list=c("mr_ivw"))
    tsmr5<-mr(dat, method_list=c("mr_raps"))
    tsmr6<-mr(dat, method_list=c("mr_wald_ratio"))
    tsmr3<-mr(dat, method_list=c("mr_egger_regression"))
    tsmr4<-mr(dat, method_list=c("mr_weighted_mode"))
    a1<-cbind(outcome, proteinname,tsmr1$nsnp,tsmr1$b,(tsmr1$b-1.96*tsmr1$se),(tsmr1$b+1.96*tsmr1$se),tsmr1$pval) #weight median
    a2<-cbind(outcome, proteinname,tsmr2$nsnp,tsmr2$b,(tsmr2$b-1.96*tsmr2$se),(tsmr2$b+1.96*tsmr2$se),tsmr2$pval) #ivw
    a5<-cbind(outcome, proteinname,tsmr5$nsnp,tsmr5$b,(tsmr5$b-1.96*tsmr5$se),(tsmr5$b+1.96*tsmr5$se),tsmr5$pval) #raps
    a6<-cbind(outcome, proteinname,tsmr6$nsnp,tsmr6$b,(tsmr6$b-1.96*tsmr6$se),(tsmr6$b+1.96*tsmr6$se),tsmr6$pval) #WR
    CI <- 0.95
    lowerCI <- function(beta,df,SE){
    return(beta - (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
    }
    upperCI <- function(beta,df,SE){
    return(beta + (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
    }
    a3 <- cbind(outcome, proteinname,tsmr3$nsnp,tsmr3$b, mapply(lowerCI, tsmr3$b, tsmr3$nsnp - 2, tsmr3$se), mapply(upperCI, tsmr3$b, tsmr3$nsnp - 2, tsmr3$se), tsmr3$pval) #egger slope
    a4 <- cbind(outcome, proteinname,tsmr4$nsnp,tsmr4$b, mapply(lowerCI, tsmr4$b, tsmr4$nsnp - 1, tsmr4$se), mapply(upperCI, tsmr4$b, tsmr4$nsnp - 1, tsmr4$se), tsmr4$pval) #weight mode
    write.table(a1,output1,quote=F,row.names=F,col.names=F,sep="\t") #weight median
    write.table(a2,output2,quote=F,row.names=F,col.names=F,sep="\t") #ivw
    write.table(a5,output5,quote=F,row.names=F,col.names=F,sep="\t") #raps
    write.table(a6,output6,quote=F,row.names=F,col.names=F,sep="\t") #wr
    write.table(a3,output3,quote=F,row.names=F,col.names=F,sep="\t") #egger slope
    write.table(a4,output4,quote=F,row.names=F,col.names=F,sep="\t") #weight mode

    ######################################## 4. sensitivity analysis #########################################
    output1<-paste0(outputpath, outcome,'_',proteinname,'_dat')
    write.table(dat, output1, quote=F, row.names=F,col.names=T,sep="\t")
    library(MRPRESSO)
    output1<-paste0(outputpath, outcome,'_',proteinname,'_sensitivity_MR_PRESSO_distortion') #distortion_test
    output2<-paste0(outputpath, outcome,'_',proteinname,'_sensitivity_MR_PRESSO_global') #global_test
    output3<-paste0(outputpath, outcome,'_',proteinname,'_sensitivity_MR_PRESSO_main') #main_MR_results
    output4<-paste0(outputpath, outcome,'_',proteinname,'_sensitivity_MR_PRESSO_outlierTest') #eachsnp.outlier_test
    bzx<-dat$beta.exposure
    sebzx<-dat$se.exposure
    bzy<-dat$beta.outcome
    sebzy<-dat$se.outcome
    data3<-as.data.frame(cbind(bzx,sebzx,bzy,sebzy))
    results<-mr_presso(BetaOutcome = "bzy", BetaExposure = "bzx", SdOutcome= "sebzy", SdExposure = "sebzx", OUTLIERtest = TRUE, DISTORTIONtest= TRUE, data = data3, NbDistribution = 10000, SignifThreshold = 0.05)
    a1<-results$`MR-PRESSO results`$`Distortion Test` #testing of significant distortion in the causal estimates before and after outlier removal
    a2<-results$`Main MR results`
    a3<-results$`MR-PRESSO results`$`Global Test`$Pvalue #detection of horizontal pleiotropy
    a4<-results$`MR-PRESSO results`$`Global Test`$RSSobs
    a5<-results$`MR-PRESSO results`$`Outlier Test` #correction of horizontal pleiotropy via outlier removal
    b2<-cbind(a4,a3)
    a2A<-cbind(a2[2],a2[3],a2[4],a2[5],a2[6])
    write.table(a1,output1,quote=F,row.names=F,col.names=F,sep="\t")
    write.table(b2,output2,quote=F,row.names=F,col.names=F,sep="\t")
    write.table(a2A,output3,quote=F,row.names=T,col.names=T,sep="\t")
    write.table(a5,output4,quote=F,row.names=F,col.names=T,sep="\t")

    output1 <- paste0(outputpath, outcome,'_',proteinname,'_sensitivity_MR_egger_intercept')
    Inter <- mr_egger_regression(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
    out <- cbind(Inter$b_i, Inter$se_i, Inter$pval_i) #egger intercept
    write.table(out,output1,quote=F,row.names=F,col.names=F,sep="\t") #egger intercept

    output1 <- paste0(outputpath, outcome,'_',proteinname,'_sensitivity_MR_LOO')
    loo <- mr_leaveoneout(dat)
    b <- data.frame(loo$SNP,loo$b,loo$se,loo$p)
    write.table(b,output1,quote=F,row.names=F,col.names=T,sep="\t") #leave one out results
    
}
}
