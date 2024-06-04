Args<-commandArgs(TRUE)
proteinname=Args[1]
MRversion=Args[2]
threshold=as.numeric(Args[3])
modality=Args[4]
matching_files <- list.files('/Project/UKB/UKB_pGWAS_annoProtein/', full.names = FALSE, recursive = FALSE)
proteinfile <- matching_files[grepl(paste0("^", paste(proteinname,'_',sep=''),sep=''), matching_files)]
diseasePath <- '/Project/UKB/disease/GWAS/'
diseases <- c('anxiety','adhd','BP','MDD','ASD','SCZ','PD','AD','ALS','MS')
diseaseGWASs <- c('anxiety.meta.full.fs.tbl','daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta','daner_PGC_BIP32b_mds7a_0416a','daner_pgc_mdd_meta_w2_no23andMe_rmUKBB','iPSYCHPGC_ASD_Nov2017','PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv','finngen_R9_G6_PARKINSON','finngen_R9_G6_ALZHEIMER', 'finngen_R9_G6_ALS','finngen_R9_G6_MS')
if(length(proteinfile)==1){
outputpath <- paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/GWAS/protein_disease_',modality,'/',MRversion,'/')
outcomepath <- diseasePath
bfile='/Project/UKB/GWAS/LD_ref_1kgenome/EUR'
outs <- diseaseGWASs
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
# https://mrcieu.github.io/ieugwasr/articles/local_ld.html
library('ieugwasr')
exposure_dat_0$rsid <- exposure_dat_0$SNP
exposure_dat_0$trait_id <- rep(proteinname,nrow(exposure_dat_0))
exposure_dat_0$exposure <- rep(proteinname,nrow(exposure_dat_0))
exposure_clumped <- ld_clump(
    exposure_dat_0,
    plink_bin = '/Project/UKB/GWAS/LD_ref_1kgenome/plink',
    bfile = bfile,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p = threshold,
    pop = "EUR"
)
# only keep snps with F statistic > 10
#exposure_clumped$r2 <- (exposure_clumped$beta.exposure^2) / ((exposure_clumped$beta.exposure^2) + exposure_clumped$N * exposure_clumped$se.exposure^2)
#exposure_clumped$F <- exposure_clumped$r2 * (exposure_clumped$N - 2) / (1 - exposure_clumped$r2)
#exposure_clumped <- exposure_clumped[exposure_clumped$F>10,]

for(i in seq(1,length(outs))){
    outcomefilename <- outs[i]
    outcome <- diseases[i]
    ## 1.Data sets preparation
    # Require two GWAS summary-level data with variants minor allele frequency (MAF) > 0.01; Removing the palindromic SNPs; Removing the long-range LD in the genome. (i.e. exposure_dat_0; outcome_dat_0)
    # can't done with plink, plink work for genotyping data, not summary data
    #########################################  2.MR Analysis preparation #########################################
    # a)Selection of conditional independent SNPs for exposures
   
    outcomefile <- paste0(outcomepath, outcomefilename)
    if (outcome %in% c('ALS','AD','PD','MS')){
        outcomedata <- read.table(outcomefile,sep='\t',comment.char = "@", header=TRUE)
    }else{
    	outcomedata <- read.table(outcomefile,sep='\t',comment.char = "#", header=TRUE)
    }

    snps=exposure_clumped$SNP
    if(outcome %in% c('PD','AD','ALS','MS')){
	    outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "rsids",
	        beta_col = "beta",
	        se_col = "sebeta",
	        effect_allele_col = "alt",
	        other_allele_col = "ref",
	        eaf_col = "af_alt",
	        pval_col = "pval"
	    )
    }else if(outcome %in% c('anxiety')){
    	outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "SNPID",
	        beta_col = "Effect",
	        se_col = "StdErr",
	        effect_allele_col = "Allele1",
	        other_allele_col = "Allele2",
	        eaf_col = "Freq1",
	        pval_col = "P.value"
	    )
    }else if(outcome %in% c('ASD','adhd','BP','MDD')){
    	outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "SNP",
	        beta_col = "OR",
	        se_col = "SE",
	        effect_allele_col = "A1",
	        other_allele_col = "A2",
	        pval_col = "P"
	    )
	    outcome_dat_0$beta.outcome <- log(outcome_dat_0$beta.outcome)
    }else if(outcome %in% c('SCZ')){
    	outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "ID",
	        beta_col = "BETA",
	        se_col = "SE",
	        effect_allele_col = "A1",
	        other_allele_col = "A2",
	        eaf_col = "FCAS",
	        pval_col = "PVAL"
	    )
    }
    outcome_dat_0$outcome <- rep(outcome, nrow(outcome_dat_0))
    # d) Quality control of instrument variants
    # Heterogeneity test: Cochran’s Q test and Rucker’s Q’ test;
    # F-statistics: greater than 10
    # Rscript RadialMR_v1.0_.R Input Output
    output1<-paste0(outputpath, outcome,'_',proteinname,'_QC_IVW') #output data IVW
    output2<-paste0(outputpath, outcome,'_',proteinname,'_QC_MR-egger') #output data MR-egger
    output3<-paste0(outputpath, outcome,'_',proteinname,'_QC_outlierSNPs') #output data outlier snp
    library(RadialMR)
    dat <- harmonise_data(exposure_dat = exposure_clumped,outcome_dat = outcome_dat_0,action=2)
    data3<-format_radial(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome,dat$SNP)
    res1<-ivw_radial(data3,0.05,1,0.0001)
    res2<-egger_radial(data3,0.05,1)
    write.table(as.data.frame(res1$data),output1,quote=F,row.names=F,col.names=T,sep="\t")
    write.table(as.data.frame(res2$data),output2,quote=F,row.names=F,col.names=T,sep="\t")
    if (class(res1$outliers)=='data.frame'){
        a1<-res1$outliers$SNP
    }else {
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
    datoutput<-paste0(outputpath, outcome,'_',proteinname,'_dat')
    write.table(dat, datoutput, quote=F, row.names=F,col.names=T,sep="\t")
}
}
