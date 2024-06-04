# colocalizaion for protein-imaging MR
library("coloc")
library(dplyr)
library(data.table)
library(tidyr)
Args<-commandArgs(TRUE)
MRversion=Args[1]
threshold=as.numeric(Args[2])
r2=as.numeric(Args[3])
print(paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Colocalization/pQTL_BraineQTL/pQTL_CortexeQTL_colocalization_',MRversion,'_',as.character(r2),'.csv'))

# input  GWAS data (Brain eQTL), need to align the effect allele between exposure and outcome manually, we will do this later in this script

eQTLframe_annotated <- fread("/Project/GTEx/Brain_Cortex_v8_SMR/Brain_Cortex_all_withMAF.txt",sep='\t')
fuma_mapped <- readRDS("/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/eQTL/FUMA_mapped_eQTL_GTEx_Cortex.RData")
used_proteins <- names(fuma_mapped)
finaloutput <- data.frame(matrix(NA,0,11))
colnames(finaloutput) <- c('exposure','outcome','V1','nsnps','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf','Best_causal_variant','SNP.PP.H4')

library(foreach)
library(doParallel)
# Perform linear regression for each combination of columns using foreach
num_cores <- detectCores()
registerDoParallel(60)  # Adjust the number of cores as needed
foreach(s=1:length(used_proteins)) %dopar% {
	s=used_proteins[s]
	proteinname=s
	print(proteinname)
	matching_files <- list.files('/Project/UKB/UKB_pGWAS_annoProtein/', full.names = FALSE, recursive = FALSE)
	proteinfile <- matching_files[grepl(paste0("^", paste(proteinname,'_',sep=''),sep=''), matching_files)]
	exposure<-fread(paste0("/Project/UKB/UKB_pGWAS_annoProtein/",proteinfile),header=T)
	genome_build_exposure='hg19'
	exposure$pval_nominal <- 10^(-exposure$LOG10P)
	names(exposure)[2]<-"CHR"
	names(exposure)[18]<-"BP"
	names(exposure)[17]<-"SNP"
	names(exposure)[16]<-'ALT'
	names(exposure)[15]<-'REF'
	names(exposure)[6]<-'A1freq'
	names(exposure)[10]<-"beta"
	names(exposure)[11]<-"SE"
	names(exposure)[20]<-"pval_nominal"
	exposure$varbeta<-exposure$SE^2
	exposure$MAF=ifelse(exposure$A1freq<0.5, exposure$A1freq,1- exposure$A1freq)
	exposure <- exposure[,c("SNP","CHR","BP","ALT","REF","MAF","beta","SE","pval_nominal","varbeta")]
	exposure$pval_nominal<-as.numeric(exposure$pval_nominal)

	for(j in seq(1,length(fuma_mapped[[s]][['gene']]))){
			m <- fuma_mapped[[s]][['gene']][j]
			outcome_genename <- fuma_mapped[[s]][['symbol']][j]
			print(outcome_genename)
			outcome <- eQTLframe_annotated[eQTLframe_annotated$Probe==m,]
			genome_build_outcome='hg19'
			names(outcome)[2]<-"CHR"
			names(outcome)[3]<-"BP"
			names(outcome)[1]<-"SNP"
			names(outcome)[5]<-'REF'
			names(outcome)[4]<-'ALT'
			names(outcome)[20]<-'MAF'
			names(outcome)[12]<-"beta"
			names(outcome)[13]<-"SE"
			names(outcome)[14]<-"pval_nominal"
			outcome <- outcome[,c('CHR','BP','SNP','REF','ALT','MAF','beta','SE','pval_nominal')]
			outcome$varbeta<-outcome$SE^2
			outcome <- na.omit(outcome)
			if(nrow(outcome)!=0){
					outcome <- outcome[,c('CHR',"SNP","BP","ALT","REF","MAF","beta","SE","pval_nominal","varbeta")]
					outcome <- outcome[!duplicated(outcome),]
					outcome <- outcome[order(outcome$pval_nominal),]
					outcome <- distinct(outcome, SNP, .keep_all = TRUE)

					commonSNP <- intersect(exposure$SNP, outcome$SNP)
					lead_exposure <- exposure[exposure$SNP %in% commonSNP,]
					outcome <- outcome[outcome$SNP %in% commonSNP,]

					library(TwoSampleMR)
					# get lead snps for the exposure
					lead_exposure <-format_data(
					  as.data.frame(lead_exposure),
					  pos = "BP",
					  chr = "CHR",
					  type='exposure',
					  snp_col = "SNP",
					  beta_col = "beta",
					  se_col = "SE",
					  effect_allele_col ="ALT",
					  other_allele_col = "REF",
					  pval_col = "pval_nominal"
					)
					library('ieugwasr')
					lead_exposure$rsid <- lead_exposure$SNP
					lead_exposure$pval <- lead_exposure$pval.exposure
					lead_exposure$trait_id <- rep(proteinname,nrow(lead_exposure))
					lead_exposure$exposure <- rep(proteinname,nrow(lead_exposure))
					lead_exposure_clumped <- ld_clump(
					    lead_exposure,
					    plink_bin = '/Project/UKB/GWAS/LD_ref_1kgenome/plink',
					    bfile = '/Project/UKB/GWAS/LD_ref_1kgenome/EUR',
					    clump_kb = 1000,
					    clump_r2 = r2,
					    clump_p = threshold,
					    pop = "EUR"
					)

					lead_outcome <-format_data(
					  as.data.frame(outcome),
					  pos = "BP",
					  chr = "CHR",
					  type='outcome',
					  snp_col = "SNP",
					  beta_col = "beta",
					  se_col = "SE",
					  effect_allele_col ="ALT",
					  other_allele_col = "REF",
					  pval_col = "pval_nominal"
					)
					library('ieugwasr')
					lead_outcome$rsid <- lead_outcome$SNP
					lead_outcome$pval <- lead_outcome$pval.outcome
					lead_outcome$trait_id <- rep(outcome_genename,nrow(lead_outcome))
					lead_outcome$outcome <- rep(outcome_genename,nrow(lead_outcome))
					lead_outcome_clumped <- tryCatch({
						ld_clump(
						    lead_outcome,
						    plink_bin = '/Project/UKB/GWAS/LD_ref_1kgenome/plink',
						    bfile = '/Project/UKB/GWAS/LD_ref_1kgenome/EUR',
						    clump_kb = 1000,
						    clump_r2 = r2,
						    clump_p = threshold,
						    pop = "EUR")
						},
						warning = function(w) {
							  # Assign a placeholder value to the model if you want to continue with further processing
							  # Handle warning
							  cat("Warning at clumping outcome ", as.character(j), ":", conditionMessage(w), "\n")
							  NULL
						}, 
						error = function(e) {
							  # Assign a placeholder value to the model if you want to continue with further processing
							  # Handle error
							  cat("Error at clumping outcome ", as.character(j), ":", conditionMessage(e), "\n")
							  NULL
						})

					# get lead snp window for lead SNPs
					# filter the region 500kb around the most significant SNP  of exposure, build 19/37
					windowframe <- list()
					if(exists('lead_exposure_clumped') & !is.null(lead_exposure_clumped)){
						lead_exposure_clumped <- as.data.frame(lead_exposure_clumped)
						for(i in seq(1,nrow(lead_exposure_clumped))){
							leadSNP_chr <- as.numeric(lead_exposure_clumped[i,'chr.exposure'])
							leadSNP_BP <- as.numeric(lead_exposure_clumped[i,'pos.exposure'])
							lead_regions <- c(leadSNP_BP-500000, leadSNP_BP+500000)
							windowframe[[length(windowframe)+1]] <- c(leadSNP_chr, lead_regions[1], lead_regions[2])
						}
					}
					if(exists('lead_outcome_clumped') & !is.null(lead_outcome_clumped)){
						lead_outcome_clumped <- as.data.frame(lead_outcome_clumped)
						for(i in seq(1,nrow(lead_outcome_clumped))){
							leadSNP_chr <- as.numeric(lead_outcome_clumped[i,'chr.outcome'])
							leadSNP_BP <- as.numeric(lead_outcome_clumped[i,'pos.outcome'])
							lead_regions <- c(leadSNP_BP-500000, leadSNP_BP+500000)
							windowframe[[length(windowframe)+1]] <- c(leadSNP_chr, lead_regions[1], lead_regions[2])
						}
					}
					Results_pleio=list()
					if(length(windowframe)!=0){
						windowframe <- as.data.frame(do.call(rbind, windowframe))
						colnames(windowframe) <- c('chrom','start','end')
						windowframe <- windowframe[!duplicated(windowframe),]
						for(i in seq(1,nrow(windowframe))){
							input <- exposure[exposure$CHR==windowframe[i,1] & exposure$BP>windowframe[i,2] & exposure$BP<windowframe[i,3],]
							input <- input[!duplicated(input$SNP),]
							inputcurrent <- merge(input, outcome,by="SNP", suffixes=c("_exposure","_outcome"))
							# align effect allele for exposure and outcome, if the effect allele is not aligned, take the inverse
							if(nrow(inputcurrent)!=0){
								inputcurrent[inputcurrent$ALT_exposure!=inputcurrent$ALT_outcome, 'beta_outcome'] <- -inputcurrent[inputcurrent$ALT_exposure!=inputcurrent$ALT_outcome, 'beta_outcome']
							}
							result <- tryCatch({
							  coloc.abf(dataset1=list(pvalues=inputcurrent$pval_nominal_exposure, beta=inputcurrent$beta_exposure, varbeta=inputcurrent$varbeta_exposure,type="quant", N=34049,snp = inputcurrent$SNP), dataset2=list(pvalues=inputcurrent$pval_nominal_outcome, beta=inputcurrent$beta_outcome, varbeta=inputcurrent$varbeta_outcome, type="quant", N=183, snp=inputcurrent$SNP), MAF=inputcurrent$MAF_exposure)
							}, 
							error = function(e) {
							  # Assign a placeholder value to the model if you want to continue with further processing
							  # Handle error
							  cat("Error at window ", as.character(i), ":", conditionMessage(e), "\n")
							  NULL
							})
							# Check if the model is null (indicating an error or warning occurred)
							if (!is.null(result)) {
							  # Summarize the model
							  Results_pleio$exposure[i]<-proteinname
							  Results_pleio$outcome[i]<-outcome_genename
								Results_pleio$V1[i]<-i
								Results_pleio$nsnps[i]<-result$summary[1]
								Results_pleio$PP.H0.abf[i]<-result$summary[2]
								Results_pleio$PP.H1.abf[i]<-result$summary[3]
								Results_pleio$PP.H2.abf[i]<-result$summary[4]
								Results_pleio$PP.H3.abf[i]<-result$summary[5]
								Results_pleio$PP.H4.abf[i]<-result$summary[6]
								Results_pleio$Best_causal_variant[i]<- result$results$snp[which(result$results$SNP.PP.H4==max(result$results$SNP.PP.H4))]
								Results_pleio$SNP.PP.H4[i]<- max(result$results$SNP.PP.H4)
								# plotting locuscompare plot

							} else {
							  # Handle the case where the model is null
							  Results_pleio$exposure[i]<-proteinname
							  Results_pleio$outcome[i]<-outcome_genename
								Results_pleio$V1[i]<-i
								Results_pleio$nsnps[i]<-NA
								Results_pleio$PP.H0.abf[i]<-NA
								Results_pleio$PP.H1.abf[i]<-NA
								Results_pleio$PP.H2.abf[i]<-NA
								Results_pleio$PP.H3.abf[i]<-NA
								Results_pleio$PP.H4.abf[i]<-NA
								Results_pleio$Best_causal_variant[i]<- NA
								Results_pleio$SNP.PP.H4[i]<- NA  
							}
						}
					}else{
								Results_pleio$exposure[1]<-proteinname
							  Results_pleio$outcome[1]<-outcome_genename
								Results_pleio$V1[1]<-1
								Results_pleio$nsnps[1]<-NA
								Results_pleio$PP.H0.abf[1]<-NA
								Results_pleio$PP.H1.abf[1]<-NA
								Results_pleio$PP.H2.abf[1]<-NA
								Results_pleio$PP.H3.abf[1]<-NA
								Results_pleio$PP.H4.abf[1]<-NA
								Results_pleio$Best_causal_variant[1]<- NA
								Results_pleio$SNP.PP.H4[1]<- NA  
					}
			}else{
			          Results_pleio=list()
								Results_pleio$exposure[1]<-proteinname
							  Results_pleio$outcome[1]<-outcome_genename
								Results_pleio$V1[1]<-1
								Results_pleio$nsnps[1]<-'None'
								Results_pleio$PP.H0.abf[1]<-'None'
								Results_pleio$PP.H1.abf[1]<-'None'
								Results_pleio$PP.H2.abf[1]<-'None'
								Results_pleio$PP.H3.abf[1]<-'None'
								Results_pleio$PP.H4.abf[1]<-'None'
								Results_pleio$Best_causal_variant[1]<- 'None'
								Results_pleio$SNP.PP.H4[1]<- 'None'  
			}
			temp <- t(do.call(rbind, Results_pleio))
			finaloutput <- rbind(finaloutput, temp)
	}
	
	write.table(finaloutput, paste0('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Colocalization/pQTL_BraineQTL/',proteinname,'_pQTL_CortexeQTL_colocalization_',MRversion,'_',as.character(r2),'_GTExv8.csv'), sep=',', row.names=F)

}

