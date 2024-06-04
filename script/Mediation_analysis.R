library('mediation')
library('survival')
library('data.table')
library(tidyr)
library(dplyr)

inverseNorm_transformation <- function(data){
  library(RNOmni)
  return(apply(data,2,function(x)RankNorm(x)))
}

getSigMetric <- function(protein,modality,pair){
  # get the more significant one for a pair of regions
  if(modality=='volume'){
    data <- readRDS("/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/imputated_median_scaled_bl_intersection_LRvolume_aparc_allremoved_result.rds")
    regions <- read.table("/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_regions.csv",sep=',',header=T)
  }else if(modality=='area'){
    data <- readRDS("/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/imputated_median_scaled_bl_intersection_LRarea_aparc_allremoved_result.rds")
    regions <- read.table("/Project/data/imputated_median_scaled_bl_intersection_LRarea_aparc_original_regions.csv",sep=',',header=T)
  }else if(modality=='thick'){
    data <- readRDS("/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/imputated_median_scaled_bl_intersection_LRthick_aparc_allremoved_result.rds")
    regions <- read.table("/Project/data/imputated_median_scaled_bl_intersection_LRthick_aparc_original_regions.csv",sep=',',header=T)
  }
  regions <- as.character(regions[1,2:ncol(regions)])
  index <- match(protein,data$proteinNames)
  used_regions <- data$regions_sig_oneProtein_overall_BH[[index]]
  used_regions_ID <- sapply(match(used_regions, regions),function(x) paste0('R',as.character(x)))
  used_pvalues <- data$regions_sig_oneProtein_overall_BH_pvalue[[index]]
  pvalueframe <- data.frame(used_regions_ID, used_pvalues)
  rownames(pvalueframe) <- pvalueframe$used_regions_ID
  pair <- as.character(pair)
  pair <- sapply(pair,function(x) paste0('R',x))
  pvalueframe <- pvalueframe[pvalueframe$used_regions_ID %in% pair,]
  return(pvalueframe[pvalueframe$used_pvalues==min(pvalueframe$used_pvalues),'used_regions_ID'])
}

getHighCorr <- function(data, thres){
  labels <- colnames(data)
  list_output <- list()
  for(i in seq(1,nrow(data))){
    for(j in seq(1,ncol(data))){
      if(i>j & data[i,j]>thres & data[i,j]!=1){
          list_output[[length(list_output)+1]] <- c(labels[i],labels[j])
      }
    }
  }
  return(list_output)
}

getHighestCorr <- function(data, thres){
  labels <- colnames(data)
  list_output <- list()
  data[data==1 | data<thres] <- 0
  if(max(data)!=0){ # no too high colinearity
    temp <- which(data==max(data),arr.ind=T)
      return(rownames(temp))
  }else{
      return(NULL)
  }
}

proteindata <- fread('/Project/data/imputated_median_scaled_bl_INT_AgeSexRemoved_allProteomic.csv')
proteindata <- as.data.frame(proteindata)
proteinName <- read.csv('/Project/data/imputated_median_scaled_bl_proteinNames.csv')
NewProteinname <- sapply(proteinName[,1],function(x) gsub('-','_',x))
colnames(proteindata) <- c('eid',NewProteinname)
basicInfo  <- read.csv('/Project/UKB/cova_data_more_withPC_all2visits.csv')
basicInfo <- basicInfo[!is.na(basicInfo$Smoking_status) & !is.na(basicInfo$Drinker_status) & basicInfo$Smoking_status!=-3 & basicInfo$Drinker_status!=-3,]
basicInfo <- fastDummies::dummy_cols(basicInfo, select_columns = "Drinker_status", remove_first_dummy = TRUE)
basicInfo <- fastDummies::dummy_cols(basicInfo, select_columns = "Smoking_status", remove_first_dummy = TRUE)
basicInfo <- basicInfo[,c('eid','interval','Sex','Age','eTIV','Site1','Site2','Education','Townsend_index','Smoking_status_1','Smoking_status_2','Drinker_status_1','Drinker_status_2')]

# metric data
volumedata <- fread('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/association_finalvolumeData.csv')
volumedata <- as.data.frame(volumedata)
thickdata <- fread('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/association_finalthickData.csv')
thickdata <- as.data.frame(thickdata)
areadata <- fread('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/association_finalareaData.csv')
areadata <- as.data.frame(areadata)
DTIFAdata <- fread('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/association_finalDTIFAData.csv')
DTIFAdata <- as.data.frame(DTIFAdata)
DTIMDdata <- fread('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/association_finalDTIMDData.csv')
DTIMDdata <- as.data.frame(DTIMDdata)
regions <- read.table('/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_regions_uniform.csv',sep=',',header = T)
originalregions <- read.table('/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_regions.csv',sep=',',header = T)
originalregions <- transpose(originalregions[1,2:ncol(originalregions)])
rownames(originalregions) <- regions[,1]
# associated metrics
assoVolume <- read.table('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Associatedmetric_volume_sigMRproteins_crossModality_MR_5e6_BH_loose_details.csv',header = T)
assoVolume <- separate(assoVolume, 'asso_volumes', into = c('metric','protein'),sep='-')
assoThick <- read.table('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Associatedmetric_thick_sigMRproteins_crossModality_MR_5e6_BH_loose_details.csv',header = T)
assoThick <- separate(assoThick, 'asso_volumes', into = c('metric','protein'),sep='-')
assoArea <- read.table('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Associatedmetric_area_sigMRproteins_crossModality_MR_5e6_BH_loose_details.csv',header = T)
assoArea <- separate(assoArea, 'asso_volumes', into = c('metric','protein'),sep='-')
assoFA <- read.table('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Associatedmetric_FA_sigMRproteins_crossModality_MR_5e6_BH_loose_details.csv',header = T)
assoFA <- separate(assoFA, 'asso_volumes', into = c('metric','protein'),sep='-')
assoMD <- read.table('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/Associatedmetric_MD_sigMRproteins_crossModality_MR_5e6_BH_loose_details.csv',header = T)
assoMD <- separate(assoMD, 'asso_volumes', into = c('metric','protein'),sep='-')



options(digits = 20)
options(scipen = 4)
list_cogdata <- list()
list_imagingdata <- list()
# cognition data
cogPath <- '/Project/UKB/cogMental.csv'
allcogs <- read.table(cogPath, header=T, sep=',')
cog_corr <- read.table('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/sigMR_crossModality_MR_5e6_BH_BH_cogCorr_bl_withimaging_loose.csv',sep=',',header=T)
# constrain the protein to those associated with both imaing and cognition
# cog_corr <- cog_corr[cog_corr$Bonf<0.05,c('list_proteins','list_cogdomains')]
used_proteins <- unique(cog_corr$list_proteins)
for(s in used_proteins){
  list_cogdata[[s]] <- colnames(allcogs)[2:ncol(allcogs)]
}
# constrain the imaging metrics to those associated with each protein
for(s in used_proteins){
  list_imagingdata[[s]][['volume']] <- assoVolume[assoVolume$protein==s,'metric']
  list_imagingdata[[s]][['thick']]  <- assoThick[assoThick$protein==s,'metric']
  list_imagingdata[[s]][['area']] <- assoArea[assoArea$protein==s,'metric']
  list_imagingdata[[s]][['FA']] <- assoFA[assoFA$protein==s,'metric']
  list_imagingdata[[s]][['MD']] <- assoMD[assoMD$protein==s,'metric']
}

removeCovariates <- function(X_before,C){
  C <- as.matrix(C)
  X_before <- as.matrix(X_before)
  # Perform linear regression to estimate the coefficients
  B <- solve(t(C) %*% C) %*% t(C) %*% X_before  # Coefficients matrix
  # Predict the covariate effect
  C_effect <- C %*% B
  # Remove the effect of covariates from each column of X
  X_corrected <- X_before - C_effect
  return(X_corrected)
}


library(foreach)
library(doParallel)
# Perform linear regression for each combination of columns using foreach
num_cores <- detectCores()
registerDoParallel(100)  # Adjust the number of cores as needed

cognitionframe <- data.frame(matrix(nrow = 0, ncol = 11))
colnames(cognitionframe) <- c('protein','a_coef','a_p','b_coef','b_p','c_coef','c_p','cp_coef','cp_p','ab_coef','ab_p')
mentalframe <- data.frame(matrix(nrow = 0, ncol = 11))
colnames(mentalframe) <- c('protein','a_coef','a_p','b_coef','b_p','c_coef','c_p','cp_coef','cp_p','ab_coef','ab_p')
sink('/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/output_customizedSelection_loose.txt')
for(j in seq(1,length(used_proteins))){
  rm(assometricData)
  # protein loop
  protein <- used_proteins[j]
  usedproteindata <- proteindata[,c(protein,'eid')]
  colnames(usedproteindata) <- c('protein','eid')
  # load associated metrics
  if (length(list_imagingdata[[protein]][['volume']])!=0){
    # volume
    used_volumedata <- volumedata[,c(1,as.numeric(sapply(list_imagingdata[[protein]][['volume']],function(x) gsub('R','',x)))+1)]
    colnames(used_volumedata)[2:ncol(used_volumedata)] <- sapply(list_imagingdata[[protein]][['volume']],function(x) paste0('volume_',x))
    assometricData <- used_volumedata
  }
  if (length(list_imagingdata[[protein]][['thick']])!=0){
    # thick
    used_thickdata <- thickdata[,c(1,as.numeric(sapply(list_imagingdata[[protein]][['thick']],function(x) gsub('R','',x)))+1)]
    colnames(used_thickdata)[2:ncol(used_thickdata)] <- sapply(list_imagingdata[[protein]][['thick']],function(x) paste0('thick_',x))
    if (!exists('assometricData')){
        assometricData <- used_thickdata
    }else{
        assometricData <- merge(assometricData,used_thickdata,by='eid')
    }
  }
 if (length(list_imagingdata[[protein]][['area']])!=0){
    # area
    used_areadata <- areadata[,c(1,as.numeric(sapply(list_imagingdata[[protein]][['area']],function(x) gsub('R','',x)))+1)]
    colnames(used_areadata)[2:ncol(used_areadata)] <- sapply(list_imagingdata[[protein]][['area']],function(x) paste0('area_',x))
    if (!exists('assometricData')){
        assometricData <- used_areadata
    }else{
        assometricData <- merge(assometricData,used_areadata,by='eid')
    }
  }
  if (length(list_imagingdata[[protein]][['FA']])!=0){
    # FA
    used_FAdata <- DTIFAdata[,c(1,as.numeric(sapply(list_imagingdata[[protein]][['FA']],function(x) gsub('R','',x)))+1)]
    colnames(used_FAdata)[2:ncol(used_FAdata)] <- sapply(list_imagingdata[[protein]][['FA']],function(x) paste0('FA_',x))
    if (!exists('assometricData')){
        assometricData <- used_FAdata
    }else{
        assometricData <- merge(assometricData,used_FAdata,by='eid')
    }
  }
  if (length(list_imagingdata[[protein]][['MD']])!=0){
    # MD
    used_MDdata <- DTIMDdata[,c(1,as.numeric(sapply(list_imagingdata[[protein]][['MD']],function(x) gsub('R','',x)))+1)]
    colnames(used_MDdata)[2:ncol(used_MDdata)] <- sapply(list_imagingdata[[protein]][['MD']],function(x) paste0('MD_',x))
    if (!exists('assometricData')){
        assometricData <- used_MDdata
    }else{
        assometricData <- merge(assometricData,used_MDdata,by='eid')
    }
  }
  assometricData[,2:ncol(assometricData)] <- inverseNorm_transformation(data.frame(assometricData[,2:ncol(assometricData)]))


  # handle the colinearity of the observation variables
  modality <- colnames(assometricData)[2]
  modality <- strsplit(modality,'_')[[1]][1]
  if(modality=='volume'){
    newassometricData <- assometricData
    corMat <- cor(newassometricData)
    while(!is.null(getHighestCorr(corMat, 0.55))){
        high_pairs <- getHighestCorr(corMat, 0.55)
        Moresig <- getSigMetric(protein, modality, gsub('volume_R','',high_pairs))
        deleted <- setdiff(high_pairs, paste0('volume_',Moresig))
        newassometricData[deleted] <- NULL
        corMat <- cor(newassometricData)
    }
  }else if(modality=='area'){
    newassometricData <- assometricData
    corMat <- cor(newassometricData)
    while(!is.null(getHighestCorr(corMat, 0.55))){
        high_pairs <- getHighestCorr(corMat, 0.55)
        Moresig <- getSigMetric(protein, modality, gsub('area_R','',high_pairs))
        deleted <- setdiff(high_pairs, paste0('area_',Moresig))
        newassometricData[deleted] <- NULL
        corMat <- cor(newassometricData)
    }
  }else if(modality=='thick'){
    newassometricData <- assometricData
    corMat <- cor(newassometricData)
    while(!is.null(getHighestCorr(corMat, 0.55))){
        high_pairs <- getHighestCorr(corMat, 0.55)
        Moresig <- getSigMetric(protein, modality, gsub('thick_R','',high_pairs))
        deleted <- setdiff(high_pairs, paste0('thick_',Moresig))
        newassometricData[deleted] <- NULL
        corMat <- cor(newassometricData)
    }
  }
  assometricData <- newassometricData

  tempframe <- data.frame(matrix(nrow = 0, ncol = 21))
  colnames(tempframe) <- c('protein','metric','cog','warning','sampleSize','pvalue_total','lowerCI_total','upperCI_total','coef_total','pvalue_indirect','lowerCI_indirect','upperCI_indirect','coef_indirect','pvalue_direct','lowerCI_direct','upperCI_direct','coef_direct','pvalue_propMediate','lowerCI_propMediate','upperCI_propMediate','coef_propMediate')
    warning_model <- 0    
    index_protein <- 2
    usedmetricData <- assometricData
    index_cog <- match(list_cogdata[[protein]][6:length(list_cogdata[[protein]])], colnames(allcogs))
    cogdata <- allcogs[,c(1,index_cog)]
    index_cog <- sapply(seq(1,ncol(cogdata)-1),function(x) x+2)
    cogdata <- merge(usedproteindata, cogdata, by='eid',all.x = T, all.y = T)
    index_basicInfo <- sapply(seq(1,ncol(basicInfo)-1),function(x) x+ncol(cogdata))
    cogdata <- merge(cogdata, basicInfo, by='eid',all.x=T, all.y=T)
    index_meric <- sapply(seq(1,ncol(usedmetricData)-1),function(x) x+ncol(cogdata))
    cogdata <- merge(cogdata, usedmetricData, by='eid',all.x=T, all.y=T)
    tempdata <- cogdata
    cleaned_data <- na.omit(tempdata)
    cleaned_data[,index_cog] <- inverseNorm_transformation(data.frame(cleaned_data[,index_cog]))
    # remove covariates before the plspm model
    cleaned_data[,index_meric] <- removeCovariates(cleaned_data[,index_meric], cleaned_data[,index_basicInfo])
    #cleaned_data[,index_cog] <- removeCovariates(cleaned_data[,index_cog],cleaned_data[,index_basicInfo[c(2,3,7,8,9,10,11,12)]])
    cleaned_data[,index_cog] <- removeCovariates(cleaned_data[,index_cog],cleaned_data[,index_basicInfo])
    colnames(cleaned_data)[index_meric] <- gsub(' ','_',colnames(cleaned_data)[index_meric])
    colnames(cleaned_data)[index_meric] <- gsub('\\(','',colnames(cleaned_data)[index_meric])
    colnames(cleaned_data)[index_meric] <- gsub('\\)','',colnames(cleaned_data)[index_meric])
    library(lavaan)
    metricname <- paste(colnames(cleaned_data)[index_meric], collapse='+')
    metricname <- gsub(' ','_',metricname)
    metricname <- gsub('\\(','',metricname)
    metricname <- gsub('\\)','',metricname)
    # mental health
    model <- paste("cog =~ depress_PHQ9+anxiety+mania+wellbeing+psychotic_experience+selfharm+mental_distress+trauma",paste0("metric =~ ", metricname),"metric ~ a*protein","cog ~ b*metric","cog ~ cp*protein","ab := a * b","total := cp + ab",sep='\n')
    sem_result <- sem(model, cleaned_data, se = "bootstrap", bootstrap = 1000)
    warning_model <- tryCatch(sem(model, cleaned_data, se = "bootstrap", bootstrap = 1000),warning=function(w) if (grepl("the optimizer warns that a solution has NOT been found!", conditionMessage(w))) {1}else{0},finally = {if (!exists("w")) {0}})
    print(protein)
    result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
    if(class(warning_model)=='lavaan'){
      CI <- parameterEstimates(sem_result, boot.ci.type = "bca.simple")
      list_output <- data.frame(protein,CI[CI$label=='a','est'],CI[CI$label=='a','pvalue'],CI[CI$label=='b','est'],CI[CI$label=='b','pvalue'],CI[CI$label=='total','est'],CI[CI$label=='total','pvalue'],CI[CI$label=='cp','est'],CI[CI$label=='cp','pvalue'],CI[CI$label=='ab','est'],CI[CI$label=='ab','pvalue'])
      colnames(list_output) <- c('protein','a_coef','a_p','b_coef','b_p','c_coef','c_p','cp_coef','cp_p','ab_coef','ab_p')
      mentalframe <- rbind(mentalframe, list_output)
      print(result)
      print(CI)
    }else if(is.numeric(warning_model)){
    if(warning_model!=1){
      CI <- parameterEstimates(sem_result, boot.ci.type = "bca.simple")
      list_output <- data.frame(protein,CI[CI$label=='a','est'],CI[CI$label=='a','pvalue'],CI[CI$label=='b','est'],CI[CI$label=='b','pvalue'],CI[CI$label=='total','est'],CI[CI$label=='total','pvalue'],CI[CI$label=='cp','est'],CI[CI$label=='cp','pvalue'],CI[CI$label=='ab','est'],CI[CI$label=='ab','pvalue'])
      colnames(list_output) <- c('protein','a_coef','a_p','b_coef','b_p','c_coef','c_p','cp_coef','cp_p','ab_coef','ab_p')
      mentalframe <- rbind(mentalframe, list_output)
      print(result)
      print(CI)
    }
    }
}
write.table(cognitionframe,'/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/cog_mediation_SEM_lessStrict_customizedSelection_loose.csv',sep=',',row.names=F)
write.table(mentalframe,'/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/mental_mediation_SEM_lessStrict_customizedSelection_loose.csv',sep=',',row.names=F)

sink()
