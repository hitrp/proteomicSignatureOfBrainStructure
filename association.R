getsampleIntersection <- function(x,y,c,field,excludePatients=FALSE){
  if(!excludePatients){
    # get the intersection of the three matrix by certain column
    temp <- intersect(x[,field],y[,field])
    intersection_field <- intersect(temp, c[,field])
    rownames(x) <- x[,field]
    rownames(y) <- y[,field]
    rownames(c) <- c[,field]
    x[,field] <- NULL
    intersection_x <- as.matrix(x[intersection_field,])
    y[,field] <- NULL
    intersection_y <- as.matrix(y[intersection_field,])
    c[,field] <- NULL
    intersection_c <- as.matrix(c[intersection_field,])

  }else{
    # get the intersection of the three matrix by certain column
    temp <- intersect(x[,field],y[,field])
    intersection_field <- intersect(temp, c[,field])
    patients <- read.table('/public/home/pren/Project/UKB/outcome/excludingSubjs.csv',sep=',')
    intersection_field <- setdiff(intersection_field, patients[,1])
    rownames(x) <- x[,field]
    rownames(y) <- y[,field]
    rownames(c) <- c[,field]
    x[,field] <- NULL
    intersection_x <- as.matrix(x[intersection_field,])
    y[,field] <- NULL
    intersection_y <- as.matrix(y[intersection_field,])
    c[,field] <- NULL
    intersection_c <- as.matrix(c[intersection_field,])
  }
  return(list(intersection_x,intersection_y,intersection_c))
}


inverseNorm_transformation <- function(data){
  library(RNOmni)
  return(apply(data,2,function(x)RankNorm(x)))
}


# relationship between protein and regional volume
X <- read.table('/public/home/pren/Project/data/imputated_median_scaled_bl_intersection_LRvolume_aparc_original_INT_AgeSexremoved_protemic.csv',sep = ',',header=FALSE)
colnames(X)[colnames(X)=='V1'] <- 'eid'
X_orig <- as.matrix(X[,2:dim(X)[2]])
Y <- read.table('/public/home/pren/Project/data/imputated_median_scaled_bl_intersection_LRvolume_aparc_original.csv',sep = ',',header=TRUE)
Y_orig <- as.matrix(Y[,2:dim(Y)[2]])
# check if the subjs in same order
X_order <-  read.table('/public/home/pren/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_original_subjs.csv',sep = ',',header=TRUE)
Y_order <- read.table('/public/home/pren/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_original_subjs.csv',sep = ',',header=TRUE)
identical(X_order,Y_order)
eid <- X_order[,1]
# load covariates
cov <- read.table('/public/home/pren/Project/UKB/cova_data_more_withPC_all2visits.csv',sep = ',',header=TRUE)
cov <- cov[!is.na(cov$smoking2) & !is.na(cov$drinking2) & !is.na(cov$Ethnicity) & cov$smoking2!=-3 & cov$drinking2!=-3 & cov$Ethnicity!=0,]
cov_orig <- merge(table(eid),cov,by='eid',all.x=TRUE)
# check for nan values
cov_orig <- cov_orig[,c('eid','age2','interval','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI2','smoking2','drinking2','Ethnicity')]
cov_orig <- cov_orig[!apply(cov_orig,1,function(x){any(is.nan(x))}),]
cov_orig <- cov_orig[!apply(cov_orig,1,function(x){any(is.na(x))}),]
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "drinking2", remove_first_dummy = TRUE)
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "Ethnicity", remove_first_dummy = TRUE)
cov_orig <- cov_orig[,c('eid','age2','interval','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI2','smoking2','drinking2_1','drinking2_2','Ethnicity_2','Ethnicity_3','Ethnicity_4','Ethnicity_5','Ethnicity_6')]
regionNames=read.table('/public/home/pren/Project/data/imputated_median_scaled_bl_intersection_LR_aparc_regions.csv',sep = ',',header = TRUE)
regionNames=as.character(regionNames[2:length(regionNames)])
proteinNames=read.table('/public/home/pren/Project/data/imputated_median_scaled_bl_proteinNames.csv', sep = ',', header=TRUE)
proteinNames=as.character(proteinNames[,1])
list_samples <- getsampleIntersection(X,Y,cov_orig,'eid',TRUE)
X <- list_samples[[1]]
Y <- list_samples[[2]]
Y <- inverseNorm_transformation(Y)
C <- list_samples[[3]]
modeltype <- 'gaussian'
modality <- 'volume'
prefix <- 'imputated_median_scaled_bl_intersection_LRvolume_aparc_allremoved'
rownames(cov) <- cov$eid
tempout <- cbind(rownames(C),Y)
colnames(tempout) <- c('eid',regionNames)
write.table(tempout,'/public/home/pren/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/association_finalvolumeData.csv',row.names = F,sep=',')
demography_subjs <- rownames(C)
demography_volume <- cov[demography_subjs,c('eid','age2','interval','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI2','smoking2','drinking2','Ethnicity')]
demography_volume$Group <- rep('volume',nrow(demography_volume))
saveRDS(demography_volume,'/public/home/pren/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/demography_volume.rds')


## model linear regression between multiple proteins and regional volumes
getCategories <- function(inputProteins){
  assay=read.table('/public/home/pren/Project/data/olink_assay.dat', sep='\t', header=TRUE)
  list_assay <- list()
  list_assay$oncology <- assay[grepl('Oncology',assay[,3]),1]
  list_assay$neurology <- assay[grepl('Neurology',assay[,3]),1]
  list_assay$cardioMet <-  assay[grepl('Cardiometabolic',assay[,3]),1]
  list_assay$inflammation <- assay[grepl('Inflammation',assay[,3]),1]
  intersection <- sapply(list_assay,function(x) length(intersect(x, inputProteins)))
  return(intersection)
}
# Define the function to perform linear regression
perform_regression <- function(i, j, modeltype='gaussian',modality) {
  y <- Y[, j]
  x <- X[, i]
  if(modeltype=='gaussian'){
    # BMI2 is not included, it's colinear with eTIV and reduce too much asssociations
    #model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+BMI2+smoking2+drinking2_1+drinking2_2+Ethnicity_1+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    if (modality=='volume'){
      model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }else if (modality=='thick'){
      model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }else if (modality=='area'){
      model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }else if (modality=='FC'){
      model <- glm(y~x+age2+interval+Sex+headmotion2+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }else if (modality=='DTI_FA'){
      model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }else if(modality=='DTI_MD'){
      model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }else if(modality=='wholebrain'){
      model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=gaussian)
    }
  }else if(modeltype=='quasipoisson'){
    if (modality=='SC'){
    model <- glm(y~x+age2+interval+Sex+eTIV+Site1+Site2+Education+Townsend_index+smoking2+drinking2_1+drinking2_2+Ethnicity_2+Ethnicity_3+Ethnicity_4+Ethnicity_5+Ethnicity_6, cbind(data.frame(x,y),C),family=quasipoisson)
    }
  }
  result <- summary(model)
  return(paste(as.character(result$coefficients[2, 1]),as.character(result$coefficients[2,2]),as.character(result$coefficients[2,3]),as.character(result$coefficients[2,4]),sep='_'))
}
library(foreach)
library(doParallel)
# Perform linear regression for each combination of columns using foreach
num_cores <- detectCores()
registerDoParallel(120)  # Adjust the number of cores as needed
results <- foreach(i = 1:ncol(X), .export = c("X", "Y", 'C','modeltype','modality'), .combine = rbind) %:%
  foreach(j = 1:ncol(Y), .export = c("X", "Y",'C', 'modeltype','modality'), .combine = c) %dopar% {
    perform_regression(i, j,modeltype,modality)
  }
# Apply the function to each element in the matrix
results_beta <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][1])))
results_se <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][2])))
results_tvalue <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][3])))
results_pvalue <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][4])))
padj_BH_overall_vector <- p.adjust(as.vector(results_pvalue), method='BH')
padj_BH_overall <- matrix(padj_BH_overall_vector,nrow = nrow(results_pvalue))
padj_BH_oneProtein <- t(apply(results_pvalue, 1, function(x) p.adjust(x, method='BH')))
padj_Bonf_oneProtein <- t(apply(results_pvalue, 1, function(x) p.adjust(x, method='bonferroni')))
padj_BH_oneRegion <- apply(results_pvalue, 2, function(x) p.adjust(x, method='BH'))
padj_Bonf_oneRegion <- apply(results_pvalue, 2, function(x) p.adjust(x, method='bonferroni'))
statistics_result <- function(pmatrix, dim_index, threshold=0.05){
  # get the number of significant relationships for each protein
  sig_bool <- pmatrix < threshold
  temp <- results_beta
  temp[!sig_bool]=0
  if(dim_index==1){
    count_sig <- rowSums(sig_bool)
    count_sig_pos <- rowSums(temp>0)
    count_sig_neg <- rowSums(temp<0)
    # get sig regions for each protein
    result_sig <- apply(sig_bool,1,function(x){regionNames[x]})
    beta_sig <- apply(temp*sig_bool,1,function(x) x[x!=0])
    se_sig <- apply(results_se*sig_bool,1,function(x) x[x!=0])
    pvalue_sig <- apply(results_pvalue*sig_bool,1,function(x) x[x!=0])
    adjpvalue_sig <- apply(pmatrix*sig_bool,1,function(x) x[x!=0])
  }else if(dim_index==2){
    count_sig <- colSums(sig_bool)
    count_sig_pos <- colSums(temp>0)
    count_sig_neg <- colSums(temp<0)
    # get sig proteins for each column
    result_sig <- apply(sig_bool,2,function(x){proteinNames[x]})
    beta_sig <- apply(temp*sig_bool,2,function(x){x[x!=0]})
    se_sig <- apply(results_se*sig_bool,2,function(x) x[x!=0])
    pvalue_sig <- apply(results_pvalue*sig_bool,2,function(x) x[x!=0])
    adjpvalue_sig <- apply(pmatrix*sig_bool,2,function(x) x[x!=0])
  }
  names(count_sig) <- NULL
  return(list(count_sig, count_sig_pos, count_sig_neg, result_sig,beta_sig,se_sig, pvalue_sig, adjpvalue_sig))
}
count_sig_oneProtein_overall_BH <- statistics_result(padj_BH_overall, 1)[[1]]
count_sig_oneProtein_overall_BH_pos <- statistics_result(padj_BH_overall, 1)[[2]]
count_sig_oneProtein_overall_BH_neg <- statistics_result(padj_BH_overall, 1)[[3]]
regions_sig_oneProtein_overall_BH <- statistics_result(padj_BH_overall, 1)[[4]]
regions_sig_oneProtein_overall_BH_beta <- statistics_result(padj_BH_overall, 1)[[5]]
regions_sig_oneProtein_overall_BH_se <- statistics_result(padj_BH_overall, 1)[[6]]
regions_sig_oneProtein_overall_BH_pvalue <- statistics_result(padj_BH_overall, 1)[[7]]
regions_sig_oneProtein_overall_BH_adjpvalue <- statistics_result(padj_BH_overall, 1)[[8]]
sig_index_oneProtein_overall_BH <- which(count_sig_oneProtein_overall_BH>0)
count_sig_oneRegion_overall_BH <- statistics_result(padj_BH_overall, 2)[[1]]
count_sig_oneRegion_overall_BH_pos <- statistics_result(padj_BH_overall, 2)[[2]]
count_sig_oneRegion_overall_BH_neg <- statistics_result(padj_BH_overall, 2)[[3]]
proteins_sig_oneRegion_overall_BH <- statistics_result(padj_BH_overall, 2)[[4]]
proteins_sig_oneRegion_overall_BH_beta <- statistics_result(padj_BH_overall, 2)[[5]]
proteins_sig_oneRegion_overall_BH_se <- statistics_result(padj_BH_overall, 2)[[6]]
proteins_sig_oneRegion_overall_BH_pvalue <- statistics_result(padj_BH_overall, 2)[[7]]
proteins_sig_oneRegion_overall_BH_adjpvalue <- statistics_result(padj_BH_overall, 2)[[8]]
sig_index_oneRegion_overall_BH <- which(count_sig_oneRegion_overall_BH>0)

# show distribution across categories
sig_oneProtein_overall_BH <- data.frame(count=count_sig_oneProtein_overall_BH[sig_index_oneProtein_overall_BH], proteinName=proteinNames[sig_index_oneProtein_overall_BH])
sig_oneProtein_overall_BH_all <- data.frame(count=count_sig_oneProtein_overall_BH, count_pos=count_sig_oneProtein_overall_BH_pos, count_neg=count_sig_oneProtein_overall_BH_neg, proteinName=proteinNames)
sig_oneProtein_overall_BH <- sig_oneProtein_overall_BH[order(sig_oneProtein_overall_BH$count, decreasing = TRUE), ]
sig_oneProtein_overall_BH_category <- getCategories(sig_oneProtein_overall_BH[,'proteinName'])
sig_oneRegion_overall_BH <- data.frame(count=count_sig_oneRegion_overall_BH[sig_index_oneRegion_overall_BH], regionName=regionNames[sig_index_oneRegion_overall_BH])
sig_oneRegion_overall_BH <- sig_oneRegion_overall_BH[order(sig_oneRegion_overall_BH$count, decreasing = TRUE), ]
sig_oneRegion_overall_BH_all <- data.frame(count=count_sig_oneRegion_overall_BH, count_pos=count_sig_oneRegion_overall_BH_pos, count_neg=count_sig_oneRegion_overall_BH_neg,regionName=regionNames)


my_data_list <- list(sig_oneProtein_overall_BH, sig_oneProtein_overall_BH_all, sig_oneRegion_overall_BH,sig_oneRegion_overall_BH_all, regions_sig_oneProtein_overall_BH,regions_sig_oneProtein_overall_BH_beta,regions_sig_oneProtein_overall_BH_se,regions_sig_oneProtein_overall_BH_pvalue,regions_sig_oneProtein_overall_BH_adjpvalue,proteins_sig_oneRegion_overall_BH,proteins_sig_oneRegion_overall_BH_beta,proteins_sig_oneRegion_overall_BH_se,proteins_sig_oneRegion_overall_BH_pvalue,proteins_sig_oneRegion_overall_BH_adjpvalue,regionNames,proteinNames,sig_oneProtein_overall_BH_category, results_beta, results_se, results_tvalue, padj_BH_overall)
names(my_data_list) <- c('sig_oneProtein_overall_BH','sig_oneProtein_overall_BH_all','sig_oneRegion_overall_BH','sig_oneRegion_overall_BH_all','regions_sig_oneProtein_overall_BH','regions_sig_oneProtein_overall_BH_beta','regions_sig_oneProtein_overall_BH_se','regions_sig_oneProtein_overall_BH_pvalue','regions_sig_oneProtein_overall_BH_adjpvalue','proteins_sig_oneRegion_overall_BH','proteins_sig_oneRegion_overall_BH_beta','proteins_sig_oneRegion_overall_BH_se','proteins_sig_oneRegion_overall_BH_pvalue','proteins_sig_oneRegion_overall_BH_adjpvalue','regionNames','proteinNames','sig_oneProtein_overall_BH_category','results_beta','results_se','results_tvalue','padj_BH_overall')
saveRDS(my_data_list, file = paste('/public/home/pren/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/',prefix,"_result.rds",sep=''))
save.image(file = paste('/public/home/pren/Project/PLS_matlab_job_vol_mean_reg_z_UnivariateRegression/',prefix,"_environment.RData",sep=''))
