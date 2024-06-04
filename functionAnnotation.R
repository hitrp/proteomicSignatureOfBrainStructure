## functional annotation:
# GO enrichment : Webgestalt online
# Tissue enrichment: FUMA GENE2FUN module
# physical distribution on chromosomes, in R language
getPosition <- function(proteins){
  library(topr)
  library(D3GB)
  outframe <- data.frame(chrom=c(),gene_start=c(),gene_end=c(),gene_symbol=c(),biotype=c())
  for(p in proteins){
   outframe <- rbind(outframe, get_gene(p, chr = NULL, build = 37))
  }
  return(outframe)
}
plotframe <- data.frame(chrom=c(),freq=c(),metric=c())
plotframe_density <- data.frame(chrom=c(),freq=c(),metric=c())
for(i in seq(1,length(list_proteins))){
  proteins <- list_proteins[[i]]
  posframe = getPosition(proteins)
  temp = table(posframe$chrom)
  density = as.data.frame(temp)
  colnames(density) <- c('chrom','count')
  chrom_length <- GRCh37
  chrom_length$length  <- chrom_length$end-chrom_length$start
  colnames(chrom_length) <- c('chrom','start','end','length')
  density <- merge(density, chrom_length[,c('chrom','length')], by='chrom', all.x=T)
  density$count <- density$count/density$length
  density <- density[order(density$count),]
  density$metric <- rep(metrics[i],nrow(density))
  freq = as.vector(unname(temp))
  posframe = data.frame(chrom=names(temp),count=freq)
  posframe <- posframe[order(posframe$count),]
  posframe$metric <- rep(metrics[i],nrow(posframe))
  plotframe <- rbind(plotframe, posframe)
  plotframe_density <- rbind(plotframe_density, density)
}

# cell type enrichment, in python language
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from scipy.stats import fisher_exact
# https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
transmitters=pd.read_csv('/data/CellTypes.gmt',sep='^',header=None)
transmitters=transmitters[0].str.split('\t', expand=True)
Astro=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='Astro',:].iloc[0]) if isinstance(s,str)]
Astro=Astro[2:len(Astro)]
Astro=[s for s in Astro if s]
Endo=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='Endo',:].iloc[0]) if isinstance(s,str)]
Endo=Endo[2:len(Endo)]
Endo=[s for s in Endo if s]
Micro=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='Micro',:].iloc[0]) if isinstance(s,str)]
Micro=Micro[2:len(Micro)]
Micro=[s for s in Micro if s]
Neuro_Ex=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='Neuro_Ex',:].iloc[0]) if isinstance(s,str)]
Neuro_Ex=Neuro_Ex[2:len(Neuro_Ex)]
Neuro_Ex=[s for s in Neuro_Ex if s]
Neuro_In=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='Neuro_In',:].iloc[0]) if isinstance(s,str)]
Neuro_In=Neuro_In[2:len(Neuro_In)]
Neuro_In=[s for s in Neuro_In if s]
Oligo=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='Oligo',:].iloc[0]) if isinstance(s,str)]
Oligo=Oligo[2:len(Oligo)]
Oligo=[s for s in Oligo if s]
OPC=[s for s in list(transmitters.loc[transmitters.iloc[:,0]=='OPC',:].iloc[0]) if isinstance(s,str)]
OPC=OPC[2:len(OPC)]
OPC=[s for s in OPC if s]

background = pd.read_csv('proteinNames_2920.txt',header=None)
background = background.iloc[:,0]
background = background.to_list()

def Cellenrichment(assoProteins, metric):
    M=len(background)
    N=len(set(Neuro_Ex) & set(background))
    n=len(assoProteins)
    x=len(set(Neuro_Ex) & set(assoProteins))
    hypergeom.sf(x-1, M, N, n)
    oddsr_Neuro_Ex, p_Neuro_Ex = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')
    
    M=len(background)
    N=len(set(Neuro_In) & set(background))
    n=len(assoProteins)
    x=len(set(Neuro_In) & set(assoProteins))
    hypergeom.sf(x-1, M, N, n)
    oddsr_Neuro_In, p_Neuro_In = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')

    M=len(background)
    N=len(set(Astro) & set(background))
    n=len(assoProteins)
    x=len(set(Astro) & set(assoProteins))
    hypergeom.sf(x-1, M, N, n)
    oddsr_Astro, p_Astro = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')

    M=len(background)
    N=len(set(Micro) & set(background))
    n=len(assoProteins)
    x=len(set(Micro) & set(assoProteins))
    hypergeom.sf(x-1, M, N, n)
    oddsr_Micro, p_Micro = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')

    M=len(background)
    N=len(set(Endo) & set(background))
    n=len(assoProteins)
    x=len(set(Endo) & set(assoProteins))
    hypergeom.sf(x - 1, M, N, n)
    hypergeom.sf(x-1, M, N, n)
    oddsr_Endo, p_Endo = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')

    M=len(background)
    N=len(set(Oligo) & set(background))
    n=len(assoProteins)
    x=len(set(Oligo) & set(assoProteins))
    hypergeom.sf(x-1, M, N, n)
    oddsr_Oligo, p_Oligo = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')

    M=len(background)
    N=len(set(OPC) & set(background))
    n=len(assoProteins)
    x=len(set(OPC) & set(assoProteins))
    hypergeom.sf(x-1, M, N, n)
    oddsr_OPC, p_OPC = fisher_exact(np.array([[x,n-x],[N-x,M-(n+N)+x]]), alternative='greater')

    pvalues=[p_Neuro_Ex,p_Neuro_In,p_Endo,p_Micro,p_Astro,p_Oligo,p_OPC]
    import statsmodels.api as sm
    fdr_corrected_p_values = sm.stats.multipletests(pvalues, method='fdr_bh')[1]

    EnrichRatio=[oddsr_Neuro_Ex,oddsr_Neuro_In,oddsr_Endo,oddsr_Micro,oddsr_Astro,oddsr_Oligo,oddsr_OPC]
    label=['Neuro_Ex','Neuro_In','Endo','Micro','Astro','Oligo','OPC']
    metrics=[metric]*7
    return pd.DataFrame({'pvalue':pvalues,'FDR':fdr_corrected_p_values, 'EnrichRatio':EnrichRatio, 'label':label,
    'metric':metrics})

files = ['volume_proteins.txt','thick_proteins.txt','area_proteins.txt','DTIFA_proteins.txt','DTIMD_proteins.txt']
outputframe = pd.DataFrame({'pvalue':[],'FDR':[],'EnrichRatio':[],'label':[],'metric':[]})
for s in files:
    assoProteins = pd.read_csv('H:/'+s,header=None)
    assoProteins = assoProteins.iloc[:,0]
    assoProteins = assoProteins.to_list()
    metric = s.replace('temp_','')
    metric = metric.replace('.txt','')
    outputframe = pd.concat([outputframe, Cellenrichment(assoProteins, metric)])

outputframe.to_csv('EachModality_enrichment_Celltype.csv')
