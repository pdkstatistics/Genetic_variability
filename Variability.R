## BASIC CLEANING
rm(list=ls(all=TRUE))

# RUN ALL THE  LINES
options(max.print = 99999, # PRINT MAXIMUM ROWS
        scipen = TRUE)     # NOT IN SCIENTIFIC FORM



##################################################
## PACKAGES
# genetic variability
library("agricolae")
library("variability")

# Clustering
library('biotools')
library("factoextra")
library("clv")


##################################################
## DATA IMPORT
# Location: D:\WORKS\OTHERS WORKS\KP21 Sherin GPB

# With Replication
#DFR = 
    read.table("clipboard", 
                 header=TRUE,
                 check.names = FALSE)

head(DFR)

# Without Replication
#DF = 
    read.table("clipboard", 
               header=TRUE, 
               check.names = FALSE)
head(DF)


##################################################
# GENETIC VARIABILITY
DF_GV = DFR

# As factors
DF_GV[,1] = as.factor(DF_GV[,1])
DF_GV[,2] = as.factor(DF_GV[,2])
str(DF_GV) # verify
dim(DF_GV)
# ANOVA and genetic characteristics
cat("\nANOVA and genetic characteristics",
    "\n---------------------------------\n")
gen.var(data=DF_GV[,-(1:2)],
        genotypevector = DF_GV[,2],
        replicationvector = DF_GV[,1])

# Genotypic correlation
cat("\nGenotypic correlation",
    "\n---------------------\n")
geno.corr(data = DF_GV[,-(1:2)],
          genotypes = DF_GV[,2],
          replication = DF_GV[,1]
          )

# Phenotypic correlation
cat("\nPhenotypic correlation",
    "\n----------------------\n")
pheno.corr(data = DF_GV[,-(1:2)],
          genotypes = DF_GV[,2],
          replication = DF_GV[,1]
          )

geno.path(dependent.var = DF_GV$SYP,
          independent.var = DF_GV[,3:10],
          genotypes = DF_GV$G,
          replication = DF_GV$R)

pheno.path(dependent.var = DF_GV$SYP,
          independent.var = DF_GV[,3:10],
          genotypes = DF_GV$G,
          replication = DF_GV$R)

sink("D:/WORKS/OTHERS WORKS/KKM Ramchander GPB/23.05.04 Root data/GV_Results.txt",
     #append = TRUE
     )
sink()

################################################
## CORRELATION
DF_COR = DFR[,-(1:2)]
#DF_COR = DF[,-1]
correlation(DF_COR)

# plot
COR = metan::corr_coef(DF_COR)
plot(COR)

?plot()
################################################
## CLUSTER ANALYSIS
# Scaling
str(DF)
DF2 = DF[,c(1,2,3,6,11,13,17,18)]
str(DF2)
DF1=scale(DF2,
          center = TRUE,
          scale = TRUE)
DF1

# Mahanalobis Distance
MD = D2.dist(DF1, cov(DF1))
M_MD = as.matrix(MD)  # as matrix
M_MD
clipr::write_clip(M_MD)

# Hierarchical Clustering
MD_HC = hclust(d = MD, method = "complete")
?hclust()

    
# Silhouette method of finding clusters
fviz_nbclust(DF1, 
             hcut, 
             nstart = 25, 
             method = "silhouette", 
             nboot = 50,
             print.summary = TRUE
             )+
    labs(subtitle = "Gap statistic method")   


MD_HC_CG = cutree(MD_HC, k=3)  # CG - cluster group
MD_HC_CG


# DENDROGRAM
fviz_dend(MD_HC, 
          k = 6, # Cut in k groups
          #h= 40, #cluster based on distance
          #show_labels =TRUE,
          color_labels_by_k = TRUE, # color labels by groups
          cex = 0.70, # label size
          lwd = 0.75, # line width
          type = "rectangle",
          #horiz = FALSE,
          rect = TRUE, # Add rectangle around groups
          rect_fill=TRUE,
          #xlab = "",
          #ylab = "Mahalanobis distance",
          theme = theme(
              plot.margin = margin(20, 20, 20, 20)
          ),
          )

?fviz_dend()

# Inter and Intra Cluster Distance
cls.scatt.diss.mx(M_MD, MD_HC_CG)



###################################################
# PERCENTAGE CONTRIBUTION
DF_PER = scale(DF,
               center = TRUE,
               scale = TRUE)

biotools::singh(DF_PER,
                cov(DF_PER))


###################################################
## PATH ANALYSIS
variable.names(DF)
Y = DF[,20] # 14:Yld
X = DF[,-c(20)] 
head(Y);head(X)
COR.Y<-correlation(Y,X)$correlation
COR.X<-correlation(X)$correlation
path.analysis(COR.X,COR.Y)


##################################################
sink("D:/WORKS/OTHERS WORKS/KP21 Sathya Varsha GPB//PathResults.txt", 
     append = TRUE)
sink()
