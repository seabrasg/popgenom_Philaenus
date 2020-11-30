# Code for performing and plotting side-byside Principal components scatterplots
# Principal Component Analysis done on aedeagus measurements and done on SNP RADseq genotypes
##############################################################################

#Setting working directory
setwd("/Users/seabrasg/Dropbox/0_Philaenus_RAD/Paper_RAD/Submission_PeerJ/docs_submitted/raw_data_Fig3")


#Installing required packages
#Remove "#" from the next code lines, if packages are not installed
#install.packages("ggplot2") 
#install. packages("gridExtra")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")
#Loading required packages
require(ggplot2)
require(gridExtra)

require(gdsfmt)
require("SNPRelate")


#### PCA of aedeagus data
###################################################
#Importing data
dados_aedeagus <- read.delim("./Aedeagus_data.txt")

#Principal Component Analysis
pca_aedeagus <- prcomp(~ TotLen+UpMeanLen+MidMeanLen+LowMeanLen+UpMeanCur,
                data = dados_aedeagus, center = TRUE, scale = TRUE, na.action = na.omit) 
  
#Creating a dataframe with PCA scores
PCA_aedeagus_scores <- data.frame(na.omit(dados_aedeagus)[,1:4], pca_aedeagus$x)
  
#Creating a dataframe PCA loadings
PCA_aedeagus_loadings <- data.frame(Variables = rownames(pca_aedeagus$rotation), pca_aedeagus$rotation)
write.csv(PCA_aedeagus_loadings, "./PCA_aedeagus_loadings.csv")

#Getting standard deviation of PCs 
PCA_aedeagus_sdev <- pca_aedeagus$sdev 

#Calculating proportion of explained variance by PC1, PC2 and PC3
prop_var_PC1_aedeagus <- (PCA_aedeagus_sdev^2)[1]/sum((PCA_aedeagus_sdev^2)[1:5])*100
prop_var_PC2_aedeagus <- (PCA_aedeagus_sdev^2)[2]/sum((PCA_aedeagus_sdev^2)[1:5])*100
prop_var_PC3_aedeagus <- (PCA_aedeagus_sdev^2)[3]/sum((PCA_aedeagus_sdev^2)[1:5])*100

#### PCA of RADseq data
###################################################
#Importing data
vcf.fn<-"./RADseq_data.vcf"#load the data
snpgdsVCF2GDS(vcf.fn,"test.gds",method="biallelic.only")#import to class. use method= biallelic.only for biallelic snps

#summary of the vcf imported
snpgdsSummary("test.gds")

#open the translated vcf imported
genofile<-snpgdsOpen("test.gds")

#import population names (file with a list of names of the population to which each individual belong - one per line and same order as in vcf file) )
pop_code <- scan("./popmap.csv", what=character())

#PCA
pca_RADseq<-snpgdsPCA(genofile)

#get sample ids
sample.id<-read.gdsn(index.gdsn(genofile,"sample.id"))

#data frame of eigenvectors
PCA_RADseq_scores<-data.frame(sample.id= pca_RADseq$sample.id,
                Population=factor((pop_code)[match(pca_RADseq$sample.id, sample.id)], levels=c("MOR","USA","UK","AZO","POR","FRAN","GRE","TURK","FIN")),
                PC1= pca_RADseq$eigenvect[,1],# first eigenvector 
                PC2= pca_RADseq$eigenvect[,2], 
                PC3= pca_RADseq$eigenvect[,3],
                stringsAsFactors=FALSE)

# variance proportion (%)
pc.percent<-pca_RADseq$varprop*100

snpgdsClose(genofile)

## MULTIPLE PLOT OF PCAs
#############################################################################
svg(file = "./plot_PCA_aedeagus_RADseq.svg", bg="white",width = 6, height = 5)

# Colours to use
c_grey<- "#bfbfbf" #MOR
c_pink<- "#fb9a99" # USA
c_darkblue<- "#0c5eba" # UK
c_purple<- "#4B0082" # Azores
c_green<- "#0cd497" #Portugal
c_olivegreen<- "#8fc637" # FRAN
c_red<- "#f24e4e" #Greece
c_yellow<- "#fdd200" # Turkey
c_blue<- "#10aae2" #Finland

# Shapes to use
filled_circle<-19 #MOR
times<-4 # USA
plus<-3 # UK
filled_losango<-18 # Azores
empty_circle<-1 #Portugal
star<-8 # FRAN
empty_triangle<-2 #Greece
filled_triangle<-17 # Turkey
filled_square<-15 #Finland


# AEDEAGUS ANALYSIS - population levels
# Levels: Finland, Greece, Morocco, Portugal, Turkey, UK, USA

colours_aedeagus<-c(c_blue, c_red,c_grey,c_green,c_yellow,c_darkblue,c_pink)
shapes_aedeagus<-c(filled_square,filled_triangle,filled_circle,empty_circle,filled_triangle,plus,times)

# RADseq ANALYSIS - population levels
# Levels: MOR USA UK AZO POR FRAN GRE TURK FIN
colours_RADseq<-c(c_grey, c_pink, c_darkblue,c_purple, c_green, c_olivegreen, c_red, c_yellow, c_blue)
shapes_RADseq<-c(filled_circle,times,plus,filled_losango,empty_circle,star,filled_triangle,filled_triangle,filled_square)


# PLOT A - AEDEAGUS PC1 vs PC2 ######################################################################
PCAplot1<- ggplot(data = PCA_aedeagus_scores,
           aes(x = PCA_aedeagus_scores$PC1, y = PCA_aedeagus_scores$PC2))+
    geom_hline(yintercept=0, linetype="dashed", color = "gray", size=0.3)+
    geom_vline(xintercept=0, linetype="dashed", color = "gray", size=0.3)+
    #change parameter for colour and shape in "geom_point"
    #if symbol shape is contoured, "fill" argument changes colour on symbol filling colour and "colour" changes the symbol contour colour
    geom_point(aes(colour = Population,shape=Population), size = 3)+
    #change the colours in "scale_colour_manual"
    scale_colour_manual(values = colours_aedeagus)+ #specify the colours to use
    #change symbols in "scale_shape_manual"
    
    scale_shape_manual(values = shapes_aedeagus)+ #specify the shapes to use
    #geom_text(aes(label = ID), hjust = 0.5, vjust = 0.3, size = 1.5, alpha = 0.7)+
     geom_segment(data = PCA_aedeagus_loadings, x = 0, y = 0,
                 xend = PCA_aedeagus_loadings$PC1*4,
                 yend = PCA_aedeagus_loadings$PC2*4,
                 arrow = arrow(length = unit(0.25,"cm")), size = 0.3)+
  
    annotate("text", label = PCA_aedeagus_loadings[,1], size = 2, alpha = 0.7,
             x = PCA_aedeagus_loadings$PC1*0.45*sqrt(nrow(PCA_aedeagus_scores)-1),
             y = PCA_aedeagus_loadings$PC2*0.45*sqrt(nrow(PCA_aedeagus_scores)-1))+
    xlab(paste0("PC1", " (", round(prop_var_PC1_aedeagus, 2), "%)"))+
    ylab(paste0("PC2", " (", round(prop_var_PC2_aedeagus, 2), "%)"))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.3, fill = NA),
          plot.title = element_text(hjust = 0.5),
          text = element_text(colour = "black", size = 8),
          axis.title = element_text(colour = "black", size = 8),
          axis.text = element_text(colour = "black", size = 6),
          axis.ticks = element_line(size = 0.3),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm")) +
    annotate("text", x = -3.2, y = 3, label = "A") 
  
# PLOT B - AEDEAGUS PC1 vs PC3 ######################################################################

PCAplot2<- ggplot(data = PCA_aedeagus_scores,
       aes(x = PCA_aedeagus_scores$PC1, y = PCA_aedeagus_scores$PC3))+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", size=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=0.3)+
  geom_point(aes(colour = Population,shape=Population), size = 3)+
  scale_colour_manual(values = colours_aedeagus)+ #specify the colours to use

  scale_shape_manual(values = shapes_aedeagus)+ #specify the shapes to use
  #geom_text(aes(label = ID), hjust = 0.5, vjust = 0.3, size = 1.5, alpha = 0.7)+
  geom_segment(data = PCA_aedeagus_loadings, x = 0, y = 0,
               xend = PCA_aedeagus_loadings$PC1*3,
               yend = PCA_aedeagus_loadings$PC3*3,
               arrow = arrow(length = unit(0.20,"cm")), size = 0.3)+

  annotate("text", label = PCA_aedeagus_loadings[,1], size = 2, alpha = 0.7,
           x = PCA_aedeagus_loadings$PC1*0.45*sqrt(nrow(PCA_aedeagus_scores)-1),
           y = PCA_aedeagus_loadings$PC3*0.45*sqrt(nrow(PCA_aedeagus_scores)-1))+
  xlab(paste0("PC1", " (", round(prop_var_PC1_aedeagus, 2), "%)"))+
  ylab(paste0("PC3", " (", round(prop_var_PC3_aedeagus, 2), "%)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.3, fill = NA),
        plot.title = element_text(hjust = 0.5),
        text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 8),
        axis.text = element_text(colour = "black", size = 6),
        axis.ticks = element_line(size = 0.3),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm")) +
  annotate("text", x = -3.2, y = 2.1, label = "B")


#### RADseq 
# Levels: MOR USA UK AZO POR FRAN GRE TURK FIN
# RADseq PC1 vs PC2

# PLOT C - RADseq PC1 vs PC2 ######################################################################

PCAplot3<- ggplot(data = PCA_RADseq_scores,
                  aes(x = PCA_RADseq_scores$PC1, y = PCA_RADseq_scores$PC2))+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", size=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=0.3)+
  geom_point(aes(colour = Population,shape=Population), size = 3)+
  scale_colour_manual(values = colours_RADseq)+ #specify the colours to use
  scale_shape_manual(values = shapes_RADseq)+
  xlab(paste0("PC1", " (", round(pc.percent[1], 2), "%)"))+
  ylab(paste0("PC2", " (", round(pc.percent[2], 2), "%)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.3, fill = NA),
        plot.title = element_text(hjust = 0.5),
        text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 8),
        axis.text = element_text(colour = "black", size = 6),
        axis.ticks = element_line(size = 0.3),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm")) +
  annotate("text", x = -0.45, y = 0.25, label = "C")

# PLOT D - RADseq PC1 vs PC3 ######################################################################

PCAplot4<- ggplot(data = PCA_RADseq_scores,
                 aes(x = PCA_RADseq_scores$PC1, y = PCA_RADseq_scores$PC3))+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", size=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=0.3)+
  geom_point(aes(colour = Population,shape=Population), size = 3)+
  scale_colour_manual(values = colours_RADseq)+ #specify the colours to use
  scale_shape_manual(values = shapes_RADseq)+
  xlab(paste0("PC1", " (", round(pc.percent[1], 2), "%)"))+
  ylab(paste0("PC3", " (", round(pc.percent[3], 2), "%)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.3, fill = NA),
        plot.title = element_text(hjust = 0.5),
        text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 8),
        axis.text = element_text(colour = "black", size = 6),
        axis.ticks = element_line(size = 0.3),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"))  +
  annotate("text", x = -0.43, y = 0.24, label = "D")

# COMMON LEGEND ######################################################################

#Creating a function to extract the legend to be used as common legend of the four plots
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

g_legend<-function(plot_to_use){
  tmp <- ggplot_gtable(ggplot_build(plot_to_use))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Extracting the common legend
mylegend <- g_legend(PCAplot3)

#Arrange plots side by side
PCAplots <- grid.arrange(
  arrangeGrob(PCAplot1 + theme(legend.position="none"),
              PCAplot2 + theme(legend.position="none"),
              PCAplot3 + theme(legend.position="none"),
              PCAplot4 + theme(legend.position="none"),
              mylegend,
              nrow=2,widths = c(2,2,1), layout_matrix = rbind(c(1,2,NA), c(3,4,5))))

dev.off()
