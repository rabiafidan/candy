

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(version = '3.16')
BiocManager::install("maftools")

install.packages("RColorBrewer")
library('maftools')
library("RColorBrewer")


badBQ=c("UME_UME003_N1d1.maf","UME_UME003_R3d1.maf","UME_UME003_R3d1.maf","UME_UME003_R3d1.maf","UME_UME004_R2d1.maf","UME_UME005_R2d1.maf","UME_UME005_R2d1.maf","UME_UME006_R1d1.maf","UME_UME007_R2d1.maf","UME_UME008_R1d1.maf","UME_UME008_R1d1.maf","UME_UME009_N1d1.maf","UME_UME011_N1d1.maf","UME_UME012_N1d1.maf","UME_UME012_N1d1.maf","UME_UME012_N1d1.maf","UME_UME012_N1d1.maf","UME_UME013_R1d1.maf","UME_UME014_N1d1.maf","UME_UME015_R2d1.maf","UME_UME021_B1d1.maf","UME_UME021_R2d1.maf","UME_UME021_R2d1.maf","UME_UME022_R1d1.maf","UME_UME022_R1d1.maf","UME_UME023_B1d1.maf","UME_UME025_B1d1.maf","UME_UME025_B1d1.maf","UME_UME025_B1d1.maf","UME_UME026_R1d1.maf","UME_UME026_R1d1.maf","UME_UME026_R1d1.maf","UME_UME028_R1d1.maf","UME_UME029_B1d1.maf","UME_UME029_B1d1.maf","UME_UME029_B1d1.maf","UME_UME029_R1d1.maf")
badBQ=unique(badBQ)
badSamp=gsub(".maf","",badBQ)
setwd("/nemo/project/proj-tracerX/working/SRM/RABIA/MAF/multi")
maf_files <- Sys.glob(file.path("/nemo/project/proj-tracerX/working/SRM/RABIA/MAF/multi", "*.maf"))
maf_single=Sys.glob(file.path("/nemo/project/proj-tracerX/working/SRM/RABIA/MAF/single", "*.maf"))

merged = merge_mafs(mafs = c(maf_files,maf_single),
                    vc_nonSyn = c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins",
                                  "Frame_Shift_Del",  "In_Frame_Ins", "In_Frame_Del", "Splice_Site"),
                    clinicalData = clinical)

merged_good=merge_mafs(mafs = maf_files[!maf_files %in% badBQ],
                       vc_nonSyn = c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins",
                                     "Frame_Shift_Del",  "In_Frame_Ins", "In_Frame_Del", "Splice_Site"),
                       clinicalData = clinical)
merged_bad=merge_mafs(mafs = maf_files[maf_files %in% badBQ],
                      vc_nonSyn = c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins",
                                    "Frame_Shift_Del",  "In_Frame_Ins", "In_Frame_Del", "Splice_Site"),
                      clinicalData = clinical)

genes_of_interest = c("VHL", "PBRM1", "SETD2", "PIK3CA", "MTOR", "PTEN", "KDM5C", "CSMD3", "BAP1", "TP53", "TSC1", "TSC2", "ARID1A", "TCEB1")


dataf=merged@data
sum(dataf$Reference_Allele != dataf$Tumor_Seq_Allele1) #0
SNP=dataf[dataf$Variant_Type=="SNP"]
SNP$mut=paste(SNP$Reference_Allele,SNP$Tumor_Seq_Allele2,sep = "-")

CT=SNP[SNP$mut%in%c("C-T","G-A"),]$t_VAF
CG=SNP[SNP$mut%in%c("C-G","G-C"),]$t_VAF
CA=SNP[SNP$mut%in%c("C-A","G-T"),]$t_VAF
TA=SNP[SNP$mut%in%c("T-A","A-T"),]$t_VAF
TG=SNP[SNP$mut%in%c("T-G","A-C"),]$t_VAF
TC=SNP[SNP$mut%in%c("T-C","A-G"),]$t_VAF


par(mfcol=c(2,3))
par(mar=c(2,2,2,2))
hist(CT,breaks = 40,col = "darkred",xlim=c(0,1),xlab='',ylab="",main = "C>T")
hist(CG,breaks = 40,col = "darkred",xlim=c(0,1),xlab='',ylab="",main = "C>G")
hist(CA,breaks = 40,col = "darkred",xlim=c(0,1),xlab='',ylab="",main = "C>A")
hist(TA,breaks = 40,col = "darkred",xlim=c(0,1),xlab='',ylab="",main = "T>A")
hist(TG,breaks = 40,col = "darkred",xlim=c(0,1),xlab='',ylab="",main = "T>G")
hist(TC,breaks = 40,col = "darkred",xlim=c(0,1),xlab='',ylab="",main = "T>C")

### Individual samples ######

ind="UME_UME007_R2d1.maf"
nrow(SNP[SNP$Source_MAF==ind,])

nrow(dataf[dataf$Source_MAF==ind])
CT=SNP[SNP$mut%in%c("C-T","G-A")&SNP$Source_MAF==ind,]$t_VAF
CG=SNP[SNP$mut%in%c("C-G","G-C")&SNP$Source_MAF==ind,]$t_VAF
CA=SNP[SNP$mut%in%c("C-A","G-T")&SNP$Source_MAF==ind,]$t_VAF
TA=SNP[SNP$mut%in%c("T-A","A-T")&SNP$Source_MAF==ind,]$t_VAF
TG=SNP[SNP$mut%in%c("T-G","A-C")&SNP$Source_MAF==ind,]$t_VAF
TC=SNP[SNP$mut%in%c("T-C","A-G")&SNP$Source_MAF==ind,]$t_VAF


par(mfcol=c(2,3))
par(mar=c(4,2,2,2))
hist(CT,breaks = 40,col = "darkred",xlim=c(0,1),main=ind,ylab="",xlab = "C>T")
hist(CG,breaks = 40,col = "darkred",xlim=c(0,1),main='',ylab="",xlab = "C>G")
hist(CA,breaks = 40,col = "darkred",xlim=c(0,1),main='',ylab="",xlab = "C>A")
hist(TA,breaks = 40,col = "darkred",xlim=c(0,1),main='',ylab="",xlab = "T>A")
hist(TG,breaks = 40,col = "darkred",xlim=c(0,1),main='',ylab="",xlab = "T>G")
hist(TC,breaks = 40,col = "darkred",xlim=c(0,1),main='',ylab="",xlab = "T>C")


tab=read.table('/nemo/project/proj-tracerX/working/SRM/RABIA/SAREK_1/samples_1.tsv',sep='\t')
colnames(tab)=c('patN','sex','sampleT','sampleN','lane','forward','reverse')
tab1=tab[-c(grep('TOR',tab$patN)),]
tab2=read.csv('/nemo/project/proj-tracerX/working/SRM/RABIA/SAREK_1/20221111_patient_data.csv',header = T)
tab2$crick_ID=gsub('UME_','',tab2$crick_ID)
grep('21',tab2$crick_ID)
tab2[66,'crick_ID']
tab2[66,'crick_ID']='UME021'
tab3=merge(tab1[tab1$sampleT==1,c('patN','sex','sampleT','sampleN')],tab2,by.x = 'patN',by.y='crick_ID')
tab3=unique(tab3)
head(tab3)
#87 new
#8 TOR

clinical=tab3[c(4,1:3,5:ncol(tab3))]
colnames(clinical)=c("Tumor_Sample_Barcode",colnames(clinical)[-1])
clinical$tumour_size=factor(clinical$tumour_size,levels = sort(unique(clinical$tumour_size)))

sizeColors=c(brewer.pal(n = 9, name ="Blues"),brewer.pal(n = 9, name ="Greens"),brewer.pal(n =5, name ="YlOrRd"))
names(sizeColors)=levels(clinical$tumour_size)

oncoplot(maf = merged,draw_titv = F, showTumorSampleBarcodes=T,SampleNamefontSize=0.7, clinicalFeatures = c('tumour_size',"Age","sex"),
         genes=genes_of_interest,barcode_mar=6)
oncoplot(maf = merged_good,draw_titv = F, showTumorSampleBarcodes=T,SampleNamefontSize=0.7, clinicalFeatures = c('tumour_size',"Age","sex"), 
         genes=genes_of_interest,barcode_mar=6,, annotationOrder=levels(clinical$tumour_size), sortByAnnotation=T) #, annotationOrder=levels(clinical$tumour_size), sortByAnnotation=T
oncoplot(maf = merged_bad,draw_titv = TRUE, showTumorSampleBarcodes=T,SampleNamefontSize=0.7, clinicalFeatures = c('tumour_size',"Age","sex"),
         genes=genes_of_interest,barcode_mar=6)

class(merged_good)
merged_good
getFields(merged_good)
plotVaf(maf = merged_good, vafCol = 't_VAF', flip = TRUE)

#Examples  
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
laml <- read.maf(maf = laml.maf, clinicalData = laml.clin)
#Basic onocplot
oncoplot(maf = laml, top = 3)
#Changing colors for variant classifications (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
#Color coding for FAB classification; try getAnnotations(x = laml) to see available annotations.
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)
oncoplot(maf = laml, colors = col, clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, annotationColor = fabcolors)







summary(tab)
pats=unique(tab1$patN)
sexes=tapply(tab1$sex, tab1$patN,unique )
ages=tapply(tab3$Age, tab3$patN,unique )
stage=tapply(tab3$pT_stage, tab3$patN, unique)
tsize=tapply(tab3$tumour_size, tab3$patN, unique)
tums=tapply(tab1$sampleN, tab1$patN, function (x) {length(unique(x))})
tums=tums-1
stg=data.frame(table(stage))
sxs=data.frame(table(sexes))
tns=data.frame(table(tums))

barplot(Freq~stage,data = stg,col='indianred4')
barplot(Freq~tums,data = tns,col='indianred4')
barplot(Freq~sexes,data=sxs,col='indianred4')
hist(ages,col = 'indianred4',breaks = 20,xlab = 'patient age',main = '')
abline(v=median(ages),lty=2,col='blue')
hist(tsize,col = 'indianred4',breaks = 20,xlab = 'tumour size',main = '')
abline(v=median(tsize),lty=2,col='blue')
hist(tums)


maf_tab=read.table("/nemo/project/proj-tracerX/working/SRM/RABIA/MAF/single/RF_RFH010_R1d1.maf",header = T,sep="\t")



################ Mutational signature  ################################


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MutationalPatterns")

library(MutationalPatterns)
library(BSgenome)
available.genomes()


ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(ref_genome, character.only = TRUE)

vcffiles=read.table("/nemo/project/proj-tracerX/working/SRM/RABIA/VEP/name_file_df.tsv", sep = "\t")

grl <- read_vcfs_as_granges(vcffiles$V2, vcffiles$V1, ref_genome)
muts <- mutations_from_vcf(grl[[1]])
head(muts, 12)

#conventional use
types <- mut_type(grl[[1]])
head(types, 12)

#context
context <- mut_context(grl[[1]], ref_genome)
head(context, 12)

#
type_occurrences <- mut_type_occurrences(grl, ref_genome)
type_occurrences
totals=type_occurrences[,1]+type_occurrences[,2]+type_occurrences[,3]+type_occurrences[,4]+
  type_occurrences[,5]+type_occurrences[,6]+type_occurrences[,7]+type_occurrences[,8]
p2 <- plot_spectrum(type_occurrences[totals>15&totals<65,], CT = TRUE, 
                    indv_points = F, by=vcffiles$V1[totals>15&totals<65],error_bars = 'none')
p2
p3 <- plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE,error_bars = 'none')
p3

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_sums=colSums(mut_mat)
head(mut_mat)
plot_96_profile(mut_mat[,mut_sums>10 & mut_sums<65])
pooled_mut_mat <- pool_mut_mat(mut_mat[,mut_sums<65], grouping = rep("All Samples",188))
plot_96_profile(pooled_mut_mat)

                               