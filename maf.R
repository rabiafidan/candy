
setwd("/nemo/project/proj-tracerX/working/SRM/RABIA/MAF/multi")
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(version = '3.16')
BiocManager::install("maftools")

install.packages("RColorBrewer")
library('maftools')
library("RColorBrewer")


badBQ=c("UME_UME003_N1d1.maf","UME_UME003_R3d1.maf","UME_UME003_R3d1.maf","UME_UME003_R3d1.maf","UME_UME004_R2d1.maf","UME_UME005_R2d1.maf","UME_UME005_R2d1.maf","UME_UME006_R1d1.maf","UME_UME007_R2d1.maf","UME_UME008_R1d1.maf","UME_UME008_R1d1.maf","UME_UME009_N1d1.maf","UME_UME011_N1d1.maf","UME_UME012_N1d1.maf","UME_UME012_N1d1.maf","UME_UME012_N1d1.maf","UME_UME012_N1d1.maf","UME_UME013_R1d1.maf","UME_UME014_N1d1.maf","UME_UME015_R2d1.maf","UME_UME021_B1d1.maf","UME_UME021_R2d1.maf","UME_UME021_R2d1.maf","UME_UME022_R1d1.maf","UME_UME022_R1d1.maf","UME_UME023_B1d1.maf","UME_UME025_B1d1.maf","UME_UME025_B1d1.maf","UME_UME025_B1d1.maf","UME_UME026_R1d1.maf","UME_UME026_R1d1.maf","UME_UME026_R1d1.maf","UME_UME028_R1d1.maf","UME_UME029_B1d1.maf","UME_UME029_B1d1.maf","UME_UME029_B1d1.maf","UME_UME029_R1d1.maf")
badBQ=unique(badBQ)

maf_files <- list.files("/nemo/project/proj-tracerX/working/SRM/RABIA/MAF_v2/multi")

merged = merge_mafs(mafs = maf_files,
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

oncoplot(maf = merged,draw_titv = TRUE, showTumorSampleBarcodes=T,SampleNamefontSize=0.7, clinicalFeatures = c('tumour_size',"Age","sex"))
oncoplot(maf = merged_good,draw_titv = TRUE, showTumorSampleBarcodes=T,SampleNamefontSize=0.7, clinicalFeatures = c('tumour_size',"Age","sex"), 
         sortByAnnotation=T, annotationOrder=levels(clinical$tumour_size),genes=genes_of_interest,barcode_mar=6)
oncoplot(maf = merged_bad,draw_titv = TRUE, showTumorSampleBarcodes=T,SampleNamefontSize=0.7, clinicalFeatures = c('tumour_size',"Age","sex"))




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

