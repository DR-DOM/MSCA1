#16S MSCA I correction after first exploration (remove: IX_B4; change: VIII_I7 vs VII_B4, VIII_B4 vs VII_I6 )

#################################################################################
#Load libraries and custom functions

#Load libraries
library(vegan)
library(ggplot2)
library(pairwiseAdonis)

#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#Load functions
source(file="cull.otu.script.R")
source(file="generate.long.format.script.R")
source(file="generate.square.dist.script.07.27.2020.R")
source(file="zscore.calculation.script.R")

#################################################################################
#Load data and metadata

#load unifrac data
unifrac=read.csv(file="otu_repr_100.tre1.wsummary_v1.csv")

#load relabund
relabund=read.csv(file="abundance_table_100.database.csv")

#load metadata
metadata=read.csv(file="metadata_corrected.csv")

#################################################################################
#Work up metadata

#change SampleType, Experiment,Treatment, Species, Organism, and Time factor levels
metadata$Set=factor(metadata$Set,levels=c("MSCA", " "))
metadata$SampleType=factor(metadata$SampleType,levels=c("I", "B"))
metadata$Experiment=factor(metadata$Experiment,levels=c("1", "2"))
metadata$Treatment=factor(metadata$Treatment, levels=c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"))
metadata$Species=factor(metadata$Species, levels=c("SW", "Dictyota", "Lobophora", "Acropora", "Diploria", "SWday", "SWnight", "TurfDay", "TurfNight"))
metadata$Organism=factor(metadata$Organism, levels=c("SW", "Algae", "Coral", "Turf"))
metadata$Time=factor(metadata$Time, levels=c("Day", "Night"))

#Subset MSCA data
metadata1=na.omit(metadata,)

#Subset Inoculum
metadata.inoculum=metadata1[metadata1$SampleType=="I",]

#Subset Bioassay
metadata.bioassay=metadata1[metadata1$SampleType=="B",]

#################################################################################
#work up unifrac

#generate squre dist MSCA
generate.square.dist.matrix(unifrac,metadata1,1) 
unifrac.matrix.msca=dist.matrix1 #rename output
unifrac.dist.msca=square.dist.matrix #rename output

#generate squre dist Inoculum
generate.square.dist.matrix(unifrac,metadata.inoculum,1)
unifrac.matrix.inoculum=dist.matrix1 #rename output
unifrac.dist.inoculum=square.dist.matrix #rename output

#generate squre dist Bioassay
generate.square.dist.matrix(unifrac,metadata.bioassay,1)
unifrac.matrix.bioassay=dist.matrix1 #rename output
unifrac.dist.bioassay=square.dist.matrix #rename output



#################################################################################
#generate and work up nmds

#nmds MSCA
nmds.msca=metaMDS(unifrac.dist.msca,k=2,trymax=100)
nmds.msca.scores=as.data.frame(scores(nmds.msca)) #extra scores
nmds.msca.scores$SampleName=rownames(nmds.msca.scores) #add sample column
nmds.msca.scores1=merge(nmds.msca.scores,metadata1,by.x="SampleName",by.y="SampleName")  #merge with metadata
        #Solution reached (20 Runs)

#nmds Inoculum
nmds.inoculum=metaMDS(unifrac.dist.inoculum,k=2,trymax=100)
nmds.inoculum.scores=as.data.frame(scores(nmds.inoculum)) #extra scores
nmds.inoculum.scores$SampleName=rownames(nmds.inoculum.scores) #add sample column
nmds.inoculum.scores1=merge(nmds.inoculum.scores,metadata.inoculum,by.x="SampleName",by.y="SampleName")  #merge with metadata
        #Warning: stress (nearly zero: you may have insufficient data)

#nmds Bioassay
nmds.bioassay=metaMDS(unifrac.dist.bioassay,k=2,trymax=100)
nmds.bioassay.scores=as.data.frame(scores(nmds.bioassay)) #extra scores
nmds.bioassay.scores$SampleName=rownames(nmds.bioassay.scores) #add sample column
nmds.bioassay.scores1=merge(nmds.bioassay.scores,metadata.bioassay,by.x="SampleName",by.y="SampleName")  #merge with metadata
        #Solution reached (20 runs)


#################################################################################
#Plot nmds in ggplot2.

#plot nmds MSCA
ggplot(nmds.msca.scores1,aes(x=NMDS1,y=NMDS2,color=Species,shape=SampleType))+
  geom_point(size=6)+
  geom_text(aes(label=SampleName),hjust=1.5,vjust=0)
        #Inoculums (on right) separate from bioassays (on left)
        #Bioassays clearly cluster together per treatment
        #Inoculum SWnight (VII_B4) seems to be off

#plot nmds Inoculum
ggplot(nmds.inoculum.scores1,aes(x=NMDS1,y=NMDS2,color=Species,shape=Experiment))+
  geom_point(size=2)+
  geom_text(aes(label=Species),hjust=1.5,vjust=0)
        #NMDS warning for low stress: Most inoculums cluster together, TurfNight (IX) and SWnight are off

#plot nmds Bioassay
ggplot(nmds.bioassay.scores1,aes(x=NMDS1,y=NMDS2,color=Organism,shape=Experiment))+
  geom_point(size=5)+
  geom_text(aes(label=Species),hjust=1.5,vjust=0)
        #Experiment 1 more spread out than Experiment 2
        #No clear separation coral vs algae
        #Turfs cluster together
        #SW cluster together

#################################################################################
#################################################################################
#Permanovas

#subset metadata per treatment

# Subset Treatment I -> inoculum got lost during sample processing
metadata.I=metadata1[metadata1$Treatment=="I",]

# Subset Treatment II
metadata.II=metadata1[metadata1$Treatment=="II",]

# Subset Treatment III
metadata.III=metadata1[metadata1$Treatment=="III",]

# Subset Treatment IV
metadata.IV=metadata1[metadata1$Treatment=="IV",]

# Subset Treatment V
metadata.V=metadata1[metadata1$Treatment=="V",]

# Subset Treatment VI
metadata.VI=metadata1[metadata1$Treatment=="VI",]

# Subset Treatment VII
metadata.VII=metadata1[metadata1$Treatment=="VII",]

# Subset Treatment VIII
metadata.VIII=metadata1[metadata1$Treatment=="VIII",]

# Subset Treatment IX
metadata.IX=metadata1[metadata1$Treatment=="IX",]

#################################################################################
#subset unifrac data per treatment

#generate unifrac dist I
generate.square.dist.matrix(unifrac,metadata.I,1)
unifrac.matrix.I=dist.matrix1 #rename output
unifrac.dist.I=square.dist.matrix #rename output

#generate unifrac dist II
generate.square.dist.matrix(unifrac,metadata.II,1)
unifrac.matrix.II=dist.matrix1 #rename output
unifrac.dist.II=square.dist.matrix #rename output

#generate unifrac dist III
generate.square.dist.matrix(unifrac,metadata.III,1)
unifrac.matrix.III=dist.matrix1 #rename output
unifrac.dist.III=square.dist.matrix #rename output

#generate unifrac dist IV
generate.square.dist.matrix(unifrac,metadata.IV,1)
unifrac.matrix.IV=dist.matrix1 #rename output
unifrac.dist.IV=square.dist.matrix #rename output

#generate unifrac dist V
generate.square.dist.matrix(unifrac,metadata.V,1)
unifrac.matrix.V=dist.matrix1 #rename output
unifrac.dist.V=square.dist.matrix #rename output

#generate unifrac dist VI
generate.square.dist.matrix(unifrac,metadata.VI,1)
unifrac.matrix.VI=dist.matrix1 #rename output
unifrac.dist.VI=square.dist.matrix #rename output

#generate unifrac dist VII
generate.square.dist.matrix(unifrac,metadata.VII,1)
unifrac.matrix.VII=dist.matrix1 #rename output
unifrac.dist.VII=square.dist.matrix #rename output

#generate unifrac dist VIII
generate.square.dist.matrix(unifrac,metadata.VIII,1)
unifrac.matrix.VIII=dist.matrix1 #rename output
unifrac.dist.VIII=square.dist.matrix #rename output

#generate unifrac dist IX
generate.square.dist.matrix(unifrac,metadata.IX,1)
unifrac.matrix.IX=dist.matrix1 #rename output
unifrac.dist.IX=square.dist.matrix #rename output

#################################################################################
#Permanovas per treatment

#Permanova I
#adonis2(unifrac.dist.I~SampleType, permutations=999, data=metadata.I)

#Permanova II
adonis2(unifrac.dist.II~SampleType, permutations=999, data=metadata.II)

#Permanova III
adonis2(unifrac.dist.III~SampleType, permutations=999, data=metadata.III)

#Permanova IV
adonis2(unifrac.dist.IV~SampleType, permutations=999, data=metadata.IV)

#Permanova V
adonis2(unifrac.dist.V~SampleType, permutations=999, data=metadata.V)

#Permanova VI
adonis2(unifrac.dist.VI~SampleType, permutations=999, data=metadata.VI)

#Permanova VII
adonis2(unifrac.dist.VII~SampleType, permutations=999, data=metadata.VII)

#Permanova VIII
adonis2(unifrac.dist.VIII~SampleType, permutations=999, data=metadata.VIII)

#Permanova IX
adonis2(unifrac.dist.IX~SampleType, permutations=999, data=metadata.IX)

    #None of the Permanovas significant (0.2 or 0.25), despite high R2 (0.81-0.98).
    #NUmber of permutations: 119 -> why not 119?

#################################################################################
#Permanovas on entire data set

adonis2(unifrac.dist.msca~Treatment*SampleType, permutations = 999, data=metadata1)

    #Treatment (R2: 0.43), SampleType (R2: 0.43), and Interaction (R2: 0.11) highly significant (0.0001)

#################################################################################
#Remove rare taxa

#work up relabund df.
relabund1=relabund[,-234:-236] #remove taxonomy columns
rownames(relabund1)=relabund1$OTUNumber #update rownames
relabund2=relabund1[,-1] #remove otu names column
relabund2.t=as.data.frame(t(relabund2)) #transpose relabund2

#cull otu
cull.otu(relabund2.t,3,.005,.01) #451 ASVs remain.

#visualize ASV idstibution.
hist(unlist(relabund2.t)) #with zeros
hist(unlist(relabund2.t)[unlist(relabund2.t)!=0]) #without zeros
hist(unlist(relabund.df.cull)[unlist(relabund.df.cull)!=0]) #after culling

#save relabund.df.cull as new name.
relabund2.t.cull=relabund.df.cull

#################################################################################
#plot a stacked barchart of the innoculum and t0 bottle

#First, prep data frames for long format script.

#work up relabund2.t.cull
relabund2.cull=as.data.frame(t(relabund2.t.cull)) #transpose

#work up taxonomy df
taxonomy=relabund[,234:236] #extract taxonomy coumns.
taxonomy1=cbind(taxonomy,t(as.data.frame(strsplit(as.character(taxonomy$OTUConTaxonomy),split=";")))) #split contaxonomy by ";", add to taxonomy df
rownames(taxonomy1)=taxonomy1$repSeqName #change rownames
colnames(taxonomy1)[c(1,4:9)]=c("OTUNumber","Domain","Phylum","Class","Order","Family","Genus") #update colnames
taxonomy1.cull=taxonomy1[taxonomy1$OTUNumber %in% rownames(relabund2.cull),] #subset for only culled ASVs


generate.long.format(relabund2.cull,metadata.t0.v1,taxonomy1.cull)

#export as csv, generate stacked bar charts in JMP
write.csv(abund.longformat,file="abund.longformat.csv")

#################################################################################

