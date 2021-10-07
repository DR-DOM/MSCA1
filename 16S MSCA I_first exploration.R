#16S MSCA I first exploration

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
metadata=read.csv(file="metadata.csv")

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
  geom_point(size=5)+
  geom_text(aes(label=SampleName),hjust=1.5,vjust=0)
        #Mostly clear separation between inoculums (right) and bioassays left
            #Exceptions:TurfNight (middle up), SWnight (with TurfDay & TurfNight bioassays), and SWday (with SWnight bioassays)
        #Bioassays cluster generally together per "Treatment", but no clear pattern regarding "Organism" or "Time"
            #Exceptions: One SW night close to inoculums -> Mixing up SWnight_Bioassay with SWday_Inoculum?
            #Exceptions: One TurfDay horizontally alligning with TurfNight_Inoculum and SWnight inoculum clustering with TurfDay -> Mixing up?

    #Next step adjust metadata table and re-run
        # -> switch VIII_I7 with VII_B4
        # -> switch VIII_B4 with VII_I6
        # -> remove IX_B4

#plot nmds Inoculum
ggplot(nmds.inoculum.scores1,aes(x=NMDS1,y=NMDS2,color=Species,shape=Experiment))+
  geom_point(size=6)+
  geom_text(aes(label=Species),hjust=1.5,vjust=0)
        #NMDS warning for low stress: Most inoculums cluster together, I6 and I7 far off, Turf night in between

#plot nmds Bioassay
ggplot(nmds.bioassay.scores1,aes(x=NMDS1,y=NMDS2,color=Species,shape=Experiment))+
  geom_point(size=6)+
  geom_text(aes(label=SampleName),hjust=1.5,vjust=0)
        #Experiment 1 more spread out than Experiment 2
        #No clear separation coral vs algae
        #turfs and SWnight appear to cluster together
        #few inoculum-bioassay mix ups will be corrected

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

    #None of the Permanovas significant (0.2), despite high R2 (95-98).
    # VIII and IX low p-value (0.4) and R2-values (14-31), but should improve after corrections

#################################################################################
#Permanovas on entire data set

adonis2(unifrac.dist.msca~Treatment*SampleType, permutations = 999, data=metadata1)

    #Treatment (R2: 0.41), SampleType (R2: 0.18), and Interaction (R2: 0.24) highly significant (0.0001)


