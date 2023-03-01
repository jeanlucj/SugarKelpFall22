setwd("~/Documents/GitRepo/SugarKelpBreedingMao/SugarKelpBreeding/TraitAnalyses201003/GS_ms")
allBLUPs0 <- read.csv("./data_in_GS_ms/AllBLUPs_SP_GP_Supp Fig 3_inputdata.csv", header=T)

datafdr <- "./output_in_GS_ms/"
#### Supp Fig 3 BV plots vs different category
data<-as.data.frame(allBLUPs0)  # This is the AllBLUPs excel file attached
data$Trait<-data$DWpM

yr<-"Yr1920_AddYr21BV"
diphap<-"MixedPloidy"
traitname<-"DWpM"
library(ggplot2)
pdf(paste0(datafdr,"Supp Fig 2",traitname,"vs Individual",yr,"_PhotoScore23_with_",diphap,"_Supp_Fig_2_062022.pdf"))

ggplot(data=data,mapping=aes(x=1:nrow(data),y=Trait,color=Category,shape=Group))+
  geom_point()+
  scale_color_manual(name="Category",
                     values=c("Founders"="black",
                              "Founder_SP_without_progeny"="grey",
                              "GP_from_Founders"="darkgoldenrod2",
                              "SP_2019"="blue1",

                              "SP_2020"="purple",

                              "SP_2021"="brown",

                              "GP_from_Farm"="tomato1"))+
  scale_shape_manual(values=c("SPs"=15,"SPs_used_GP_from_Farm"=3,"GPs"=17))+
  labs(x="Individual",y=paste0("BVs for ",traitname))+
  theme_bw()
dev.off()

p <- ggplot(allBLUPs0, aes(x=Category, y=DWpM)) +
  geom_boxplot(notch=T)
p
