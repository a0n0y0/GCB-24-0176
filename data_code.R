####
####
####
####
####
setwd("F:/JGS-SMT微生物测序数据/R分析/Abundance")
###
BacPly <- read.csv("JGSbac_All_others.csv",header = T,check.names=F)
head(BacPly)
## 
Bac_spe <- BacPly[-c(1:11)]
head(Bac_spe)
tax <- BacPly[1:11]
head(tax)
##
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F)
head(group)
group$treat <- factor(group$treat,levels = c("CK","CALN","CAHN","UALN","UAHN"))
group$Sampling_time <- factor(group$Sampling_time,levels = c("Jan., 2020","July, 2020","Jan., 2021","July, 2021"))
str(group)
##
## 
library(tidyverse)
phy <- Bac_spe %>%
  group_by(tax$Phylum) %>% #
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="Sample",value = "abun",-Phylum) %>% #
  left_join(group,by=c("Sample"="Sample1")) %>%
  select(treat,Sampling_time,Phylum,abun) %>%
  group_by(treat,Sampling_time,Phylum) %>% #
  summarise_all(mean)
dim(phy)
head(phy)
#
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.7)(20)
col3 = pal_d3("category20",alpha = 0.5)(4)
mypal=c(col,col2,col3)
mypal
library(ggplot2)
B <- ggplot()+
  geom_bar(data=phy,
           aes(x=treat,
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "fill",
           width=0.5)+
  facet_grid(.~Sampling_time)+
  scale_fill_manual(values = mypal)+
  theme_bw()+
  scale_y_continuous(
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  guides(fill=guide_legend(title = "Phylum",ncol = 1))+
  labs(x=NULL)+
  theme(legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 12,colour = "black",family = "serif"))+
  theme(axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.text.y = element_text(face = "bold", 
                                   size = 10,color="black",family = "serif"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black",family = "serif"))
B
#axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
ggsave('Bac_Abundance_others.pdf',B,width=8,height=6)
ggsave("Bac_Abundance_others.TIFF",B,width = 8,height = 6)
###
###
library(tidyverse)
phyO <- Bac_spe %>%
  group_by(tax$Phylum) %>% #
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="Sample",value = "abun",-Phylum) %>% # 
  left_join(group,by=c("Sample"="Sample1")) %>%
  select(treat,season_2,Phylum,abun) %>%
  group_by(treat,season_2,Phylum) %>% # 
  summarise_all(mean)
dim(phyO)
head(phyO)
#颜色
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.7)(20)
col3 = pal_d3("category20",alpha = 0.5)(4)
mypal=c(col,col2,col3)
mypal

#
library(ggplot2)
#pdf("abs_stack.pdf",width = 12,height = 10,family="Times")
D <- ggplot()+
  geom_bar(data=phyO,
           aes(x=treat,
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "fill",# 
           width=0.5)+
  facet_grid(.~season_2)+
  scale_fill_manual(values = mypal)+
  theme_bw()+
  scale_y_continuous(#
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  guides(fill=guide_legend(title = "Phylum",ncol = 1))+
  labs(x=NULL)+
  theme(legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 12,colour = "black",family = "serif"))+
  theme(axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.text.y = element_text(face = "bold", 
                                   size = 10,color="black",family = "serif"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black",family = "serif"))
D
#axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
ggsave('Bac_Abundance_O_others.pdf',D,width=4.05,height=6)
ggsave("Bac_Abundance_O_others.TIFF",D,width = 4.05,height = 6)

#####Family#####
#
setwd("F:/JGS-SMT微生物测序数据/R分析/Abundance")
###
BacPly <- read.csv("JGSbac_All.csv",header = T,check.names=F)
head(BacPly)
## 
Bac_spe <- BacPly[-c(1:11)]
head(Bac_spe)
tax <- BacPly[1:11]
head(tax)
#
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F)
head(group)

group$treat <- factor(group$treat,`levels`  = c("CK","CALN","CAHN","UALN","UAHN"))
group$Sampling_time <- factor(group$Sampling_time,levels = c("Jan., 2020","July, 2020","Jan., 2021","July, 2021"))
str(group)
#
library(tidyverse)
Family01 <- Bac_spe %>%
  group_by(tax$Family) %>% 
  summarise_all(sum) %>%
  rename(Family = `tax$Family`)
Family01
write.csv(Family01,'Bac_Family_data_site.csv')
#
Family01 <- t(Family01)
Family01
write.csv(Family01,'Bac_Family_data_site01.csv')
#
data1 <- read.csv("Bac_Family_data_site.csv",row.names=1,header = T,check.names=F)
data2 <- t(data1)
data2
phy03 <- data2/rowSums(data2)
phy03
rowSums(phy03)
write.csv(phy03,'Bac_Family_data_relative.csv')
##
#samples
library(tidyverse)
dat1 <- read.csv("Bac_Family_data_site01.csv",header = T,check.names=F)
dat_Bac <- merge(group,dat1 ,by="Sample")
dat_Bac
write.csv(dat_Bac,'Bac_Family_group.csv')
#
library(tidyverse)
Family <- Bac_spe %>%
  group_by(tax$Family) %>% #
  summarise_all(sum) %>%
  rename(Family = `tax$Family`) %>%
  gather(key="Sample",value = "abun",-Family) %>% #
  left_join(group,by=c("Sample"="Sample")) %>%
  select(treat,Sampling_time,Family,abun) %>%
  group_by(treat,Sampling_time,Family) %>% # 
  summarise_all(mean)
dim(Family)
head(Family)
write.csv(Family,'Bac_Family_data.csv')
#
library(RColorBrewer)
display.brewer.all()##
mycol <- brewer.pal(n=5,name='Set1')
mycol
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.7)(20)
col3 = pal_d3("category20",alpha = 0.5)(4)
mypal=c(col,col2,col3)
mypal
#
##
library(reshape2)
library(ggplot2)
library(vegan)
library(picante)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
#
shapiro.test(data2$f__Acetobacteraceae)
data2 %>% select(treat,f__Acetobacteraceae) %>% group_by(treat) %>%
  summarise_all(.funs=funs(statistic=shapiro.test(.)$statistic,
                           p.value=shapiro.test(.)$p.value))###

display.brewer.all()##颜色面板
mycol <- brewer.pal(n=5,name='Set1')#
mycol
col1 <- brewer.pal(n=8,name='Greys')##
col1
#"#FFFFFF" "#F0F0F0" "#D9D9D9" "#BDBDBD" "#969696" "#737373" "#525252" "#252525"
#Shannon
#Overview
data2$treat<-factor(data2$treat,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_1 <- ggboxplot(data2,x="treat",y="f__norank_o__Vicinamibacterales",color="treat",
                 palette=mycol,legend="NA",#
                 add="jitter",#
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="season_2")+
  stat_compare_means(method = "kruskal.test",label.y = 5,label.x = 0.8,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='wilcox',label="p.signif")+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
A_1
###
write.csv(data2,'Bac_Family_t_test.csv')#
#
library(tidyverse)
data3 <- data2 %>%  rows_delete(tibble(Approach="CK"))
data3
library(stats) 
#Shapiro-Wilk
#1
with(data3, shapiro.test(f__Acidothermaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Acidothermaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Acidothermaceae~Approach, data=data3)
##t
t.test(f__Acidothermaceae~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) #
#2
with(data3, shapiro.test(f__Solirubrobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Solirubrobacteraceae[Approach == "UAN"]))
#levene方
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Solirubrobacteraceae~Approach, data=data3)
##t
wilcox.test(f__Solirubrobacteraceae~Approach, data=data3, var.equal=T) #
#3
with(data3, shapiro.test(f__Xanthobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Xanthobacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Xanthobacteraceae~Approach, data=data3)
##t
t.test(f__Xanthobacteraceae~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) #
#4
with(data3, shapiro.test(f__Burkholderiaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Burkholderiaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Burkholderiaceae~Approach, data=data3)
##
wilcox.test(f__Burkholderiaceae~Approach, data=data3, var.equal=T) #
#5
with(data3, shapiro.test(f__norank_o__Subgroup_2[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__Subgroup_2[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__Subgroup_2~Approach, data=data3)
##检验
t.test(f__norank_o__Subgroup_2~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) 
#6
with(data3, shapiro.test(f__norank_o__Gaiellales[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__Gaiellales[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__Gaiellales~Approach, data=data3)
##
wilcox.test(f__norank_o__Gaiellales~Approach, data=data3, var.equal=T) #
#7
with(data3, shapiro.test(f__Acidobacteriaceae_Subgroup_1[Approach == "CAN"]))
with(data3, shapiro.test(f__Acidobacteriaceae_Subgroup_1[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Acidobacteriaceae_Subgroup_1~Approach, data=data3)
##
t.test(f__Acidobacteriaceae_Subgroup_1~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) #
#8
with(data3, shapiro.test(f__norank_o__Elsterales[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__Elsterales[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__Elsterales~Approach, data=data3)
##
wilcox.test(f__norank_o__Elsterales~Approach, data=data3, var.equal=T) #
#9
with(data3, shapiro.test(f__unclassified_o__Gammaproteobacteria_Incertae_Sedis[Approach == "CAN"]))
with(data3, shapiro.test(f__unclassified_o__Gammaproteobacteria_Incertae_Sedis[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__unclassified_o__Gammaproteobacteria_Incertae_Sedis~Approach, data=data3)
##
t.test(f__unclassified_o__Gammaproteobacteria_Incertae_Sedis~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) #
#10
with(data3, shapiro.test(f__norank_o__Acidobacteriales[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__Acidobacteriales[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__Acidobacteriales~Approach, data=data3)
##
t.test(f__norank_o__Acidobacteriales~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) 
#11
with(data3, shapiro.test(f__Rhodanobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Rhodanobacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Rhodanobacteraceae~Approach, data=data3)
##检验
wilcox.test(f__Rhodanobacteraceae~Approach, data=data3, var.equal=T) 
#12
with(data3, shapiro.test(f__norank_o__norank_c__norank_p__WPS_2[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__norank_c__norank_p__WPS_2[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__norank_c__norank_p__WPS_2~Approach, data=data3)
##
wilcox.test(f__norank_o__norank_c__norank_p__WPS_2~Approach, data=data3, var.equal=T) 
#13
with(data3, shapiro.test(f__norank_o__norank_c__Acidimicrobiia[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__norank_c__Acidimicrobiia[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__norank_c__Acidimicrobiia~Approach, data=data3)
##
t.test(f__norank_o__norank_c__Acidimicrobiia~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) 
#14
with(data3, shapiro.test(f__norank_o__norank_c__AD3[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__norank_c__AD3[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__norank_c__AD3~Approach, data=data3)
##
wilcox.test(f__norank_o__norank_c__AD3~Approach, data=data3, var.equal=T) 
#15
with(data3, shapiro.test(f__Gemmatimonadaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Gemmatimonadaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Gemmatimonadaceae~Approach, data=data3)
##
wilcox.test(f__Gemmatimonadaceae~Approach, data=data3, var.equal=T) 
#16
with(data3, shapiro.test(f__Gemmataceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Gemmataceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Gemmataceae~Approach, data=data3)
##
wilcox.test(f__Gemmataceae~Approach, data=data3, var.equal=T) 
#17
with(data3, shapiro.test(f__Micropepsaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Micropepsaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Micropepsaceae~Approach, data=data3)
##
wilcox.test(f__Micropepsaceae~Approach, data=data3, var.equal=T) 
#18
with(data3, shapiro.test(f__Ktedonobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Ktedonobacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Ktedonobacteraceae~Approach, data=data3)
##
wilcox.test(f__Ktedonobacteraceae~Approach, data=data3, var.equal=T) 
#19
with(data3, shapiro.test(f__Beijerinckiaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Beijerinckiaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Beijerinckiaceae~Approach, data=data3)
##
wilcox.test(f__Beijerinckiaceae~Approach, data=data3, var.equal=T)
#20
with(data3, shapiro.test(f__unclassified_k__norank_d__Bacteria[Approach == "CAN"]))
with(data3, shapiro.test(f__unclassified_k__norank_d__Bacteria[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__unclassified_k__norank_d__Bacteria~Approach, data=data3)
##
wilcox.test(f__unclassified_k__norank_d__Bacteria~Approach, data=data3, var.equal=T)
#21
with(data3, shapiro.test(f__Mycobacteriaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Mycobacteriaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Mycobacteriaceae~Approach, data=data3)
##
wilcox.test(f__Mycobacteriaceae~Approach, data=data3, var.equal=T)
#22
with(data3, shapiro.test(f__Acetobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Acetobacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Acetobacteraceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Acetobacteraceae~Approach, data=data3, var.equal=T) #
#23
with(data3, shapiro.test(f__Pedosphaeraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Pedosphaeraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Pedosphaeraceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Pedosphaeraceae~Approach, data=data3, var.equal=T) #
#24
with(data3, shapiro.test(f__Caulobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Caulobacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Caulobacteraceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Caulobacteraceae~Approach, data=data3, var.equal=F) #
#25
with(data3, shapiro.test(f__SC_I_84[Approach == "CAN"]))
with(data3, shapiro.test(f__SC_I_84[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__SC_I_84~Approach, data=data3)
##wilcoxon
wilcox.test(f__SC_I_84~Approach, data=data3, var.equal=F) #
#26
with(data3, shapiro.test(f__norank_o__IMCC26256[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__IMCC26256[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__IMCC26256~Approach, data=data3)
##wilcoxon
t.test(f__norank_o__IMCC26256~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) 
#27
with(data3, shapiro.test(f__Sphingobacteriaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Sphingobacteriaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Sphingobacteriaceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Sphingobacteriaceae~Approach, data=data3, var.equal=F) #
with(data3, shapiro.test(f__Isosphaeraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Isosphaeraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Isosphaeraceae~Approach, data=data3)
##wilcoxon
t.test(f__Isosphaeraceae~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) 
#29
with(data3, shapiro.test(f__Bryobacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Bryobacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Bryobacteraceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Bryobacteraceae~Approach, data=data3, var.equal=T) #
#30
with(data3, shapiro.test(f__Pseudonocardiaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Pseudonocardiaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Pseudonocardiaceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Pseudonocardiaceae~Approach, data=data3, var.equal=T) #
#31
with(data3, shapiro.test(f__WD2101_soil_group[Approach == "CAN"]))
with(data3, shapiro.test(f__WD2101_soil_group[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__WD2101_soil_group~Approach, data=data3)
##wilcoxon
wilcox.test(f__WD2101_soil_group~Approach, data=data3, var.equal=T) #
#32
with(data3, shapiro.test(f__Solibacteraceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Solibacteraceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Solibacteraceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Solibacteraceae~Approach, data=data3, var.equal=T) #
#33
with(data3, shapiro.test(f__Paenibacillaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Paenibacillaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Paenibacillaceae~Approach, data=data3)
##wilcoxon
wilcox.test(f__Paenibacillaceae~Approach, data=data3, var.equal=T) #
with(data3, shapiro.test(f__Pirellulaceae[Approach == "CAN"]))
with(data3, shapiro.test(f__Pirellulaceae[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__Pirellulaceae~Approach, data=data3)
##wilcoxon
t.test(f__Pirellulaceae~Approach, data=data3, paired = FALSE, alternative = 'two.sided', var.equal=T) 
#35
with(data3, shapiro.test(f__norank_o__Vicinamibacterales[Approach == "CAN"]))
with(data3, shapiro.test(f__norank_o__Vicinamibacterales[Approach == "UAN"]))
#levene
library(car)
data3$Approach <- as.factor(data3$Approach)
leveneTest(f__norank_o__Vicinamibacterales~Approach, data=data3)
##wilcoxon
wilcox.test(f__norank_o__Vicinamibacterales~Approach, data=data3, var.equal=T) #
#
#
#

getwd()
setwd("F:/JGS-SMT微生物测序数据/R分析/Abundance")
###
data1_1 <- read.csv("Fun_Family_data_relative_单个科.csv",header = T,check.names=F)
head(data1_1)
##
group <- read.csv("jgsgroupFuntotalNew.csv",header = T,check.names=F)
head(group)
#
group$treat <- factor(group$treat,`levels`  = c("CK","CALN","CAHN","UALN","UAHN"))
group$Sampling_time <- factor(group$Sampling_time,levels = c("Jan., 2020","July, 2020","Jan., 2021","July, 2021"))
str(group)
#
data2 <- merge(group,data1_1 ,by="Sample")
data2
head(data2)
##
######
######
######
######
######
setwd("F:/JGS-SMT微生物测序数据/R分析/Abundance")
funPly <- read.csv("JGSfun_All.csv",header = T,check.names=F)
head(funPly)
## 
fun_spe <- funPly[-c(1:11)]
head(fun_spe)
tax <- funPly[1:11]
head(tax)
###
group <- read.csv("jgsgroupFuntotalNew.csv",header = T,check.names=F)
head(group)
group$treat <- factor(group$treat,`levels`  = c("CK","CALN","CAHN","UALN","UAHN"))
group$Sampling_time <- factor(group$Sampling_time,levels = c("Jan., 2020","July, 2020","Jan., 2021","July, 2021"))
str(group)
#phy01 <- fun_spe %>%
  group_by(tax$Phylum) %>% # 
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`)
phy01
write.csv(phy01,'Fun_phy_data_site.csv')
#
phy01 <- t(phy01)
phy01
write.csv(phy01,'Fun_phy_data_site01.csv')
#
data1 <- read.csv("Fun_phy_data_site.csv",row.names=1,header = T,check.names=F)
data2 <- t(data1)
data2
phy03 <- data2/rowSums(data2)
phy03
rowSums(phy03)
write.csv(phy03,'Fun_phy_data_relative.csv')

###
FunPly <- read.csv("JGSfun_All_others.csv",header = T,check.names=F)
head(FunPly)
## 
Fun_spe <- FunPly[-c(1:11)]
head(Fun_spe)
tax <- FunPly[1:11]
head(tax)
##
group <- read.csv("jgsgroupFuntotalNew.csv",header = T,check.names=F)
head(group)
group$treat <- factor(group$treat,levels = c("CK","CALN","CAHN","UALN","UAHN"))
group$Sampling_time <- factor(group$Sampling_time,levels = c("Jan., 2020","July, 2020","Jan., 2021","July, 2021"))
str(group)
##
## spe和tax
library(tidyverse)
Funphy <- Fun_spe %>%
  group_by(tax$Phylum) %>% # 
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="Sample",value = "abun",-Phylum) %>% # 
  left_join(group,by=c("Sample"="Sample1")) %>%
  select(treat,Sampling_time,Phylum,abun) %>%
  group_by(treat,Sampling_time,Phylum) %>% # 
  summarise_all(mean)
dim(Funphy)
head(Funphy)
#
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.7)(20)
col3 = pal_d3("category20",alpha = 0.5)(4)
mypal=c(col,col2,col3)
mypal
#pdf("abs_stack.pdf",width = 12,height = 10,family="Times")
#dev.off()#
library(ggplot2)
#pdf("abs_stack.pdf",width = 12,height = 10,family="Times")
E <- ggplot()+
  geom_bar(data=Funphy,
           aes(x=treat,
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "fill",# 
           width=0.5)+
  facet_grid(.~Sampling_time)+
  scale_fill_manual(values = mypal)+
  theme_bw()+
  scale_y_continuous(
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  guides(fill=guide_legend(title = "Phylum",ncol = 1))+
  labs(x=NULL)+
  theme(legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 12,colour = "black",family = "serif"))+
  theme(axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.text.y = element_text(face = "bold", 
                                   size = 10,color="black",family = "serif"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black",family = "serif"))
E
#axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
ggsave('Fun_Abundance_others.pdf',E,width=8,height=6)
ggsave("Fun_Abundance_others.TIFF",E,width = 8,height = 6)
###

###
library(tidyverse)
phyO <- Fun_spe %>%
  group_by(tax$Phylum) %>% # 
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="Sample",value = "abun",-Phylum) %>% # 
  left_join(group,by=c("Sample"="Sample1")) %>%
  select(treat,season_4,Phylum,abun) %>%
  group_by(treat,season_4,Phylum) %>% # 
  summarise_all(mean)
dim(phyO)
head(phyO)
#
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.7)(20)
col3 = pal_d3("category20",alpha = 0.5)(4)
mypal=c(col,col2,col3)
mypal

#
library(ggplot2)
#pdf("abs_stack.pdf",width = 12,height = 10,family="Times")
G <- ggplot()+
  geom_bar(data=phyO,
           aes(x=treat,
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "fill",# 
           width=0.5)+
  facet_grid(.~season_4)+
  scale_fill_manual(values = mypal)+
  theme_bw()+
  scale_y_continuous(#expand=c(0,0), 
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  guides(fill=guide_legend(title = "Phylum",ncol = 1))+
  labs(x=NULL)+
  theme(legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 12,colour = "black",family = "serif"))+
  theme(axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.text.y = element_text(face = "bold", 
                                   size = 10,color="black",family = "serif"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black",family = "serif"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black",family = "serif"))
G
#axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
ggsave('Fun_Abundance_O_others.pdf',G,width=5,height=6)
ggsave("Fun_Abundance_O_others.TIFF",G,width = 4.2,height = 6)
##JGS##
##JGS##
##JGS##
##############JGSalpha#####
#################################
###
getwd()
setwd("F:/JGS-SMT微生物测序数据/R分析/20230119Bacteria")
library(reshape2)
library(ggplot2)
library(vegan)
library(picante)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
Index <- read.csv("BacAlphLMM.csv",header = T,row.names=1,check.names=F)
Index$samples <- rownames(Index)###
groups <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F)
colnames(groups)[1:3] <- c('samples','group','approach')
df <- merge(Index,groups,by='samples')
df
#write.csv(df,'df_new.csv')
shapiro.test(df$Shannon)
df %>% select(group,Shannon) %>% group_by(group) %>%
  summarise_all(.funs=funs(statistic=shapiro.test(.)$statistic,
                           p.value=shapiro.test(.)$p.value))###

display.brewer.all()##
mycol <- brewer.pal(n=5,name='Set1')##也可手动设置颜色mycol <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00')
mycol
col1 <- brewer.pal(n=8,name='Greys')##也可手动设置颜色mycol <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00')
col1
#"#FFFFFF" "#F0F0F0" "#D9D9D9" "#BDBDBD" "#969696" "#737373" "#525252" "#252525"
#Shannon
#Overview
df$group<-factor(df$group,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_1 <- ggboxplot(df,x="group",y="Shannon",color="group",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="season_2.y")+
  stat_compare_means(method = "anova",label.y = 8.55,label.x = 0.8,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
A_1
A_1 <- A_1+labs(x=NULL,y="Shannon index")
A_1
A_1 <- A_1+scale_y_continuous(limits=c(5.7,9))
A_1
A_1 <- A_1 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
A_1
#A <- p_1+labs(fill="Treatment")##
#A <- p_1+scale_fill_discrete(name="Condition")##
#Sampling_time#
#Sampling_time#
#Sampling_time#
#Sampling_time#
#Sampling_time#
getwd()
setwd("F:/JGS-SMT微生物测序数据/R分析/20230119Bacteria")
library(reshape2)
library(ggplot2)
library(vegan)
library(picante)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
df <- read.csv("BacAlphLMM.csv",header = T,row.names=1,check.names=F)
df$samples <- rownames(df)###
df_1 <- df %>% filter(Sampling_time=="Jan., 2020")
df_1
df_2 <- df %>% filter(Sampling_time=="July, 2020")
df_2
df_3 <- df %>% filter(Sampling_time=="Jan., 2021")
df_3
df_4 <- df %>% filter(Sampling_time=="July, 2021")
df_4
#write.csv(df,'df_new.csv')
shapiro.test(df$Shannon)
df %>% select(Treatment,Shannon) %>% group_by(Treatment) %>%
  summarise_all(.funs=funs(statistic=shapiro.test(.)$statistic,
                           p.value=shapiro.test(.)$p.value))###

display.brewer.all()##
mycol <- brewer.pal(n=5,name='Set1')##
mycol
col1 <- brewer.pal(n=8,name='Greys')##
col1
#"#FFFFFF" "#F0F0F0" "#D9D9D9" "#BDBDBD" "#969696" "#737373" "#525252" "#252525"
#Shannon
#Overview
df$Treatment<-factor(df$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_1 <- ggboxplot(df,x="Treatment",y="Shannon",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=2,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="season_2")+
  stat_compare_means(method = "anova",label.y = 8.55,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
A_1
A_1 <- A_1+labs(x=NULL,y="Shannon index")
A_1
A_1 <- A_1+scale_y_continuous(limits=c(5.7,9))
A_1
A_1 <- A_1 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
A_1


#Jan., 2020
df_1$Treatment<-factor(df_1$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_2 <- ggboxplot(df_1,x="Treatment",y="Shannon",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=2,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 8.55,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")),#
  )#
A_2
A_2 <- A_2+labs(x=NULL,y=NULL)
A_2
A_2 <- A_2+scale_y_continuous(limits=c(5.7,9))
A_2

A_2 <- A_2 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
A_2

#July, 2020
df_2$Treatment<-factor(df_2$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_3 <- ggboxplot(df_2,x="Treatment",y="Shannon",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=2,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 8.55,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
A_3
A_3 <- A_3+labs(x=NULL,y=NULL)
A_3
A_3 <- A_3+scale_y_continuous(limits=c(5.7,9))
A_3
A_3 <- A_3 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
A_3
#
#Jan., 2021
df_3$Treatment<-factor(df_3$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_4 <- ggboxplot(df_3,x="Treatment",y="Shannon",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=2,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 8.55,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
A_4
A_4 <- A_4+labs(x=NULL,y=NULL)
A_4
A_4 <- A_4+scale_y_continuous(limits=c(5.7,9))
A_4
A_4 <- A_4 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
A_4
#
#
#July, 2021
df_4$Treatment<-factor(df_4$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
A_5 <- ggboxplot(df_4,x="Treatment",y="Shannon",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=2,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 8.55,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
A_5
A_5 <- A_5+labs(x=NULL,y=NULL)
A_5
A_5 <- A_5+scale_y_continuous(limits=c(5.7,9))
A_5
A_5 <- A_5 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
A_5


library(patchwork)
A <- A_1+A_2+A_3+A_4+A_5+plot_layout(guides="collect",ncol=5,nrow=1)&
  theme(legend.position = "top",
        legend.title=element_text(size=15,face="bold",colour="black",family="serif"),
        legend.text = element_text(family = "serif",colour = "black", size = 14))


A
ggsave('ABS01.pdf',A,width=10,height=4)
ggsave("ABS01.TIFF",A,width = 10,height = 4)
#Richness
#Richness
#Richness
df$Treatment<-factor(df$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
B_1 <- ggboxplot(df,x="Treatment",y="Richness",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="season_2")+
  stat_compare_means(method = "anova",label.y = 6000,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
B_1
B_1 <- B_1+labs(x=NULL,y="Richness index")
B_1
B_1 <- B_1+scale_y_continuous(limits=c(1000,6500))
B_1
B_1 <- B_1 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
B_1


#Jan., 2020
df_1$Treatment<-factor(df_1$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
B_2 <- ggboxplot(df_1,x="Treatment",y="Richness",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#添加数据点，#shape="group",
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 6000,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")),#
  )#
B_2
B_2 <- B_2+labs(x=NULL,y=NULL)
B_2
B_2 <- B_2+scale_y_continuous(limits=c(1000,6500))
B_2

B_2 <- B_2 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
B_2

#July, 2020
df_2$Treatment<-factor(df_2$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
B_3 <- ggboxplot(df_2,x="Treatment",y="Richness",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#添加数据点，#shape="group",
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 6000,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
B_3
B_3 <- B_3+labs(x=NULL,y=NULL)
B_3
B_3 <- B_3+scale_y_continuous(limits=c(1000,6500))
B_3
B_3 <- B_3 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
B_3
#
#Jan., 2021
df_3$Treatment<-factor(df_3$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
B_4 <- ggboxplot(df_3,x="Treatment",y="Richness",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 6000,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
B_4
B_4 <- B_4+labs(x=NULL,y=NULL)
B_4
B_4 <- B_4+scale_y_continuous(limits=c(1000,6500))
B_4
B_4 <- B_4 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
B_4
#
#
#July, 2021
df_4$Treatment<-factor(df_4$Treatment,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
B_5 <- ggboxplot(df_4,x="Treatment",y="Richness",color="Treatment",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=3,jitter=0.1,alpha=0.5),outlaier.shape=NA,
                 facet.by="Sampling_time")+
  stat_compare_means(method = "anova",label.y = 6000,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('UALN','UAHN')),
                     method='t.test',label="p.signif")+###
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,angle =45,family = "serif"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=1,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))##
B_5
B_5 <- B_5+labs(x=NULL,y=NULL)
B_5
B_5 <- B_5+scale_y_continuous(limits=c(1000,6500))
B_5
B_5 <- B_5 + theme(legend.position="right",
                   legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
B_5


library(patchwork)
B <- B_1+B_2+B_3+B_4+B_5+plot_layout(guides="collect",ncol=5,nrow=1)&
  theme(legend.position = "top",
        legend.title=element_text(size=15,face="bold",colour="black",family="serif"),
        legend.text = element_text(family = "serif",colour = "black", size = 14))
B
ggsave('BBR.pdf',B,width=10,height=4)
ggsave("BBR.TIFF",B,width = 10,height = 4)

library(patchwork)
BacOverview <- A_1+A_2+A_3+A_4+A_5+B_1+B_2+B_3+B_4+B_5+plot_layout(guides="collect",ncol=5,nrow=2)&
  theme(legend.position = "top",
        legend.title=element_text(size=15,face="bold",colour="black",family="serif"),
        legend.text = element_text(family = "serif",colour = "black", size = 14))
BacOverview
ggsave('BacOverview.pdf',BacOverview,width=10,height=4)
ggsave("BacOverview.TIFF",BacOverview,width = 12,height = 8)

###beta
###beta
###beta
#####beta
###PCoA
###PCoA
###PCoA
setwd("F:/JGS-SMT微生物测序数据/R分析/20230119Bacteria")
Bac_O <- read.csv("JGSBacASV_order.csv",header = T,row.names=1,check.names=F)
head(Bac_O)
Bac_O_t <-t(Bac_O)
library(vegan)
asv.dist_Bac_O <- vegdist(Bac_O_t)
asv.dist_Bac_O
Bac_O_pcoa <- cmdscale(asv.dist_Bac_O,eig=TRUE)
Bac_O_pcoa12 <- Bac_O_pcoa$points[,1:2]
Bac_O_pc_importance <- round(Bac_O_pcoa$eig/sum(Bac_O_pcoa$eig)*100,digits=2)
write.csv(Bac_O_pcoa12,'PC_Bac_order.csv')
library(ggplot2)
Bac_O_pcoa12 <- as.data.frame(Bac_O_pcoa12)
Bac_O_pcoa12$samples <- row.names(Bac_O_pcoa12)
head(Bac_O_pcoa12)
group_o <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F)
group_o
colnames(group_o) <- c("samples","group","group1","Season1","Season2")
head(group_o)
library(tidyverse)
df_O <- merge(Bac_O_pcoa12,group_o,by="samples")
df_O
#overview
N_treatment<-factor(df_O$group,levels = c("CK","CALN","CAHN","UALN","UAHN"))
N_approach<-factor(df_O$group1,levels = c("CK","CAN","UAN"))
library(showtext)
windowsFonts()###
Font <- c('STKaiti.TTF','simhei.TTF') #
for (i in Font) {
  font_path = i
  font_name = tools::file_path_sans_ext(basename(font_path))
  font_add(font_name, font_path)
}
font_families() ### 
showtext_auto(enable=TRUE) 
###
Bac_pc_O <- ggplot(df_O,aes(x=V1,y=V2,colour=N_approach,shape=N_treatment))+#
  geom_point(size=3)+#
  geom_text(aes(label=samples,y=V2+0.01),size=4,vjust=0)#
Bac_pc_O
###
Fon <- "serif"
Bac_pc_O <- ggplot(df_O,aes(x=V1,y=V2))+#
  geom_point(size=3,aes(shape=N_treatment,color=N_approach))#+#

Bac_pc_O
##
Bac_pc_O <- Bac_pc_O + stat_ellipse(linetype = 1,level = 0.95,aes(group  = N_approach, colour =  N_approach)) 
Bac_pc_O
Bac_pc_O <- Bac_pc_O+facet_wrap(~Season2 )#
Bac_pc_O
#
Bac_pc_O <- Bac_pc_O+theme(axis.title.x=element_text(size=15,family="serif"),
                           axis.title.y=element_text(size=15,family="serif",angle=90),#
                           axis.text.y=element_text(size=13,family="serif"),#
                           axis.text.x=element_text(size=13,family="serif"),#
                           legend.title = element_text(size=13,family="serif"),
                           legend.text = element_text(size=11,family="serif"),
                           strip.text = element_text(family = "serif",face='bold',size=13),
                           strip.background = element_rect(fill=c("#BDBDBD")))
Bac_pc_O
Bac_pc_O <- Bac_pc_O+xlab(paste0("PCo1(",Bac_O_pc_importance[1],"%)"))+ylab(paste0("PCo2(",Bac_O_pc_importance[2],"%)"))
Bac_pc_O
#
set.seed(1)
asv.div_O <- adonis2(Bac_O_t ~ group, data = group_o, permutations = 999, method="bray")
asv.div_O
###
set.seed(1)
asv.div_O1 <- adonis2(Bac_O_t ~ group1, data = group_o, permutations = 999, method="bray")
asv.div_O1
#
Bac_O_adonis <- paste0("N_treatment adonis R2: ",round(asv.div_O$R2,2), "; P: ", asv.div_O$`Pr(>F)`)
Bac_O_adonis1 <- paste0("N_approach adonis R2: ",round(asv.div_O1$R2,2), "; P: ", asv.div_O1$`Pr(>F)`)
# install.packages("ggalt")
#BiocManager::install("ggalt")
library(ggalt)
Bac_pc_O <- Bac_pc_O+labs(title=Bac_O_adonis,subtitle=Bac_O_adonis1)+
  theme(plot.title = element_text(family = "serif", #
                                  face = "bold",     #
                                  color = "black",      
                                  size = 12,          #
                                  hjust = 0,          #
                                  vjust = 0.5,          #
                                  angle = 0,          #
  ),
  plot.subtitle = element_text(family = "serif", #
                               face = "bold",     #
                               color = "black",      
                               size = 12,          
                               hjust = 0,          #
                               vjust = 0.5,          #
                               angle = 0, )
  )# 

Bac_pc_O

###保存图片
ggsave("PB_O_1_20230308.pdf",Bac_pc_O,width=5,height=4)
library(pairwiseAdonis)
# This is a wrapper function for multilevel pairwise comparison 
# using adonis() from package 'vegan'. 
# The function returns adjusted p-values using p.adjust().
#??pairwiseAdonis
data1.pairwise.adonis1 <- pairwise.adonis(x=Bac_O_t, factors=group_o$group, sim.function = "vegdist",
                                          sim.method = "bray",
                                          p.adjust.m = "BH",
                                          reduce = NULL,
                                          perm = 999)
data1.pairwise.adonis1
#结果导出
write.csv(data1.pairwise.adonis1,'BacBeta_Pairwise_treat.csv')
####
data1.pairwise.adonis2 <- pairwise.adonis(x=Bac_O_t, factors=group_o$group1, sim.function = "vegdist",
                                          sim.method = "bray",
                                          p.adjust.m = "BH",
                                          reduce = NULL,
                                          perm = 999)
data1.pairwise.adonis2
write.csv(data1.pairwise.adonis2,'BacBeta_Pairwise_Approa.csv')
#
library(ggpubr)
library(patchwork)
tab2 <- ggtexttable(data1.pairwise.adonis1[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                    theme = ttheme("blank")) %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
  tab_add_hline(at.row = nrow(data1.pairwise.adonis1)+1, row.side = "bottom", linewidth = 1)

tab3 <- ggtexttable(data1.pairwise.adonis2[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                    theme = ttheme("blank")) %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
  tab_add_hline(at.row = nrow(data1.pairwise.adonis2)+1, row.side = "bottom", linewidth = 1)

PB_O_Pair = PB_O + tab2 / tab3 + plot_layout(ncol = 2, heights=c(100,3,1))
PB_O_Pair
ggsave("PB_O_Pair.pdf",PB_O_Pair,width=12,height=6)
############
setwd("F:/JGS-SMT微生物测序数据/R分析/20230119Bacteria")
Bac <- read.csv("JGSBacASV.csv",header = T,row.names=1,check.names=F)
head(Bac)
SampleBac<-c("E1","M22","M20","I5","M15", "M6", "M4", "M5","A1" , "M1","M25","I1","M18", "M13", "A5",
             "A4","A13","E19","E28" ,"I25","I27","M8","I19","A11", "A8" ,"E27", "M11","M27","M19","A7",
             "E18", "A12", "A19", "I26", "I28", "M12", "M26", "A20", "E12", "E26", "I12","E4","E8","E11","M28",
             "A18" ,"A6",  "A14", "E20", "I8"  ,"I15" ,"I4"  ,"E22", "E25" ,"E6","E13" ,"I6" , "I13" ,"I11" ,"I7",
             "A15" ,"A21", "E15", "I22", "I18" ,"I14", "I21", "E5",  "E7","E14", "E21" ,"M7"  ,"M14" ,"M21", "I20")
##CK#
#CK
#CK
#CK
#CK
#CK
Bac_CK <- c("E1","M22","M20","I5","M15", "M6", "M4", "M5","A1" , "M1","M25","I1","M18", "M13", "A5")
Bac_CK <- Bac[,Bac_CK]
head(Bac_CK)
###对数据进行转置
Bac_CK_t <-t(Bac_CK)
###加载计算距离时所需的R包
library(vegan)
###计算bray_curtis距离
dist_Bac_CK <- vegdist(Bac_CK_t)
dist_Bac_CK
###PCoA分析
Bac_CK_pcoa <- cmdscale(dist_Bac_CK,eig=TRUE)
Bac_CK_pc12 <- Bac_CK_pcoa$points[,1:2]
Bac_CK_pc12
#write.csv(Bac_CK_pc12,"Bac_CK_pc12")
Bac_CK_importance <- round(Bac_CK_pcoa$eig/sum(Bac_CK_pcoa$eig)*100,digits=2)
##
library(ggplot2)
###pc12
Bac_CK_pc12 <- as.data.frame(Bac_CK_pc12)
###
Bac_CK_pc12$samples <- row.names(Bac_CK_pc12)
head(Bac_CK_pc12)
###
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F,row.names=1)
group
group$samples <- row.names(group)#
group
#
library(tidyverse)
Group_CK <- group %>% filter(samples==c("E1" , "M22", "M20" ,"I5" , "M15", "M6" , "M4" , "M5" , "A1" , "M1" , "M25", "I1" , "M18" ,"M13" ,"A5"))
Group_CK
Group_CK$Year <- as.factor(Group_CK$Year)
str(Group_CK)
###
colnames(Group_CK) <- c("treat","Approach","Season1","Season2","Sampling time","Year","samples")
head(Group_CK)
###
library(tidyverse)
df <- merge(Bac_CK_pc12,Group_CK,by="samples")
df
###
Season<-factor(df$Season1,levels = c("Dry season","Wet season"))
Year<-factor(df$Year,levels = c("2020","2021"))
####
library(showtext)
windowsFonts()###
#
Font <- c('STKaiti.TTF','simhei.TTF') #
for (i in Font) {
  font_path = i
  font_name = tools::file_path_sans_ext(basename(font_path))
  font_add(font_name, font_path)
}
font_families() ###
showtext_auto(enable=TRUE) #

###
Bac_PC_CK <- ggplot(df,aes(x=V1,y=V2,colour=Season,shape=factor(Year)))+#
  geom_point(size=3)+#
  geom_text(aes(label=samples,y=V2),size=4,vjust=0)#
Bac_PC_CK

###
Fon <- "serif"
Bac_PC_CK <- ggplot(df,aes(x=V1,y=V2))+#
  geom_point(size=3,aes(colour=Season,shape=factor(Year)))#+#
#stat_chull(geom="polygon",aes(group=N_approach,color=N_approach,fill=N_approach),alpha=0.1)#
Bac_PC_CK
##
Bac_PC_CK <- Bac_PC_CK + stat_ellipse(linetype = 1,level = 0.95,aes(group  = Season, colour =  Season)) 
Bac_PC_CK
Bac_PC_CK <- Bac_PC_CK+facet_wrap(~treat )#
Bac_PC_CK
#
Bac_PC_CK <- Bac_PC_CK+theme(axis.title.x=element_text(size=12,family="serif"),#
                             axis.title.y=element_text(size=12,family="serif",angle=90),#
                             axis.text.y=element_text(size=10,family="serif"),#
                             axis.text.x=element_text(size=10,family="serif"),#
                             legend.title = element_text(size=10,family="serif"),
                             legend.text = element_text(size=8,family="serif"),
                             strip.text = element_text(family = "serif",face='bold',size=13),
                             strip.background = element_rect(fill=c("#BDBDBD")))
Bac_PC_CK
Bac_PC_CK <- Bac_PC_CK+xlab(paste0("PCo1(",Bac_CK_importance[1],"%)"))+ylab(paste0("PCo2(",Bac_CK_importance[2],"%)"))
Bac_PC_CK

#??adonis2
set.seed(1)
asv.div_O <- adonis2(Bac_CK_t ~ Season1, data = Group_CK, permutations = 999, method="bray")
asv.div_O
##
set.seed(1)
asv.div_O1 <- adonis2(Bac_CK_t ~ Year, data = Group_CK, permutations = 999, method="bray")
asv.div_O1

####
Fun_O_adonis <- paste0("Season R2: ",round(asv.div_O$R2,2), "; P: ", asv.div_O$`Pr(>F)`)
Fun_O_adonis1 <- paste0("Year R2: ",round(asv.div_O1$R2,2), "; P: ", asv.div_O1$`Pr(>F)`)
# install.packages("ggalt")
#BiocManager::install("ggalt")
library(ggalt)
Bac_PC_CK <- Bac_PC_CK+labs(title=Fun_O_adonis,subtitle=Fun_O_adonis1)+
  theme(plot.title = element_text(family = "serif", #
                                  face = "bold",     #
                                  color = "black",      #
                                  size = 11,          #
                                  hjust = 0,          #
                                  vjust = 0,          #
                                  angle = 0,          #
  ),
  plot.subtitle = element_text(family = "serif", #
                               face = "bold",     #
                               color = "black",      #
                               size = 11,          #
                               hjust = 0,          #
                               vjust = 0,          #
                               angle = 0, )
  )#
Bac_PC_CK

###
ggsave("Bac_PC_CK.pdf",Bac_PC_CK,width=4,height=3)
ggsave("Bac_PC_CK.TIFF",Bac_PC_CK,width =4,height = 3)

#CALN
Bac_CALN <- c("A4","A13","E19","E28" ,"I25","I27","M8","I19","A11", "A8" ,"E27", "M11","M27","M19","A7")
Bac_CALN <- Bac[,Bac_CALN]
head(Bac_CALN)
###
Bac_CALN_t <-t(Bac_CALN)
###
library(vegan)
##
dist_Bac_CALN <- vegdist(Bac_CALN_t)
dist_Bac_CALN
###
Bac_CALN_pcoa <- cmdscale(dist_Bac_CALN,eig=TRUE)
Bac_CALN_pc12 <- Bac_CALN_pcoa$points[,1:2]
Bac_CALN_pc12
#write.csv(Bac_CK_pc12,"Bac_CK_pc12")
Bac_CALN_importance <- round(Bac_CALN_pcoa$eig/sum(Bac_CALN_pcoa$eig)*100,digits=2)
###
library(ggplot2)
###
Bac_CALN_pc12 <- as.data.frame(Bac_CALN_pc12)
###
Bac_CALN_pc12$samples <- row.names(Bac_CALN_pc12)
head(Bac_CALN_pc12)
###
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F,row.names=1)
group
group$samples <- row.names(group)#
group
#
Bac_CALN <- c("A4","A13","E19","E28" ,"I25","I27","M8","I19","A11", "A8" ,"E27", "M11","M27","M19","A7")
#
library(tidyverse)
Group_CALN <- group %>% filter(samples==c("A4","A13","E19","E28" ,"I25","I27","M8","I19","A11", "A8" ,"E27", "M11","M27","M19","A7"))
Group_CALN
Group_CALN$Year <- as.factor(Group_CALN$Year)
str(Group_CALN)
#
colnames(Group_CALN) <- c("treat","Approach","Season1","Season2","Sampling time","Year","samples")
head(Group_CALN)
##
library(tidyverse)
df <- merge(Bac_CALN_pc12,Group_CALN,by="samples")
df
##
Season<-factor(df$Season1,levels = c("Dry season","Wet season"))
Year<-factor(df$Year,levels = c("2020","2021"))
##
library(showtext)
windowsFonts()###
#
Font <- c('STKaiti.TTF','simhei.TTF') ##
for (i in Font) {
  font_path = i
  font_name = tools::file_path_sans_ext(basename(font_path))
  font_add(font_name, font_path)
}
font_families() ### 
showtext_auto(enable=TRUE) #

###
Bac_PC_CALN <- ggplot(df,aes(x=V1,y=V2,colour=Season,shape=factor(Year)))+#
  geom_point(size=3)+#
  geom_text(aes(label=samples,y=V2),size=4,vjust=0)#
Bac_PC_CALN

###
Fon <- "serif"
Bac_PC_CALN <- ggplot(df,aes(x=V1,y=V2))+#
  geom_point(size=3,aes(colour=Season,shape=factor(Year)))#+#
#stat_chull(geom="polygon",aes(group=N_approach,color=N_approach,fill=N_approach),alpha=0.1)#
Bac_PC_CALN
##
Bac_PC_CALN <- Bac_PC_CALN + stat_ellipse(linetype = 1,level = 0.95,aes(group  = Season, colour =  Season)) 
Bac_PC_CALN
Bac_PC_CALN <- Bac_PC_CALN+facet_wrap(~treat )#
Bac_PC_CALN
#
Bac_PC_CALN <- Bac_PC_CALN+theme(axis.title.x=element_text(size=12,family="serif"),
                                 axis.title.y=element_text(size=12,family="serif",angle=90),#
                                 axis.text.y=element_text(size=10,family="serif"),#
                                 axis.text.x=element_text(size=10,family="serif"),#
                                 legend.title = element_text(size=10,family="serif"),
                                 legend.text = element_text(size=8,family="serif"),
                                 strip.text = element_text(family = "serif",face='bold',size=13),
                                 strip.background = element_rect(fill=c("#BDBDBD")))
Bac_PC_CALN
Bac_PC_CALN <- Bac_PC_CALN+xlab(paste0("PCo1(",Bac_CALN_importance[1],"%)"))+ylab(paste0("PCo2(",Bac_CALN_importance[2],"%)"))
Bac_PC_CALN

#####
#??adonis2
set.seed(1)
asv.div_O <- adonis2(Bac_CALN_t ~ Season1, data = Group_CALN, permutations = 999, method="bray")
asv.div_O
###
set.seed(1)
asv.div_O1 <- adonis2(Bac_CALN_t ~ Year, data = Group_CALN, permutations = 999, method="bray")
asv.div_O1
####
Fun_O_adonis <- paste0("Season R2: ",round(asv.div_O$R2,2), "; P: ", asv.div_O$`Pr(>F)`)
Fun_O_adonis1 <- paste0("Year R2: ",round(asv.div_O1$R2,2), "; P: ", asv.div_O1$`Pr(>F)`)
# install.packages("ggalt")
#BiocManager::install("ggalt")
library(ggalt)
Bac_PC_CALN <- Bac_PC_CALN+labs(title=Fun_O_adonis,subtitle=Fun_O_adonis1)+#
  theme(plot.title = element_text(family = "serif", #
                                  face = "bold",     #
                                  color = "black",      #
                                  size = 11,          #
                                  hjust = 0,          #
                                  vjust = 0,          #
                                  angle = 0,          #
  ),
  plot.subtitle = element_text(family = "serif", #
                               face = "bold",     #
                               color = "black",      #
                               size = 11,          #
                               hjust = 0,          #
                               vjust = 0,          #
                               angle = 0, )
  )# 
Bac_PC_CALN

###
ggsave("Bac_PC_CALN.pdf",Bac_PC_CALN,width=4,height=3)
ggsave("Bac_PC_CALN.TIFF",Bac_PC_CALN,width =4,height = 3)

#CAHN
Bac_CAHN <- c("E18", "A12", "A19", "I26", "I28", "M12", "M26", "A20", "E12", "E26", "I12","E4","E8","E11","M28")
Bac_CAHN <- Bac[,Bac_CAHN]
head(Bac_CAHN)
###
Bac_CAHN_t <-t(Bac_CAHN)
###
library(vegan)
###
dist_Bac_CAHN <- vegdist(Bac_CAHN_t)
dist_Bac_CAHN
###PCoA
Bac_CAHN_pcoa <- cmdscale(dist_Bac_CAHN,eig=TRUE)
Bac_CAHN_pc12 <- Bac_CAHN_pcoa$points[,1:2]
Bac_CAHN_pc12
#write.csv(Bac_CK_pc12,"Bac_CK_pc12")
Bac_CAHN_importance <- round(Bac_CAHN_pcoa$eig/sum(Bac_CAHN_pcoa$eig)*100,digits=2)
###
library(ggplot2)
###pc12原来是matrix,转化为data.frame
Bac_CAHN_pc12 <- as.data.frame(Bac_CAHN_pc12)
###
Bac_CAHN_pc12$samples <- row.names(Bac_CAHN_pc12)
head(Bac_CAHN_pc12)
###
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F,row.names=1)
group
group$samples <- row.names(group)#
group
#
Bac_CAHN <- c("E18", "A12", "A19", "I26", "I28", "M12", "M26", "A20", "E12", "E26", "I12","E4","E8","E11","M28")
#
library(tidyverse)
Group_CAHN <- group %>% filter(samples==c("E18", "A12", "A19", "I26", "I28", "M12", "M26", "A20", "E12", "E26", "I12","E4","E8","E11","M28"))
Group_CAHN
Group_CAHN$Year <- as.factor(Group_CAHN$Year)
str(Group_CAHN)
###
colnames(Group_CAHN) <- c("treat","Approach","Season1","Season2","Sampling time","Year","samples")
head(Group_CAHN)
###
library(tidyverse)
df <- merge(Bac_CAHN_pc12,Group_CAHN,by="samples")
df
###
Season<-factor(df$Season1,levels = c("Dry season","Wet season"))
Year<-factor(df$Year,levels = c("2020","2021"))
####
library(showtext)
windowsFonts()##
#
Font <- c('STKaiti.TTF','simhei.TTF') ##
for (i in Font) {
  font_path = i
  font_name = tools::file_path_sans_ext(basename(font_path))
  font_add(font_name, font_path)
}
font_families() ###
showtext_auto(enable=TRUE) #

###
Bac_PC_CAHN <- ggplot(df,aes(x=V1,y=V2,colour=Season,shape=factor(Year)))+#
  geom_point(size=3)+#
  geom_text(aes(label=samples,y=V2),size=4,vjust=0)#
Bac_PC_CAHN

###
Fon <- "serif"
Bac_PC_CAHN <- ggplot(df,aes(x=V1,y=V2))+#
  geom_point(size=3,aes(colour=Season,shape=factor(Year)))#+#
#stat_chull(geom="polygon",aes(group=N_approach,color=N_approach,fill=N_approach),alpha=0.1)#
Bac_PC_CAHN
##
Bac_PC_CAHN <- Bac_PC_CAHN + stat_ellipse(linetype = 1,level = 0.95,aes(group  = Season, colour =  Season)) 
Bac_PC_CAHN
Bac_PC_CAHN <- Bac_PC_CAHN+facet_wrap(~treat )#
Bac_PC_CAHN
#
Bac_PC_CAHN <- Bac_PC_CAHN+theme(axis.title.x=element_text(size=12,family="serif"),#
                                 axis.title.y=element_text(size=12,family="serif",angle=90),#
                                 axis.text.y=element_text(size=10,family="serif"),#
                                 axis.text.x=element_text(size=10,family="serif"),#
                                 legend.title = element_text(size=10,family="serif"),
                                 legend.text = element_text(size=8,family="serif"),
                                 strip.text = element_text(family = "serif",face='bold',size=13),
                                 strip.background = element_rect(fill=c("#BDBDBD")))
Bac_PC_CAHN
Bac_PC_CAHN <- Bac_PC_CAHN+xlab(paste0("PCo1(",Bac_CAHN_importance[1],"%)"))+ylab(paste0("PCo2(",Bac_CAHN_importance[2],"%)"))
Bac_PC_CAHN

######
#??adonis2
set.seed(123)
asv.div_O <- adonis2(Bac_CAHN_t ~ Season1, data = Group_CAHN, permutations = 999, method="bray")
asv.div_O
###
set.seed(123)
asv.div_O1 <- adonis2(Bac_CAHN_t ~ factor(Year), data = Group_CAHN, permutations = 999, method="bray")
asv.div_O1
str(Year)
#
####
Fun_O_adonis <- paste0("Season R2: ",round(asv.div_O$R2,2), "; P: ", asv.div_O$`Pr(>F)`)
Fun_O_adonis1 <- paste0("Year R2: ",round(asv.div_O1$R2,2), "; P: ", asv.div_O1$`Pr(>F)`)
# install.packages("ggalt")
#BiocManager::install("ggalt")
library(ggalt)
Bac_PC_CAHN <- Bac_PC_CAHN+labs(title=Fun_O_adonis,subtitle=Fun_O_adonis1)+#
  theme(plot.title = element_text(family = "serif", #
                                  face = "bold",     #
                                  color = "black",      #
                                  size = 11,          #
                                  hjust = 0,          #
                                  vjust = 0,          #
                                  angle = 0,          #
  ),
  plot.subtitle = element_text(family = "serif", 
                               face = "bold",     #
                               color = "black",      #
                               size = 11,          #
                               hjust = 0,          #
                               vjust = 0,          #
                               angle = 0, )
  )# 
Bac_PC_CAHN

###
ggsave("Bac_PC_CAHN.pdf",Bac_PC_CAHN,width=4,height=3)
ggsave("Bac_PC_CAHN.TIFF",Bac_PC_CAHN,width =4,height = 3)

#UALN
Bac_UALN <- c("A18" ,"A6",  "A14", "E20", "I8"  ,"I15" ,"I4"  ,"E22", "E25" ,"E6","E13" ,"I6" , "I13" ,"I11" ,"I7")
Bac_UALN <- c("A18" ,"A6",  "A14", "E20", "I8"  ,"I15" ,"I4"  ,"E22", "E25" ,"E6","E13" ,"I6" , "I13" ,"I11" ,"I7")
Bac_UALN <- Bac[,Bac_UALN]
head(Bac_UALN)
###
Bac_UALN_t <-t(Bac_UALN)
###
library(vegan)
###
dist_Bac_UALN <- vegdist(Bac_UALN_t)
dist_Bac_UALN
###
Bac_UALN_pcoa <- cmdscale(dist_Bac_UALN,eig=TRUE)
Bac_UALN_pc12 <- Bac_UALN_pcoa$points[,1:2]
Bac_UALN_pc12
#write.csv(Bac_CK_pc12,"Bac_CK_pc12")
Bac_UALN_importance <- round(Bac_UALN_pcoa$eig/sum(Bac_UALN_pcoa$eig)*100,digits=2)
###
library(ggplot2)
###pc12
Bac_UALN_pc12 <- as.data.frame(Bac_UALN_pc12)
###
Bac_UALN_pc12$samples <- row.names(Bac_UALN_pc12)
head(Bac_UALN_pc12)
###
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F,row.names=1)
group
group$samples <- row.names(group)#
group
#
Bac_UALN <- c("A18" ,"A6",  "A14", "E20", "I8"  ,"I15" ,"I4"  ,"E22", "E25" ,"E6","E13" ,"I6" , "I13" ,"I11" ,"I7")
#
library(tidyverse)
Group_UALN <- group %>% filter(samples==c("A18" ,"A6",  "A14", "E20", "I8"  ,"I15" ,"I4"  ,"E22", "E25" ,"E6","E13" ,"I6" , "I13" ,"I11" ,"I7"))
Group_UALN
Group_UALN$Year <- as.factor(Group_UALN$Year)
str(Group_UALN)
###
colnames(Group_UALN) <- c("treat","Approach","Season1","Season2","Sampling time","Year","samples")
head(Group_UALN)
###
library(tidyverse)
df <- merge(Bac_UALN_pc12,Group_UALN,by="samples")
df
###
Season<-factor(df$Season1,levels = c("Dry season","Wet season"))
Year<-factor(df$Year,levels = c("2020","2021"))
###
library(showtext)
windowsFonts()###
#
Font <- c('STKaiti.TTF','simhei.TTF') ##
for (i in Font) {
  font_path = i
  font_name = tools::file_path_sans_ext(basename(font_path))
  font_add(font_name, font_path)
}
font_families() ### 
showtext_auto(enable=TRUE) #

###
Bac_PC_UALN <- ggplot(df,aes(x=V1,y=V2,colour=Season,shape=factor(Year)))+#
  geom_point(size=3)+#
  geom_text(aes(label=samples,y=V2),size=4,vjust=0)#
Bac_PC_UALN

###
Fon <- "serif"
Bac_PC_UALN <- ggplot(df,aes(x=V1,y=V2))+
  geom_point(size=3,aes(colour=Season,shape=factor(Year)))#+
#stat_chull(geom="polygon",aes(group=N_approach,color=N_approach,fill=N_approach),alpha=0.1)#
Bac_PC_UALN
##
Bac_PC_UALN <- Bac_PC_UALN + stat_ellipse(linetype = 1,level = 0.95,aes(group  = Season, colour =  Season)) 
Bac_PC_UALN
Bac_PC_UALN <- Bac_PC_UALN+facet_wrap(~treat )#
Bac_PC_UALN
#
Bac_PC_UALN <- Bac_PC_UALN+theme(axis.title.x=element_text(size=12,family="serif"),#
                                 axis.title.y=element_text(size=12,family="serif",angle=90),#
                                 axis.text.y=element_text(size=10,family="serif"),#
                                 axis.text.x=element_text(size=10,family="serif"),#
                                 legend.title = element_text(size=10,family="serif"),
                                 legend.text = element_text(size=8,family="serif"),
                                 strip.text = element_text(family = "serif",face='bold',size=13),
                                 strip.background = element_rect(fill=c("#BDBDBD")))
Bac_PC_UALN
Bac_PC_UALN <- Bac_PC_UALN+xlab(paste0("PCo1(",Bac_UALN_importance[1],"%)"))+ylab(paste0("PCo2(",Bac_UALN_importance[2],"%)"))
Bac_PC_UALN

######
#??adonis2
set.seed(1)
asv.div_O <- adonis2(Bac_UALN_t ~ Season1, data = Group_UALN, permutations = 999, method="bray")
asv.div_O
###
set.seed(1)
asv.div_O1 <- adonis2(Bac_UALN_t ~ Year, data = Group_UALN, permutations = 999, method="bray")
asv.div_O1

#
####
Fun_O_adonis <- paste0("Season R2: ",round(asv.div_O$R2,2), "; P: ", asv.div_O$`Pr(>F)`)
Fun_O_adonis1 <- paste0("Year R2: ",round(asv.div_O1$R2,2), "; P: ", asv.div_O1$`Pr(>F)`)
# install.packages("ggalt")
#BiocManager::install("ggalt")
library(ggalt)
Bac_PC_UALN <- Bac_PC_UALN+labs(title=Fun_O_adonis,subtitle=Fun_O_adonis1)+#
  theme(plot.title = element_text(family = "serif", #
                                  face = "bold",     #
                                  color = "black",      #
                                  size = 11,          #
                                  hjust = 0,          #
                                  vjust = 0,          #
                                  angle = 0,          #
  ),
  plot.subtitle = element_text(family = "serif", #
                               face = "bold",     #
                               color = "black",      #
                               size = 11,          #
                               hjust = 0,          #
                               vjust = 0,          #
                               angle = 0, )
  )# 
Bac_PC_UALN

###
ggsave("Bac_PC_UALN.pdf",Bac_PC_UALN,width=4,height=3)
ggsave("Bac_PC_UALN.TIFF",Bac_PC_UALN,width =4,height = 3)

#UAHN
Bac_UAHN <- c("A15" ,"A21", "E15", "I22", "I18" ,"I14", "I21", "E5",  "E7","E14", "E21" ,"M7"  ,"M14" ,"M21", "I20")
Bac_UAHN <- Bac[,Bac_UAHN]
head(Bac_UAHN)
###
Bac_UAHN_t <-t(Bac_UAHN)
###
library(vegan)
###
dist_Bac_UAHN <- vegdist(Bac_UAHN_t)
dist_Bac_UAHN
###
Bac_UAHN_pcoa <- cmdscale(dist_Bac_UAHN,eig=TRUE)
Bac_UAHN_pc12 <- Bac_UAHN_pcoa$points[,1:2]
Bac_UAHN_pc12
#write.csv(Bac_CK_pc12,"Bac_CK_pc12")
Bac_UAHN_importance <- round(Bac_UAHN_pcoa$eig/sum(Bac_UAHN_pcoa$eig)*100,digits=2)
###
library(ggplot2)
###
Bac_UAHN_pc12 <- as.data.frame(Bac_UAHN_pc12)
###
Bac_UAHN_pc12$samples <- row.names(Bac_UAHN_pc12)
head(Bac_UAHN_pc12)
###
group <- read.csv("jgsgrouptotalBac(new).csv",header = T,check.names=F,row.names=1)
group
group$samples <- row.names(group)
group
#
Bac_UAHN <- c("A15" ,"A21", "E15", "I22", "I18" ,"I14", "I21", "E5",  "E7","E14", "E21" ,"M7"  ,"M14" ,"M21", "I20")
#
library(tidyverse)
Group_UAHN <- group %>% filter(samples==c("A15" ,"A21", "E15", "I22", "I18" ,"I14", "I21", "E5",  "E7","E14", "E21" ,"M7"  ,"M14" ,"M21", "I20"))
Group_UAHN
Group_UAHN$Year <- as.factor(Group_UAHN$Year)
str(Group_UAHN)
###
colnames(Group_UAHN) <- c("treat","Approach","Season1","Season2","Sampling time","Year","samples")
head(Group_UAHN)
###
library(tidyverse)
df <- merge(Bac_UAHN_pc12,Group_UAHN,by="samples")
df
###
Season<-factor(df$Season1,levels = c("Dry season","Wet season"))
Year<-factor(df$Year,levels = c("2020","2021"))
####
library(showtext)
windowsFonts()###
#
Font <- c('STKaiti.TTF','simhei.TTF') ##
for (i in Font) {
  font_path = i
  font_name = tools::file_path_sans_ext(basename(font_path))
  font_add(font_name, font_path)
}
font_families() ### 
showtext_auto(enable=TRUE) #

###
Bac_PC_UAHN <- ggplot(df,aes(x=V1,y=V2,colour=Season,shape=factor(Year)))+#
  geom_point(size=3)+#
  geom_text(aes(label=samples,y=V2),size=4,vjust=0)#
Bac_PC_UAHN

###
Fon <- "serif"
Bac_PC_UAHN <- ggplot(df,aes(x=V1,y=V2))+
  geom_point(size=3,aes(colour=Season,shape=factor(Year)))#+#
#stat_chull(geom="polygon",aes(group=N_approach,color=N_approach,fill=N_approach),alpha=0.1)#
Bac_PC_UAHN
##
Bac_PC_UAHN <- Bac_PC_UAHN + stat_ellipse(linetype = 1,level = 0.95,aes(group  = Season, colour =  Season)) 
Bac_PC_UAHN
Bac_PC_UAHN <- Bac_PC_UAHN+facet_wrap(~treat )#
Bac_PC_UAHN
#
Bac_PC_UAHN <- Bac_PC_UAHN+theme(axis.title.x=element_text(size=12,family="serif"),#
                                 axis.title.y=element_text(size=12,family="serif",angle=90),#
                                 axis.text.y=element_text(size=10,family="serif"),#
                                 axis.text.x=element_text(size=10,family="serif"),#
                                 legend.title = element_text(size=10,family="serif"),
                                 legend.text = element_text(size=8,family="serif"),
                                 strip.text = element_text(family = "serif",face='bold',size=13),
                                 strip.background = element_rect(fill=c("#BDBDBD")))
Bac_PC_UAHN
Bac_PC_UAHN <- Bac_PC_UAHN+xlab(paste0("PCo1(",Bac_UAHN_importance[1],"%)"))+ylab(paste0("PCo2(",Bac_UAHN_importance[2],"%)"))
Bac_PC_UAHN

######
#??adonis2
set.seed(123)
asv.div_O <- adonis2(Bac_UAHN_t ~ Season1, data = Group_UAHN, permutations = 999, method="bray")
asv.div_O
##
set.seed(123)
asv.div_O1 <- adonis2(Bac_UAHN_t ~ Year, data = Group_UAHN, permutations = 999, method="bray")
asv.div_O1

#
####
Fun_O_adonis <- paste0("Season R2: ",round(asv.div_O$R2,2), "; P: ", asv.div_O$`Pr(>F)`)
Fun_O_adonis1 <- paste0("Year R2: ",round(asv.div_O1$R2,2), "; P: ", asv.div_O1$`Pr(>F)`)
# install.packages("ggalt")
#BiocManager::install("ggalt")
library(ggalt)
Bac_PC_UAHN <- Bac_PC_UAHN+labs(title=Fun_O_adonis,subtitle=Fun_O_adonis1)+#
  theme(plot.title = element_text(family = "serif", #
                                  face = "bold",     #
                                  color = "black",      
                                  size = 11,          #
                                  hjust = 0,          #
                                  vjust = 0,          #
                                  angle = 0,          #
  ),
  plot.subtitle = element_text(family = "serif", #
                               face = "bold",     #
                               color = "black",      #
                               size = 11,          #
                               hjust = 0,          #
                               vjust = 0,          #
                               angle = 0, )
  )# 
Bac_PC_UAHN

###
ggsave("Bac_PC_UAHN.pdf",Bac_PC_UAHN,width=4,height=3)
ggsave("Bac_PC_UAHN.TIFF",Bac_PC_UAHN,width =4,height = 3)

######JGSbacLMM
setwd("F:/JGS-SMT微生物测序数据/R分析/20230119Bacteria")
library(tidyverse)
library(lme4)
library(tidyr)
library(Matrix)
library(nlme)
library(lmerTest)
library(glmm.hp)
data<-read.csv("BacAlphLMM.csv",header=T,row.names=1)
str(data)
shapiro.test(data$Shannon)
shapiro.test(data$Simpson)#需转换
shapiro.test(data$Richness)
shapiro.test(data$Chao)
#
Shannon<-lmer(Shannon~factor(Treatment)*factor(Year)*factor(Season)+(1|Block), data=data)
summary(Shannon)
anova(Shannon)
ranova(Shannon)
glmm.hp(Shannon)

Simpson<-lmer(Simpson~factor(Treatment)*factor(Time)+(1|Block), data=data2)
summary(Simpson)
anova(Simpson)
ranova(Simpson)
#isSingular(Simpson,tol=1e-05)
glmm.hp(Simpson)

Richness<-lmer(Richness~factor(Treatment)*factor(Year)*factor(Season)+(1|Block), data=data)
summary(Richness)
anova(Richness)
ranova(Richness)

Chao<-lmer(Chao~factor(Treatment)*factor(Time)+(1|Block), data=data2)
summary(Chao)
anova(Chao)
ranova(Chao)



######SMF
######SMF
######SMF
######SMF
######SMF
library(reshape2)
library(ggplot2)
library(vegan)
library(picante)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
setwd("F:/JGS-SMT微生物测序数据/功能基因")
df <- read.csv("JGS_Multi.csv",header = T,row.names=1, sep= ",")
df$Group <- as.factor(df$Group)
str(df)
df$Group<-factor(df$Group,levels = c("CK","CALN","CAHN","UALN","UAHN"))##
#
shapiro.test(df$LC)
df %>% select(Group,LC) %>% group_by(Group) %>%
  summarise_all(.funs=funs(statistic=shapiro.test(.)$statistic,
                           p.value=shapiro.test(.)$p.value))###
#
bartlett.test(Sucrase~Group,data=df)
#
display.brewer.all()##
mycol <- brewer.pal(n=5,name='Set1')##
mycol
#DOC
Label<-c(rep("DOC",75))
A <- data.frame(Label)
df<-mutate(df,A)
df
DOC <- ggboxplot(df,x="Group",y="DOC",color="Group",
                 palette=mycol,legend="NA",#"bottom",
                 add="jitter",#
                 add.params=list(size=1.2,jitter=0.1,alpha=0.5),outlaier.shape=NA,#jitter=0.1
                 facet.by="Label")+
  stat_compare_means(method = "kruskal.test",label.y = 4,label.x = 1,
                     label = 'p.format',family = "serif")+
  stat_compare_means(comparisons = list(c('CK','CALN'),
                                        c('CK','CAHN'),
                                        c('CK','UALN'),
                                        c('CK','UAHN'),
                                        c('CALN','CAHN'),
                                        c('CALN','UALN'),
                                        c('CALN','UAHN'),
                                        c('CAHN','UALN'),
                                        c('CAHN','UAHN'),
                                        c('CAHN','UAHN')),
                     method='wilcox',label="p.signif")+##
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color=NA),
        panel.border = element_rect(fill = NA,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black",face = "bold",hjust = 1,family = "serif",angle =45),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10,colour = "black",face = "bold",hjust = 0.5,family = "serif"),
        axis.title.y = element_text(vjust=0.2,size = 12,face = "bold",family = "serif"),
        strip.text = element_text(size = 14,family = "serif",face='bold'),
        strip.background = element_rect(fill=c("#BDBDBD")))#
DOC
DOC <- DOC+labs(x=NULL,y=bquote(~g ~kg^-1))
DOC
#C_1 <- C_1+scale_y_continuous(limits=c(1,9))
#C_1
DOC <- DOC + theme(legend.position="right",
                   legend.title=element_text(size=14,face="bold",colour="black",family="serif"),
                   legend.text = element_text(family = "serif",
                                              colour = "black", 
                                              size = 14))
DOC

#####

######
######
######
######
??multifunc
library(multifunc)
library(tidyverse)
library(gridExtra)
library(plyr)
#for plotting
library(ggplot2)
library(patchwork)

#for data
library(tidyr)
library(dplyr)
library(purrr)
library(forcats)
#for analysus
library(car)
setwd("F:/JGS-SMT微生物测序数据/功能基因")
df <- read.csv("JGS_Multi.csv",header = T,row.names=1, sep= ",")%>% as_tibble()

allVars<-qw(NH4,NO3,TC,TN,SCN,DOC,DON,DOCN,AP,TP,SNP,LC,LN,LCN, ACP,Catalase, BG, Sucrase, Urease,
            C_degra, C_fixa, Me_Meta, Denitri, Nitri, P_cycling, S_cycling)
allVars
varIdx<-which(names(df) %in% allVars)
varIdx
Soil_vars <- df[,allVars]
#Soil_vars <- scale(Soil_vars)
#str(Soil_vars)
#Soil_vars <- as.matrix(Soil_vars)
#计算方差膨胀因子
library(tidyverse)
library(caret)
#scale
dim(Soil_vars)
#
#model1 <- lm(B_Shannon~NH4+NO3+TC+TN+SCN+DOC+DON+DOCN+AP+TP+SNP+LC+LN+LCN+ ACP+ Catalase+ BG+ Sucrase+ Urease+
# C_degra+ C_fixa+ Me_Meta+ Denitri+ Nitri+ P_cycling+ S_cycling,data=df)
#summary(model1)
#car::vif(model1)
model1 <- lm(scale(B_Shannon)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TC)+scale(TN)+scale(SCN)+scale(DOC)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LC)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model1)
car::vif(model1)

#去
model1 <- lm(scale(B_Shannon)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TC)+scale(TN)+scale(DOC)+scale(DON)+scale(AP)+scale(TP)+
               scale(LC)+scale(LN)+
               scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+scale(Urease)+
               scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model1)
car::vif(model1)
#
model2 <- lm(scale(B_Richness)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TC)+scale(TN)+scale(SCN)+scale(DOC)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LC)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model2)
car::vif(model2)
#
model2 <- lm(scale(B_Richness)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TN)+scale(SCN)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model2)
car::vif(model2)
#
model3 <- lm(scale(F_Shannon)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TC)+scale(TN)+scale(SCN)+scale(DOC)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LC)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model3)
car::vif(model3)
#
#去
model3 <- lm(scale(F_Shannon)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TN)+scale(SCN)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model3)
car::vif(model3)
#
model4 <- lm(scale(F_Richness)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TC)+scale(TN)+scale(SCN)+scale(DOC)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LC)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model4)
car::vif(model4)
#
model4 <- lm(scale(F_Richness)~scale(SWC)+scale(pH)+scale(NH4)+scale(NO3)+scale(TN)+scale(SCN)+scale(DON)+scale(DOCN)+scale(AP)+
               scale(TP)+scale(SNP)+scale(LN)+scale(LCN)+scale(ACP)+ scale(Catalase)+ scale(BG)+ scale(Sucrase)+
               scale(Urease)+scale(C_degra)+ scale(C_fixa)+ scale(Me_Meta)+ scale(Denitri)+ scale(Nitri)+ scale(P_cycling)+ scale(S_cycling),data=df)
summary(model4)
car::vif(model4)


##
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(ggfun)##
library(ggpubr)
display.brewer.all()#
mycol <- brewer.pal(n=5,name='Set1')#
p1<-ggplot(df, aes(x=B_Richness,y=ACP)) + 
  geom_point(shape=21,color="black",fill="#36BED9",stroke=0.25,alpha=0.8)+
  scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"))+
  scale_size(range = c(1, 8))+
  geom_smooth(formula=y~x,method="lm",alpha=0.5,fill="grey",color="#001eef")
p1
#####
##
###
Aver_multi<-cbind(df, getStdAndMeanFunctions(df, allVars,standardizeZScore))
write.csv(Aver_multi,'Aver_multi.csv')
Aver_multi <- read.csv("Aver_multi01.csv",header = T,row.names=1)
#plot it
library(ggpmisc)
A1 <- ggplot(aes(x=scale(B_Shannon), y=scale(meanFunction)),data=Aver_multi)+geom_point(size=3,alpha=0.6)+
  theme_bw(base_size=15)+
  theme(axis.title = element_text(
    family = "serif",#
    face='bold', 
    size=1, #
    lineheight = 1),#
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=12),
    strip.text = element_text(family = "serif",face='bold'))+ #
  stat_smooth(method="lm", colour="black", size=2) +
  stat_fit_deviations(formula = y ~ x, color = NA) +
  stat_poly_eq(use_label(c("adj.R2", "p.value.label")),family = "serif")+
  facet_wrap(~group1 )+
  xlab("")+
  ylab("")
A1
#xlab("\nB_Shannon (z score)") +
#ylab("Average multifunctionality (z score)\n")


model <- lm(meanFunction~scale(B_Shannon),data=Aver_multi)
values <- summary(model)
values


A2 <- ggplot(aes(x=scale(B_Richness), y=scale(meanFunction)),data=Aver_multi)+geom_point(size=3,alpha=0.6)+
  theme_bw(base_size=15)+
  theme(axis.title = element_text(
    family = "serif",##
    face='bold', #
    size=1, ##
    lineheight = 1),#
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=12),
    strip.text = element_text(family = "serif",face='bold'))+ ##
  stat_smooth(method="lm", colour="black", size=2) +
  stat_fit_deviations(formula = y ~ x, color = NA) +
  stat_poly_eq(use_label(c("adj.R2", "p.value.label")),family = "serif")+
  facet_wrap(~group2 )+
  xlab("")+
  ylab("")#+
#scale_y_continuous(name=NULL,labels=NULL)
#xlab("\nB_Richness (z score)") +
#ylab("Average multifunctionality (z score)\n")

A2
model <- lm(meanFunction~scale(B_Shannon),data=Aver_multi)
values <- summary(model)
values

B1 <- ggplot(aes(x=scale(F_Shannon), y=scale(meanFunction)),data=Aver_multi)+geom_point(size=3,alpha=0.6)+
  theme_bw(base_size=15)+
  theme(axis.title = element_text(
    family = "serif",#
    face='bold', #
    size=1, ##
    lineheight = 1),
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=12),
    strip.text = element_text(family = "serif",face='bold'))+ ##
  stat_smooth(method="lm", colour="black", size=2) +
  stat_fit_deviations(formula = y ~ x, color = NA) +
  stat_poly_eq(use_label(c("adj.R2", "p.value.label")),family = "serif")+
  facet_wrap(~group1 )+
  scale_y_continuous(limits=c(-2.5,3))+
  #abels=c(0.4,0.5,0.6,0.7))+
  xlab("")+
  ylab("")
#xlab("\nF_Shannon (z score)") +
#ylab("Average multifunctionality (z score)\n")

B1

#####
#####
#####
#####
setwd("F:/JGS-SMT微生物测序数据/功能基因")
library(tidyverse)
library(lme4)
library(tidyr)
library(Matrix)
library(nlme)
library(lmerTest)
library(glmm.hp)
Aver_multi_LMM01 <- read.csv("Aver_multi_CAN.csv",header = T,row.names=1)
Aver_lmm01 <-lmer(meanFunction~factor(group)*factor(Time)+(1|Block), data=Aver_multi_LMM01)
summary(Aver_lmm01)
anova(Aver_lmm01)
ranova(Aver_lmm01)
glmm.hp(Aver_lmm01)


##
##
Aver_multi <- read.csv("Aver_multi_CAN.csv",header = T,row.names=1)
#plot it
#
library(RColorBrewer)
display.brewer.all()##
mycol <- brewer.pal(n=9,name='Set1')#
mycol
color=c("#E41A1C","#377EB8","#4DAF4A")
Aver_multi$group<-factor(Aver_multi$group,levels = c("CK","CALN","CAHN"))##
####F_Shannon
library(ggpmisc)
BS <- ggplot(Aver_multi,aes(x=scale(F_Shannon), y=meanFunction,colour = group,shape=group))+
  scale_color_manual(values = color)+
  geom_point(size=2,alpha=0.8)+#geom_point(size=3,alpha=0.6)+
  geom_smooth(method="lm",formula = y ~ x,se=T,linewidth=0.8,alpha=0.3) +#colour="black"
  #stat_fit_deviations(formula = y ~ x)+
  stat_poly_eq(use_label(c("adj.R2", "p.value.label"),sep = "*\", \"*"),family = "serif",size=5,
               geom = "text", label.x = -2, label.y = c(0.75,0.7,0.65), hjust = 0,position = "identity")+
  theme_bw(base_size=15)+
  theme(axis.title = element_text(
    family = "serif",##
    face='bold', ##
    size=0, ##
    lineheight = 0),#
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=12),
    strip.text = element_text(family = "serif",face='bold',size=15))+
  facet_wrap(~group1 )+
  labs(x=NULL,y=NULL)
BS <- BS + theme(legend.position="right",
                 legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                 legend.text = element_text(family = "serif",
                                            colour = "black", 
                                            size = 10))
BS
#xlab("\nB_Shannon (z score)") +
#ylab("Average multifunctionality (z score)\n")
#scale_color_manual(values = color)+
#scale_fill_manual(values = color)

??use_label
model <- lm(meanFunction~scale(B_Shannon),data=Aver_multi)
values <- summary(model)
values
###
##F_Richness##
library(ggpmisc)
BR <- ggplot(Aver_multi,aes(x=scale(F_Richness), y=meanFunction,colour = group,shape=group))+
  geom_point(size=2,alpha=0.8)+#geom_point(size=3,alpha=0.6)+
  scale_color_manual(values = color)+
  geom_smooth(method="lm" , linewidth=0.8,alpha=0.3) +#colour="black"
  #stat_fit_deviations(formula = y ~ x)+
  stat_poly_eq(use_label(c("adj.R2", "p.value.label"),sep = "*\", \"*"),family = "serif",size=5,
               geom = "text", label.x = -2, label.y = c(0.75,0.7,0.65), hjust = 0,position = "identity")+
  theme_bw(base_size=15)+
  theme(axis.title = element_text(
    family = "serif",##
    face='bold', ##
    size=0, #
    lineheight = 0),##
    axis.text = element_text(
      family = "serif",
      color="black",
      size=12),
    strip.text = element_text(family = "serif",face='bold',size=15))+
  facet_wrap(~group2 )+
  labs(x=NULL,y=NULL)
BR <- BR + theme(legend.position="right",
                 legend.title=element_text(size=12,face="bold",colour="black",family="serif"),
                 legend.text = element_text(family = "serif",
                                            colour = "black", 
                                            size = 10))
BR
#xlab("\nB_Shannon (z score)") +
#ylab("Average multifunctionality (z score)\n")
BS<-BS+theme(plot.margin = margin(0,1.5,0,1)) # 
BR<-BR+theme(plot.margin = margin(0,1.5,0,1))
library(patchwork)
A7 <- BS+BR+plot_layout(ncol=2,nrow=1,guides = "collect")&theme(legend.position = "top",
                                                                legend.title=element_text(size=15,face="bold",colour="black",family="serif"),
                                                                legend.text = element_text(family = "serif",colour = "black", size = 14))



A7
ggsave('J_Aver_multi_Appr_CAN.pdf',A7,width=9,height=3.7)
ggsave("J_Aver_multi_Appr_CAN.TIFF",A7,width = 9,height = 3.7)


##
######
library(multifunc)
library(tidyverse)
library(gridExtra)
library(plyr)
#for plotting
library(ggplot2)
library(patchwork)

#for data
library(tidyr)
library(dplyr)
library(purrr)
library(forcats)
#for analysus
library(car)
#####
#####
#####
setwd("F:/JGS-SMT微生物测序数据/功能基因")
df <- read.csv("JGS_Multi_CAN.csv",header = T,row.names=1, sep= ",")%>% as_tibble()

allVars<-qw(NH4,NO3,TC,TN,SCN,DOC,DON,DOCN,AP,TP,SNP,LC,LN,LCN, ACP, Catalase, BG, Sucrase, Urease,
            C_degra, C_fixa, Me_Meta, Denitri, Nitri, P_cycling, S_cycling)
allVars
###B—Shannon##
Thresh_BS<-getFuncsMaxed(df, allVars, threshmin=0.05, threshmax=0.99, prepend=c("B_Shannon"), maxN=7)
allVarsLinear_BS<-glm(funcMaxed ~ B_Shannon, data=subset(Thresh_BS, Thresh_BS$thresholds=="0.8"), family=quasipoisson(link="identity"))
Anova(allVarsLinear_BS, test.statistic="F")
summary(allVarsLinear_BS)
Thresh_BS$percent <- 100*Thresh_BS$thresholds

#####20%,40%,60%,80%
Plot_BS<-subset(Thresh_BS, Thresh_BS$thresholds %in% qw(0.2, 0.4, 0.6, 0.8))
Plot_BS$percent<-paste(100*Plot_BS$thresholds, "%", sep="")
#绘制阈值图
library(ggpmisc)
B_S_T <- ggplot(Plot_BS, aes(x=B_Shannon, y=funcMaxed))+geom_point(size=1.5,alpha=0.6,colour="#377EB8")+
  stat_smooth(method="glm", 
              method.args = list(family=quasipoisson(link="identity")),
              colour="#377EB8",fill="#377EB8", linewidth=1,alpha=0.3) +
  ylab(expression("Number of Functions " >=Threshold)) +
  xlab("Bacterial shannon index") +
  facet_wrap(~percent,scales = 'free_y' )+
  stat_poly_eq(use_label(c("adj.R2","p.value.label"),sep = "*\", \"*"),family = "serif",size=3.5,colour="#377EB8")+
  theme_bw(base_size=12) +
  theme(axis.title = element_text(
    family = "serif",#
    size=12, #
    lineheight = 0),#
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=12),strip.text = element_text(family = "serif",face='bold',size=12))+
  #scale_y_continuous(limits=c(0,30))+
  theme(strip.background = element_rect(fill=c("#BDBDBD")))

B_S_T
ggsave('B_S_T_CAN_1.pdf',B_S_T,width=8,height=6)
ggsave("B_S_T_CAN_1.TIFF",B_S_T,width = 5,height = 4)
####
BS_LinearSlopes<-getCoefTab(funcMaxed ~ B_Shannon,
                            data = Thresh_BS, 
                            coefVar = "B_Shannon",
                            family = quasipoisson(link="identity"))

BS_LinearSlopes

BS_IDX <- getIndices(BS_LinearSlopes, Thresh_BS, 
                     funcMaxed ~ B_Shannon)
BS_IDX

BS_LinearSlopes$estimate[which(BS_LinearSlopes$thresholds == BS_IDX$Tmde)]
Thresh_BS$IDX <- 0
Thresh_BS$IDX [which(Thresh_BS$thresholds %in%
                       c(BS_IDX$Tmin, BS_IDX$Tmax,
                         BS_IDX$Tmde))] <- 1
###
ggplot(data=Thresh_BS, aes(x=B_Shannon, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions " >= Threshold)) +
  xlab("B_Shannon") +
  stat_smooth(method="glm", family=quasipoisson(link="identity"), lwd=0.8,
              fill=NA, aes(color=percent)) +
  theme_bw(base_size=14) +
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")
###
#
#
Label<-as.data.frame(rep("Shannon (Bacteria)",4275))
colnames(Label) <- c("SB")
Label
#
Thresh_BS <- cbind(Thresh_BS,Label)
Thresh_BS
#
BS_Thres <- ggplot(data=Thresh_BS, 
                   aes(x=B_Shannon, y=funcMaxed, group=percent)) +
  geom_smooth(method="glm", se = F,
              aes(color=percent, lwd = IDX)) +
  ylab(expression("Number of Functions " >=Threshold)) +
  xlab("Index (unitless)") +
  scale_color_gradient(name="Percent of \nMaximum",
                       low="blue", high="red") +
  scale_size(range=c(0.3,2), guide="none") +
  #annotate(geom="text",size=3.5,face='bold', family = "serif",x=5.57, y=c(14.9,16.5,25.5),
  #label=c("Tmax", "Tmde", "Tmin"))+
  #annotate(geom="text",size=3.5, x=7.45, family = "serif",
  #y=c(26,22,20),#y=c(BS_IDX$Mmin, BS_IDX$Mmax,BS_IDX$Mmde)
  #label=c("Mmin", "Mmde","Mmax" ))+
  facet_wrap(~SB )+#
  theme_bw(base_size=10) +
  theme(axis.title = element_text(
    family = "serif",#
    size=10, #
    lineheight = 0),#
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=10),strip.text = element_text(family = "serif",face='bold',size=12))+
  scale_y_continuous(limits=c(-5,27))+
  scale_x_continuous(limits=c(5.7,7.3))+
  theme(#legend.position="top",
    axis.title.x=element_text(size=15,family="serif",lineheight = 0),
    axis.title.y=element_text(size=15,family="serif",margin = margin(r=0.1,unit='cm'),lineheight = 1),#
    axis.text.y=element_text(size=13,family="serif",color="black"),#
    axis.text.x=element_text(size=13,family="serif",color="black"),#
    legend.title = element_text(size=13,family="serif"),
    legend.text = element_text(size=11,family="serif"),
    strip.text = element_text(family = "serif",face='bold',size=13),
    strip.background = element_rect(fill=c("#BDBDBD")))

BS_Thres

ggsave('BS_Thres_CAN_1.pdf',BS_Thres,width=4,height=5)
ggsave("BS_Thres_CAN_1.TIFF",BS_Thres,width = 4,height = 5)
#####Change in Number of Functions per Addition of 1 Species
#####Change in Number of Functions per Addition of 1 Species
#
library(RColorBrewer)
display.brewer.all()##
mycol <- brewer.pal(n=9,name='Blues')##
mycol
color=c("#377EB8","#984EA3")
#
#
Label<-as.data.frame(rep("Shannon (Bacteria)",95))
colnames(Label) <- c("SB")
Label
#
BS_LinearSlopes <- cbind(BS_LinearSlopes,Label)
BS_LinearSlopes
#
BS_CNF <- ggplot(BS_LinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="#6BAED6",alpha=0.5,aes(x=thresholds*100, 
                                           ymin = estimate - 1.96*.data[["std.error"]],
                                           ymax= estimate + 1.96*.data[["std.error"]])) +
  geom_point(size=0.8,aes(x=thresholds*100, y = estimate)) +
  labs(x = "Threshold (%)", 
       y ="Change in Number of Functions\nper Addition of Species") +
  geom_hline(yintercept = 0,linetype=2)+
  #annotate(geom="text",family = "serif",size=4,face='bold', y=c(-0.3, -0.3, -0.3, 
  #BS_IDX$Rmde.linear+0.6),
  #x=c(BS_IDX$Tmin*100, 
  #BS_IDX$Tmde*65,
  #BS_IDX$Tmax*120, 
  #BS_IDX$Tmde*100),
  #label=c("Tmin", "Tmde", "Tmax","Rmde"), 
  #color="black")+
  theme_bw(base_size=12)+
  facet_wrap(~SB )#
BS_CNF
BS_CNF <- BS_CNF+theme(axis.title.x=element_text(size=15,family="serif",lineheight = 0),#
                       axis.title.y=element_text(size=15,family="serif",margin = margin(r=0.1,unit='cm'),lineheight = 1),#
                       axis.text.y=element_text(size=13,family="serif",color="black"),#
                       axis.text.x=element_text(size=13,family="serif",color="black"),#
                       legend.title = element_text(size=13,family="serif"),
                       legend.text = element_text(size=11,family="serif"),
                       strip.text = element_text(family = "serif",face='bold',size=13),
                       strip.background = element_rect(fill=c("#BDBDBD")))
BS_CNF
BS_CNF <- BS_CNF+
  scale_x_continuous(limits=c(0,100))
BS_CNF
library(patchwork)
BS_Thres<-BS_Thres+theme(plot.margin = margin(0,1.5,0,1)) # 
BS_CNF<-BS_CNF+theme(plot.margin = margin(0,1.5,0,1))

BS_Mul_CAN <- BS_Thres+BS_CNF+plot_layout(ncol=1,nrow=2)
BS_Mul_CAN
ggsave('BS_Mul_CAN_1.pdf',BS_Mul_CAN,width=4,height=8)
ggsave("BS_Mul_CAN_1.TIFF",BS_Mul_CAN,width = 4,height = 8)

###B—Richness####
Thresh_BR<-getFuncsMaxed(df, allVars, threshmin=0.05, threshmax=0.99, prepend=c("B_Richness"), maxN=7)
allVarsLinear_BR<-glm(funcMaxed ~ B_Richness, data=subset(Thresh_BR, Thresh_BR$thresholds=="0.8"), family=quasipoisson(link="identity"))
Anova(allVarsLinear_BR, test.statistic="F")
summary(allVarsLinear_BR)
Thresh_BR$percent <- 100*Thresh_BR$thresholds

#####20%,40%,60%,80%
Plot_BR<-subset(Thresh_BR, Thresh_BR$thresholds %in% qw(0.2, 0.4, 0.6, 0.8))
Plot_BR$percent<-paste(100*Plot_BR$thresholds, "%", sep="")
#
library(ggpmisc)
B_R_T <- ggplot(Plot_BR, aes(x=B_Richness, y=funcMaxed))+geom_point(size=1.5,alpha=0.6,colour="#377EB8")+
  stat_smooth(method="glm", 
              method.args = list(family=quasipoisson(link="identity")),
              colour="#377EB8",fill="#377EB8", linewidth=1,alpha=0.3) +
  ylab(expression("Number of Functions " >=Threshold)) +
  xlab("Bacterial richness index") +
  facet_wrap(~percent ,scales = 'free_y')+
  stat_poly_eq(use_label(c("adj.R2", "p.value.label"),sep = "*\", \"*"),family = "serif",size=3.5,colour="#377EB8")+
  theme_bw(base_size=12) +
  theme(axis.title = element_text(
    family = "serif",##
    size=12, ##
    lineheight = 0),#
    axis.text = element_text(
      family = "serif",#
      color="black",
      size=12),strip.text = element_text(family = "serif",face='bold',size=12))+
  #scale_y_continuous(limits=c(0,30))+
  theme(strip.background = element_rect(fill=c("#BDBDBD")))
B_R_T
ggsave('B_R_T_1.pdf',B_R_T,width=8,height=6)
ggsave("B_R_T_1.TIFF",B_R_T,width = 5,height = 4)
######
#####
#####
#####
#####
#####
####CAN####
setwd("F:/JGS-SMT微生物测序数据/功能基因")
df <- read.csv("JGS_Multi_CAN.csv",header = T,row.names=1, sep= ",")
df_scale <- scale(df[,5:43])
str(df_scale)
head(df_scale)
df_scale <- as.data.frame(df_scale)
group <- df[,1:4]
library(tidyverse)
df_CAN <- mutate(df_scale,group)
write.csv(df_CAN,'df_scale_CAN.csv')
#
df_CAN <- read.csv("df_scale_CAN.csv",header = T,row.names=1, sep= ",")

#SlectVars<-c('SWC','pH','NH4','NO3','TC','TN','SCN','DOC','DON','DOCN','AP','TP','SNP','LC','LN','LCN',
#'ACP', 'Catalase', 'BG', 'Sucrase', 'Urease','C_degra', 'C_fixa', 'Me_Meta', 
#'Denitri', 'Nitri', 'P_cycling', 'S_cycling')
#Soil_vars <- Index[,SoilVars]
#计算方差膨胀因子
library(tidyverse)
library(caret)
#
dim(df_CAN)
model1 <- lm(meanFunction~SWC+pH+NH4+NO3+TC+TN+SCN+DOC+DON+DOCN+AP+TP+SNP+LM+LC+LN+LCN+ FRB+B_Shannon+B_Richness+F_Shannon+F_Richness
             ,data=df_CAN)
summary(model1)
car::vif(model1)
#
dim(df_CAN)
#
model1_1 <- lm(meanFunction~SWC+pH+LM+LC+LN+DOC+NH4+SCN+B_Richness+F_Richness
               ,data=df_CAN)
summary(model1_1)
car::vif(model1_1)
#
library(piecewiseSEM)
library(tidyverse)
library(lme4)
library(tidyr)
library(Matrix)
library(nlme)
library(lmerTest)
library(dplyr)
#
library(emmeans)
CANnutr <- lmer(SCN ~ factor(N_rate)+(1|Block), data=df_CAN)
summary(CANnutr)
anova(CANnutr)
ranova(CANnutr)
emmeans(CANnutr,specs=pairwise~"N_rate")
#
CAN_SCN <- lm(SCN~N_rate,data=df_CAN)
summary(CAN_SCN)
anova(CAN_SCN)
emmeans(CAN_SCN,specs=pairwise~"N_rate")
#
#
#
#
########
#soil env
soilenvM <- lmer(meanFunction ~ SWC+pH+(1|Block), data=df_CAN)
summary(soilenvM)
anova(soilenvM)
ranova(soilenvM)
a1<-summary(soilenvM)$coefficients[2,1]
a2<-summary(soilenvM)$coefficients[3,1]
a0<-summary(soilenvM)$coefficients[1,1]
soilenvM=a1*df_CAN$SWC+a2*df_CAN$pH+a0
#
soilenvF <- lmer(F_Richness ~ SWC+pH+(1|Block), data=df_CAN)
summary(soilenvF)
anova(soilenvF)
ranova(soilenvF)
b1<-summary(soilenvF)$coefficients[2,1]
b2<-summary(soilenvF)$coefficients[3,1]
b0<-summary(soilenvF)$coefficients[1,1]
soilenvF=b1*df_CAN$SWC+b2*df_CAN$pH+b0
#
soilenvB <- lmer(B_Richness ~ SWC+pH+(1|Block), data=df_CAN)
summary(soilenvB)
anova(soilenvB)
ranova(soilenvB)
c1<-summary(soilenvB)$coefficients[2,1]
c2<-summary(soilenvB)$coefficients[3,1]
c0<-summary(soilenvB)$coefficients[1,1]
soilenvB=c1*df_CAN$SWC+c2*df_CAN$pH+c0
#Soilnutr
SoilnutrM <- lmer(meanFunction ~ SCN+DOC+NH4+(1|Block), data=df_CAN)
summary(SoilnutrM)
anova(SoilnutrM)
ranova(SoilnutrM)
d1<-summary(SoilnutrM)$coefficients[2,1]
d2<-summary(SoilnutrM)$coefficients[3,1]
d3<-summary(SoilnutrM)$coefficients[4,1]
d0<-summary(SoilnutrM)$coefficients[1,1]
SoilnutrM = d1*df_CAN$SCN+d2*df_CAN$DOC+d3*df_CAN$NH4+d0
#
SoilnutrF <- lmer(F_Richness ~ SCN+DOC+NH4+(1|Block), data=df_CAN)
summary(SoilnutrF)
anova(SoilnutrF)
ranova(SoilnutrF)
e1<-summary(SoilnutrF)$coefficients[2,1]
e2<-summary(SoilnutrF)$coefficients[3,1]
e3<-summary(SoilnutrF)$coefficients[4,1]
e0<-summary(SoilnutrF)$coefficients[1,1]
SoilnutrF = e1*df_CAN$SCN+e2*df_CAN$DOC+e3*df_CAN$NH4+e0
#
SoilnutrB <- lmer(B_Richness ~ SCN+DOC+NH4+(1|Block), data=df_CAN)
summary(SoilnutrB)
anova(SoilnutrB)
ranova(SoilnutrB)
f1<-summary(SoilnutrB)$coefficients[2,1]
f2<-summary(SoilnutrB)$coefficients[3,1]
f3<-summary(SoilnutrB)$coefficients[4,1]
f0<-summary(SoilnutrB)$coefficients[1,1]
SoilnutrB = f1*df_CAN$SCN+f2*df_CAN$DOC+f3*df_CAN$NH4+f0
#
#litter
#
LitterM <- lmer(meanFunction ~ LM+LC+LN+(1|Block), data=df_CAN)
summary(LitterM)
anova(LitterM)
ranova(LitterM)
g1<-summary(LitterM)$coefficients[2,1]
g2<-summary(LitterM)$coefficients[3,1]
g3<-summary(LitterM)$coefficients[4,1]
g0<-summary(LitterM)$coefficients[1,1]
LitterM=g1*df_CAN$LM+g2*df_CAN$LC+g3*df_CAN$LN+d0
#
LitterB <- lmer(B_Richness ~ LM+LC+LN+(1|Block), data=df_CAN)
summary(LitterB )
anova(LitterB )
ranova(LitterB)
h1<-summary(LitterB)$coefficients[2,1]
h2<-summary(LitterB)$coefficients[3,1]
h3<-summary(LitterB)$coefficients[4,1]
h0<-summary(LitterB)$coefficients[1,1]
LitterB=h1*df_CAN$LM+h2*df_CAN$LC+h3*df_CAN$LN+h0
#
LitterF <- lmer(F_Richness ~ LM+LC+LN+(1|Block), data=df_CAN)
summary(LitterF )
anova(LitterF )
ranova(LitterF)
i1<-summary(LitterF)$coefficients[2,1]
i2<-summary(LitterF)$coefficients[3,1]
i3<-summary(LitterF)$coefficients[4,1]
i0<-summary(LitterF)$coefficients[1,1]
LitterF=i1*df_CAN$LM+i2*df_CAN$LC+i3*df_CAN$LN+i0

#合并数据
str(df_CAN)
df_CAN01 <- mutate(df_CAN,SoilnutrM,SoilnutrB,SoilnutrF,soilenvM,soilenvF,soilenvB,LitterM,LitterB,LitterF)

#------------------------- fit the model using psem()
#
CAN_SEM<-psem(
  lme(meanFunction ~ SoilnutrM+soilenvM+LitterM+F_Richness+B_Richness,random=~1|Block, data=df_CAN01),
  lme(B_Richness ~ SoilnutrB+soilenvB+LitterB,random=~1|Block, data=df_CAN01),
  lme(F_Richness ~SoilnutrF+soilenvF+LitterF,random=~1|Block, data=df_CAN01),
  B_Richness %~~% F_Richness,
  data=df_CAN01
)%>%dSep(conditioning=T)
#
CAN_SEM04<-psem(
  lme(meanFunction ~ SoilnutrM+soilenvM+LitterM+F_Richness+B_Richness,random=~1|Block, data=df_CAN01),
  lme(B_Richness ~ SoilnutrB+soilenvB+LitterB,random=~1|Block, data=df_CAN01),
  lme(F_Richness ~ SoilnutrF+soilenvF+LitterF,random=~1|Block, data=df_CAN01),
  B_Richness %~~% F_Richness,
  data=df_CAN01
)
summary(CAN_SEM04)
#
CAN_SEM05<-psem(
  lme(meanFunction ~ soilenvM+SoilnutrM+LitterM+F_Richness+B_Richness,random=~1|Block, data=df_CAN01),
  lme(B_Richness ~ SoilnutrB+soilenvB+LitterB,random=~1|Block, data=df_CAN01),
  lme(F_Richness ~ SoilnutrF+LitterF,random=~1|Block, data=df_CAN01),
  B_Richness %~~% F_Richness,
  data=df_CAN01
)
summary(CAN_SEM05)
plot(CAN_SEM05)

##########
##########
##########
##########
####CAN
setwd("F:/JGS-SMT微生物测序数据/功能基因")
df <- read.csv("JGS_Multi_CAN.csv",header = T,row.names=1, sep= ",")
df_scale <- scale(df[,5:43])
str(df_scale)
head(df_scale)
df_scale <- as.data.frame(df_scale)
group <- df[,1:4]
library(tidyverse)
df_CAN <- mutate(df_scale,group)
write.csv(df_CAN,'df_scale_CAN.csv')
df_CAN <- read.csv("df_scale_CAN.csv",header = T,row.names=1)
library(ggplot2)  
library(vegan)
library(randomForest)
library(rfPermute)

RanForVars<-c("SWC", "pH", "NH4", "NO3", "TC","TN", "DOC", "DON","AP","TP","SCN","DOCN","SNP","LM","LC","LN","LCN",
              "Catalase","BG","Sucrase","Urease","ACP","C_degra","C_fixa","Me_Meta","Denitri","Nitri","P_cycling",
              "S_cycling","B_PCoA1","B_PCoA2","F_PCoA1","F_PCoA2","B_Richness","F_Richness","meanFunction","N_rate")

Soil<-which(names(df_CAN) %in% RanForVars)
RanForVars
RanForVars <- df_CAN[,RanForVars]
str(RanForVars)
RanForVars <- as.data.frame(RanForVars)
RanForVars %>% mutate(N_rate=as.factor(N_rate)) -> RanForVars_CAN
str(RanForVars_CAN)
####
#
set.seed(1234)###% Var explained: 79.63
Function_rf_CAN<-randomForest(meanFunction~., data =RanForVars_CAN,importance = TRUE, ntree = 2000,proximity=TRUE)
Function_rf_CAN
#
library(rfUtilities)
set.seed(1234)#Model signifiant at p = 0.01, Model R-square:  0.801
rf_perm_CAN <-rf.significance(Function_rf_CAN, RanForVars_CAN[,-c(36)], nperm=99, ntree = 2000)
rf_perm_CAN
#
set.seed(1234)
rf_results<-rfPermute(RanForVars_CAN$meanFunction~., data =RanForVars_CAN, nrep=99,importance = TRUE, ntree = 2000,  num.cores = 1)#
rf_results
##
rf_import_CAN <- importance(rf_results, sort.by = NULL ,decreasing = TRUE)#scale = TRUE
rf_import_CAN
#rf_import_UAN与predictor_var二选一;predictor_var<- data.frame(importance(rf_results, scale = TRUE), check.names = FALSE)
rf_CAN_plot <- rf_import_CAN %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) 
rf_CAN_plot
rf_CAN_plot <- rf_CAN_plot[order(rf_CAN_plot$`X.IncMSE`,decreasing=T),]#降序排列
rf_CAN_plot

#
Label<-as.data.frame(rep("Canopy addition of N",36))
rf_CAN_plot<-cbind(rf_CAN_plot,Label)
ID<-as.data.frame(c(1:36))
rf_CAN_plot<-cbind(rf_CAN_plot,ID)
rf_CAN_plot
str(rf_CAN_plot)

###
library(ggplot2)
library(ggrepel)
plot_CAN <- ggplot(data=rf_CAN_plot, aes(x = names, y = X.IncMSE))+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 2,label = label),family="serif",size=5)+
  labs(x = " ", y = "Increase in MSE (%)")+
  coord_flip()+
  facet_grid( ~ rf_CAN_plot$`rep("Canopy addition of N", 36)`, drop=TRUE,scale="free",space="free_x")+
  theme_bw(base_size=12)+
  theme(axis.title = element_text(
    family = "serif",#
    size=12, #
    lineheight = 0),#
    axis.text = element_text(
      family = "serif",
      color="black",
      size=12),
    strip.text = element_text(family = "serif",face='bold',size=14))+
  geom_text(aes(x = 1, y = 13,
                label = "Average multifunctionality"),stat="unique",fontface ="bold",
            size = 4,family = "serif")+
  geom_text(aes(x = 1.8, y = 18,
                label = "R^2 = 0.787"),stat="unique",fontface ="bold",
            size = 4,family = "serif")+
  geom_text(aes(x = 2.6, y = 19,
                label = "P = 0.01"),stat="unique",fontface ="bold",
            size = 4,family = "serif")

plot_CAN
ggsave('MSE_CAN.pdf',plot_CAN,width=5,height=7)
ggsave("MSE_CAN.TIFF",plot_CAN,width = 5,height =7)
##
#####
#####
#####
#####
#####
library(ggplot2)  
library(vegan)
library(randomForest)
library(rfPermute)
library(tidyverse)
setwd("F:/JGS-SMT微生物测序数据/随机森林预测物种")
###
FunPly <- read.csv("JGSfun_All.csv",header = T,check.names=F)
head(FunPly)
## 
Fun_spe <- FunPly[-c(1:11)]
head(Fun_spe)
tax <- FunPly[1:11]
head(tax)
#
##
group <- read.csv("jgsgroupFuntotalNew.csv",header = T,check.names=F)
head(group)
group$treat <- factor(group$treat,`levels`  = c("CK","CALN","CAHN","UALN","UAHN"))
group$Sampling_time <- factor(group$Sampling_time,levels = c("Jan., 2020","July, 2020","Jan., 2021","July, 2021"))
str(group)
#
phy01 <- Fun_spe %>%
  group_by(tax$Genus) %>% # 
  summarise_all(sum) %>%
  rename(Genus = `tax$Genus`)
phy01
write.csv(phy01,'JGS_Fun_Genus_data.csv')
#
Fun <- read.csv("JGS_Fun_Genus_data.csv",header = T,row.names=1)
OrderFun_CAN<-c("A13", "A4", "A11", "E1","A7","E11","A21","E7","A20","A6","I28","E8", "I15","A1", "I13", 
                "A12","M19","I12", "M14","M20","M21","M22","E19","I6","M8","M28","I4", "I11","I18","I26", 
                "A15","A18","A19","M15","M18","M25","M26","E14","E20","I8","M11","I19", "I21","I22","I25")#

FunNum<-which(names(Fun) %in% OrderFun_CAN)
FunNum
OrderFun_CAN <- Fun[,OrderFun_CAN]
OrderFun_CAN
write.csv (OrderFun_CAN, file ="OrderFunGenus_CAN.csv")
str(OrderFun_CAN)
###
FunCAN_t <-t(OrderFun_CAN)
#
Fun_group <- read.csv("JGS_CAN_微生物RFM.csv",header = T,row.names=1)
Fun_group
#
data_CAN <- cbind(Fun_group,FunCAN_t)
write.csv (data_CAN, file ="JGS_CAN_Genus_真菌预测.csv")
#
data_CAN_1 <- read.csv("JGS_CAN_Genus_真菌预测.csv",header = T,row.names=1)
data_CAN_1
#去除1-4列
data_CAN01 <- data_CAN_1[,-c(1:4)]
str(data_CAN01)
#data_UAN01 <- as.matrix(data_UAN01)
##
####
set.seed(1234)###% Var explained: 15.69
Var_rf_CAN<-randomForest(meanFunction~., data =data_CAN01,importance = TRUE, proximity=TRUE, ntree = 1000)
Var_rf_CAN
#
library(rfUtilities)
set.seed(1234)#Model signifiant at p = 0.01, Model R-square:  0.15
rf_perm_CAN <-rf.significance(Var_rf_CAN, data_CAN01[,-c(1)], nperm=99, ntree = 1000)
rf_perm_CAN
#
set.seed(1234)
rf_results<-rfPermute(data_CAN01$meanFunction~., data =data_CAN01, nperm=99,importance = TRUE, ntree = 1000,  num.cores = 1)#
rf_results
##
rf_import_CAN <- importance(rf_results, sort.by = NULL, scale = TRUE,decreasing = TRUE)
rf_import_CAN
#
rf_CAN_plot <- rf_import_CAN %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) 
rf_CAN_plot
rf_CAN_plot <- rf_CAN_plot[order(rf_CAN_plot$`X.IncMSE`,decreasing=T),]#降序排列
rf_CAN_plot
write.csv (rf_CAN_plot, file ="rf_CAN_plot绘图01.csv")
#pearman
#spearman
#spearman
#spearman
#spearman
Fun <- read.csv("CAN_Fun_spearman.csv",header = T,row.names=1)
Fun <- t(Fun)
Fun
#
Bac_group <- read.csv("JGS_CAN_微生物RFM.csv",header = T,row.names=1)
Bac_group
##相
library(psych)
res <- corr.test(Fun,Bac_group$meanFunction,use="pairwise",
                 method="spearman",adjust="fdr",
                 alpha=0.05)
res
res$p.adj
res$r
write.csv (res$p.adj, file ="CAN_Fun_p_1.adj.csv")
write.csv (res$r, file ="CAN_Fun_R_1.adj.csv")
##






