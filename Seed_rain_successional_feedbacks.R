##############################################################################
############Seed rain-successional feedbacks in wet tropical forests##########
##############################################################################
#Author Huanca Nunez created on 03/09/2020
################ 
library(readr)
library(corrplot)
library(datasets)
library(dplyr)
library(ecodist)
library(effects)
library(gridExtra)
library(grid)
library(graphics)
library(ggplot2)
library(ggpubr)
library(legendMap)
library(lm4)
library(MASS)
library(magrittr) 
library(SpadeR)
library(sjPlot)
library(sf)
require(stringr) 
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(reshape)
library(maptools)
library(picante)
library(tidyr)
library(vegan)

####################
## Study site Map###
#####################
##Map La Selva Biological Station, Costa Rica

#(devtools::install_github("dkahle/ggmap"))
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-85.15, -83.12), ylim = c(10, 11.15), expand = FALSE)


cordenates_map <- read.csv("~cordenates_map.csv")
data=cordenates_map
color=c("red","red","red","red","blue")
str(data)
library(ggmap)
citation("ggmap")
register_google(key = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
al1 = get_map(location = c(lon = -84.0780, lat = 10.40), zoom = 12, maptype = 'satellite')
al1MAP = ggmap(al1)
al1MAP + geom_point(data = data, aes(x =x, y = y), colour =color, size = 5)+
  geom_text(data = data, aes(x, y+0.006, label = Plot),size =8.5,fontface="bold",colour="white")+
  annotate("text", x=-84.02, y=10.43, label = "La Selva", colour = "white", size = 3,fontface="bold")+
  annotate("text", x=-84.02, y=10.426, label = "Biological Station", colour = "white", size = 3,fontface="bold")+
  geom_vline(xintercept=-84.16, linetype="dashed", color = "grey")+
  scale_bar(lon = -84.15, lat = 10.30,
            distance_lon = 2, distance_lat = 0.5, distance_legend = 1,
            dist_unit = "km",
            arrow_length = 2.5,arrow_distance = 1.5, arrow_north_size = 12)+
  xlab("Longitude") +
  ylab("Latitude")+
  theme(axis.text.y = element_text(size=18),axis.text.x   = element_text(size=18),
        axis.title.y =element_text(size=20),axis.title.x =element_text(size=20))

al1MAP


###################
### Figure 2 NMDS###
####################

#1. Long format data to matrix
Wholedata <- read_csv("seed_rain.csv")

WholedatabyTrap=Wholedata %>% group_by(Data,uniqueplace,TrapN,Code2) %>% summarize(NSeeds=sum(seeds_N))

WholedatabyTrap$plotYtrap=paste(WholedatabyTrap$uniqueplace,WholedatabyTrap$TrapN)
WholedatabyTrap$NSeeds=as.integer(WholedatabyTrap$NSeeds)

WholedatabyTrap2=WholedatabyTrap[,4:6]
colnames(WholedatabyTrap2)

WholedatabyTrap3=spread(WholedatabyTrap2, Code2, NSeeds, fill=0)
WholedatabyTrap4=WholedatabyTrap3[,2:179] # short form data or matrix
rownames(WholedatabyTrap4)=WholedatabyTrap3$plotYtrap 
write.csv(WholedatabyTrap4,"matrix_data.csv")

comm1 <- metaMDS(WholedatabyTrap4, dist = "chao") 

#
WholedatabyTrap_envi=unique(WholedatabyTrap[c("plotYtrap", "uniqueplace","Data")])
WholedatabyTrap_envi$group=c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10))
WholedatabyTrap_envi$color=c(rep("red",10),rep("red",10),rep("red",10),rep("red",10),rep("blue",10),rep("blue",10),rep("blue",10),rep("blue",10),rep("black",10))
WholedatabyTrap_envi$Site=c(rep("LS1",10),rep("TI1",10),rep("PS1",10),rep("CR1",10),rep("LS2",10),rep("TI2",10),rep("PS2",10),rep("CR2",10),rep("TPP2",10))
WholedatabyTrap_envi$Place=c(rep("LS",10),rep("TI",10),rep("PS",10),rep("CR",10),rep("LS",10),rep("TI",10),rep("PS",10),rep("CR",10),rep("PP",10))


fig=ordiplot(comm1, display = "sites", type = "none", main="NMDS ordination of seed rain by forests age") 


Specieslist=as.data.frame(fig$species)
Specieslist$Code2=rownames(Specieslist)
Specieslist_envi=unique(Wholedata[c("Code2", "LifeHistory")])
Specieslist_envi$Code2[duplicated(Specieslist_envi$Code2)]
Specieslist_envi2=merge(Specieslist_envi,Specieslist, by= "Code2", all=T)
Specieslist_envi2$group2=NA
Specieslist_envi2$group2[(Specieslist_envi2$LifeHistory== "ES" )] = "violet"
Specieslist_envi2$group2[(Specieslist_envi2$LifeHistory== "NES" )] = "darkgreen"

ord <- cca(WholedatabyTrap4 ~ uniqueplace, WholedatabyTrap_envi)

####
par(mai=c(0.8,0.8,0.5,0.5))
ordiplot (ord, display = 'si', type = 'n', ylim = c(-1,1),xlim = c(-1.5,1.5),xlab="NMDS1",xaxs="i", ylab="NMDS2")
ordiellipse(comm1, WholedatabyTrap_envi$group,col = colores,kind="se")
text(comm1,display="species", cex=0.5, col = Specieslist_envi2$group2)
abline(h=0,lty = "dotted"); abline(v=0,lty = "dotted") 

temp <- locator(1) # On the chart, click where you would like the text to appear
text(temp,"M",col = "Black",font = 2 )

legend(0.89,1.38,#"topright", , # position
       legend = c("  Years", "A1 = 12 (LS)","B1 = 15","C1 = 20 (LS)","D1 = 25","A2 = 32 (LS)", "B2 = 35", "C2 = 40 (LS)","D2 = 45", "M = Mature", "Forest   (LS)"), 
       text.col = c("black", "red","red","red","red","blue","blue","blue","blue","black","black"),
       cex = 0.66,
       bty = "n") # border

stressplot(comm1)


###Permanova & betadispers

sp2=WholedatabyTrap4
Env_wholetrap4<- read_csv("Env-wholetrap4.csv")
Code=rownames(sp2)
Years=c(rep(12,10), rep(15,10), rep(20,10),rep(25,10),rep(32,10),rep(35,10),rep(40,10),rep(45,10),rep(100,10))
Env_wholetrap4$code=Code
Env_wholetrap4$Years=Years
SpeEnv=cbind(sp2,Env_wholetrap4)


BC.dist=vegdist(SpeEnv[,1:178], distance="chao")
disp.age = betadisper(BC.dist, SpeEnv$site_time)
anova(disp.age)
TukeyHSD(disp.age)

# significance by site_time effect
Adonis=adonis(SpeEnv[,1:178]~SpeEnv$site_time, perm=999, by = "margin") 

#capture.output(Adonis,file="test.doc")
output_adonis=pairwise.adonis(SpeEnv[,1:178],SpeEnv$site_time, sim.method = "chao", p.adjust.m = "holm") 

################################
##### Figure 3 = Null model #####
###############################
library(picante)
library(ecodist)
#### Modified R code for null model analysis (we used as base the Stegen et al. (2013) code) ####

raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model	
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##make_null:
  
  ##looping over each pairwise community combination:
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_chao<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); 
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; 
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); 
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; 
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); # null.spXsite;
        
        ##calculate null chao
        null_chao[i] = distance(null.spXsite,method='chao');
        
      }; # end reps loop
      
      ## empirically observed chao
      obs.chao = distance(spXsite[c(null.one,null.two),],method='chao');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_chao==obs.chao);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_chao<obs.chao);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##here the modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      
      print(c(null.one,null.two,date()));
      
    }; ## end null.two loop
    
  }; ## end null.one loop
  
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }	
  
  return(results)
  
}


Nseed_by_specie <- read.csv("matrix_data.csv")

data=Nseed_by_specie[c(1:10),] ## 10 traps from site (similarity will be compared across all of them)
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
A1RC=raup_crick_abundance (data3)
A1RC2=t(as.data.frame.list(A1RC))
rownames(A1RC2)=NULL
colnames(A1RC2)="index"
A1RC2=as.data.frame(A1RC2)

data=Nseed_by_specie[c(11:20),] ## 10 traps from site 2
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
B1RC=raup_crick_abundance (data3)
B1RC2=t(as.data.frame.list(B1RC))
rownames(B1RC2)=NULL
colnames(B1RC2)="index"
B1RC2=as.data.frame(B1RC2)

data=Nseed_by_specie[c(21:30),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
C1RC=raup_crick_abundance (data3)
C1RC2=t(as.data.frame.list(C1RC))
rownames(C1RC2)=NULL
colnames(C1RC2)="index"
C1RC2=as.data.frame(C1RC2)

data=Nseed_by_specie[c(31:40),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
D1RC=raup_crick_abundance (data3)
D1RC2=t(as.data.frame.list(D1RC))
rownames(D1RC2)=NULL
colnames(D1RC2)="index"
D1RC2=as.data.frame(D1RC2)

data=Nseed_by_specie[c(41:50),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
A2RC=raup_crick_abundance (data3)
A2RC2=t(as.data.frame.list(A2RC))
rownames(A2RC2)=NULL
colnames(A2RC2)="index"
A2RC2=as.data.frame(A2RC2)

data=Nseed_by_specie[c(51:60),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
B2RC=raup_crick_abundance (data3)
B2RC2=t(as.data.frame.list(B2RC))
rownames(B2RC2)=NULL
colnames(B2RC2)="index"
B2RC2=as.data.frame(B2RC2)

data=Nseed_by_specie[c(61:70),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
C2RC=raup_crick_abundance (data3)
C2RC2=t(as.data.frame.list(C2RC))
rownames(C2RC2)=NULL
colnames(C2RC2)="index"
C2RC2=as.data.frame(C2RC2)

data=Nseed_by_specie[c(71:80),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
D2RC=raup_crick_abundance (data3)
D2RC2=t(as.data.frame.list(D2RC))
rownames(D2RC2)=NULL
colnames(D2RC2)="index"
D2RC2=as.data.frame(D2RC2)


A1RC2$plot="A1"
A1RC2$age=12
A1RC2$data=1
A1RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(A1RC2)=c("index","plot","age","data","unit")

B1RC2$plot="B1"
B1RC2$age=15
B1RC2$data=1
B1RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(B1RC2)=c("index","plot","age","data","unit")

C1RC2$plot="C1"
C1RC2$age=20
C1RC2$data=1
C1RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(C1RC2)=c("index","plot","age","data","unit")

D1RC2$plot="D1"
D1RC2$age=25
D1RC2$data=1
D1RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(D1RC2)=c("index","plot","age","data","unit")

A2RC2$plot="A2"
A2RC2$age=33
A2RC2$data=2
A2RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(A2RC2)=c("index","plot","age","data","unit")

B2RC2$plot="B2"
B2RC2$age=36
B2RC2$data=2
B2RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(B2RC2)=c("index","plot","age","data","unit")

C2RC2$plot="C2"
C2RC2$age=42
C2RC2$data=2
C2RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(C2RC2)=c("index","plot","age","data","unit")


D2RC2$plot="D2"
D2RC2$age=46
D2RC2$data=2
D2RC2$unit=c(c(1:9),c(2:9),c(3:9),c(4:9),c(5:9),c(6:9),c(7:9),c(8:9),9)
colnames(D2RC2)=c("index","plot","age","data","unit")

RAUP=rbind(A1RC2,B1RC2,C1RC2,D1RC2,A2RC2,B2RC2,C2RC2,D2RC2)
write.csv(RAUP,"RAUP.csv")

RAUP$location1=(RAUP$unit)
RAUP$location2=(RAUP$unti2)
RAUP$diff_locat=abs(RAUP$location1-RAUP$location2) ##abs = absolute value

RAUP$dist_group=1
RAUP$dist_group[RAUP$diff_locat<50]="near dist"
RAUP$dist_group[RAUP$diff_locat>50&RAUP$diff_locat<130]="mid dist"
RAUP$dist_group[RAUP$diff_locat>130]="far dist"
RAUP$dist_group=as.factor(RAUP$dist_group)


##t-test confidence interval of each point
Raup_means=read.csv("RAUP_45.csv") # there is 45 comparisons for each site ⨯ time point combination.
Raup_meansA1=Raup_means[Raup_means$plot=="A1",]
errorA1 <- qt(0.975,df=length(Raup_meansA1$index)-1)*sd(Raup_meansA1$index)/sqrt(length(Raup_meansA1$index))
leftA1 <- mean(Raup_meansA1$index)-errorA1
rightA1 <- mean(Raup_meansA1$index)+errorA1

Raup_meansB1=Raup_means[Raup_means$plot=="B1",]
errorB1 <- qt(0.975,df=length(Raup_meansB1$index)-1)*sd(Raup_meansB1$index)/sqrt(length(Raup_meansB1$index))
leftB1 <- mean(Raup_meansB1$index)-errorB1
rightB1 <- mean(Raup_meansB1$index)+errorB1

Raup_meansC1=Raup_means[Raup_means$plot=="C1",]
errorC1 <- qt(0.975,df=length(Raup_meansC1$index)-1)*sd(Raup_meansC1$index)/sqrt(length(Raup_meansC1$index))
leftC1 <- mean(Raup_meansC1$index)-errorC1
rightC1 <- mean(Raup_meansC1$index)+errorC1

Raup_meansD1=Raup_means[Raup_means$plot=="D1",]
errorD1 <- qt(0.975,df=length(Raup_meansD1$index)-1)*sd(Raup_meansD1$index)/sqrt(length(Raup_meansD1$index))
leftD1 <- mean(Raup_meansD1$index)-errorD1
rightD1 <- mean(Raup_meansD1$index)+errorD1

Raup_meansA2=Raup_means[Raup_means$plot=="A2",]
errorA2 <- qt(0.975,df=length(Raup_meansA2$index)-1)*sd(Raup_meansA2$index)/sqrt(length(Raup_meansA2$index))
leftA2 <- mean(Raup_meansA2$index)-errorA2
rightA2 <- mean(Raup_meansA2$index)+errorA2

Raup_meansB2=Raup_means[Raup_means$plot=="B2",]
errorB2 <- qt(0.975,df=length(Raup_meansB2$index)-1)*sd(Raup_meansB2$index)/sqrt(length(Raup_meansB2$index))
leftB2 <- mean(Raup_meansB2$index)-errorB2
rightB2 <- mean(Raup_meansB2$index)+errorB2

Raup_meansC2=Raup_means[Raup_means$plot=="C2",]
errorC2 <- qt(0.975,df=length(Raup_meansC2$index)-1)*sd(Raup_meansC2$index)/sqrt(length(Raup_meansC2$index))
leftC2 <- mean(Raup_meansC2$index)-errorC2
rightC2 <- mean(Raup_meansC2$index)+errorC2

Raup_meansD2=Raup_means[Raup_means$plot=="D2",]
errorD2 <- qt(0.975,df=length(Raup_meansD2$index)-1)*sd(Raup_meansD2$index)/sqrt(length(Raup_meansD2$index))
leftD2 <- mean(Raup_meansD2$index)-errorD2
rightD2 <- mean(Raup_meansD2$index)+errorD2

summary_plot=as.data.frame(c(1:8))
summary_plot$mean=c(mean(Raup_meansA1$index),mean(Raup_meansB1$index),mean(Raup_meansC1$index),mean(Raup_meansD1$index),mean(Raup_meansA2$index),mean(Raup_meansB2$index),mean(Raup_meansC2$index),mean(Raup_meansD2$index))
summary_plot$CILeft=c(leftA1,leftB1,leftC1,leftD1,leftA2,leftB2,leftC2,leftD2)
summary_plot$CIRight=c(rightA1,rightB1,rightC1,rightD1,rightA2,rightB2,rightC2,rightD2)
summary_plot$Age=as.numeric(c(12,15,20,25,32,35,40,45))
str(summary_plot)
colnames(summary_plot)=c("N","mean","CILeft","CIRight","Age")
summary_plot=as.data.frame(summary_plot)

colnames(summary_plot)
p1=ggplot(summary_plot, aes(Age, mean)) + geom_point(data=Total_data1, aes(Age,index,colour = factor(dist_group)),alpha = 0.7)+ 
  geom_errorbar(aes(ymin=summary_plot$CILeft, ymax=summary_plot$CIRight), width=0.4) + theme_bw(base_size=12) +
  scale_color_manual("Distance",breaks = c("far dist", "mid dist", "near dist"),
                     values=c("blue", "purple", "red"),labels = c("Far(> 130 m)", "Mid (50 <130 m)", "Near (< 50 m)"))+
  scale_y_continuous(limits = c(-1.1,1.25))+
  scale_x_continuous(breaks =c(12,15,20,25,32,35,40,45))+
  labs(y = "Within successional plot dissimilarity (RCab)")+
  labs(x = "\nForest successional age ")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_hline(yintercept=1, linetype="dotted", color = "black")+
  geom_hline(yintercept=-1, linetype="dotted", color = "black")+
  annotate("text", x=30, y=1.1, label= "Deterministic (divergent) community assembly",size= 7) +
  annotate("text", x=30, y=-1.1, label= "Deterministic (convergent) community assembly",size= 7)+
  annotate("text", x=10, y=1.25, label= "(a)", size= 6)+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_smooth(data=Raup_means, mapping = aes(x = age, y = index), method = "lm")

chrono=p1+theme(legend.position="top",legend.background = element_rect(fill = "gray99"),legend.title = element_text(size=20),legend.text = element_text(size=20))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x   = element_text(size=24, angle=90),axis.text.y   = element_text(size=22),axis.title=element_text(size=20))


###Raup temporal turnover #####
# only ten comparisons for each site

data=Nseed_by_specie[c(1:10,41:50),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
A1RC=raup_crick_abundance (data3)
A1RC2=t(as.data.frame.list(A1RC))
rownames(A1RC2)=NULL
colnames(A1RC2)="index"
A1RC2=as.data.frame(A1RC2)
A1RC2$row=row.names(A1RC2)
A1RC3=A1RC2[c(10,29,47,64,80,95,109,122,134,145), ] 
rownames(A1RC3)=1:10
A1RC4=subset( A1RC3, select = -row )

data=Nseed_by_specie[c(11:20,51:60),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
B1RC=raup_crick_abundance (data3)
B1RC2=t(as.data.frame.list(B1RC))
rownames(B1RC2)=NULL
colnames(B1RC2)="index"
B1RC2=as.data.frame(B1RC2)
B1RC2$row=row.names(B1RC2)
B1RC3=B1RC2[c(10,29,47,64,80,95,109,122,134,145), ] 
rownames(B1RC3)=1:10
B1RC4=subset( B1RC3, select = -row )

data=Nseed_by_specie[c(21:30,61:70),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
C1RC=raup_crick_abundance (data3)
C1RC2=t(as.data.frame.list(C1RC))
rownames(C1RC2)=NULL
colnames(C1RC2)="index"
C1RC2=as.data.frame(C1RC2)
C1RC2$row=row.names(C1RC2)
C1RC3=C1RC2[c(10,29,47,64,80,95,109,122,134,145), ] 
rownames(C1RC3)=1:10
C1RC4=subset( C1RC3, select = -row )

data=Nseed_by_specie[c(31:40,71:80),]
data2<- as.data.frame(t(data[-1])) 
data2=data2[apply(data2[-1], 1, function(x) !all(x==0)),]
data3=as.data.frame(t(data2))
data3=cbind(a=data$a,data3)
D1RC=raup_crick_abundance (data3)
D1RC2=t(as.data.frame.list(D1RC))
rownames(D1RC2)=NULL
colnames(D1RC2)="index"
D1RC2=as.data.frame(D1RC2)
D1RC2$row=row.names(D1RC2)
D1RC3=D1RC2[c(10,29,47,64,80,95,109,122,134,145), ] 
rownames(D1RC3)=1:10
D1RC4=subset( D1RC3, select = -row )


A1RC4$plot="A1"
A1RC4$age=12
A1RC4$data=1
A1RC4$unit=1:10
colnames(A1RC4)=c("index","plot","age","data","unit")

B1RC4$plot="B1"
B1RC4$age=15
B1RC4$data=1
B1RC4$unit=1:10
colnames(B1RC4)=c("index","plot","age","data","unit")

C1RC4$plot="C1"
C1RC4$age=20
C1RC4$data=1
C1RC4$unit=1:10
colnames(C1RC4)=c("index","plot","age","data","unit")

D1RC4$plot="D1"
D1RC4$age=25
D1RC4$data=1
D1RC4$unit=1:10
colnames(D1RC4)=c("index","plot","age","data","unit")

RAUP=rbind(A1RC4,B1RC4,C1RC4,D1RC4)
write.csv(RAUP,"RAUP_temporal.csv")

##t-test confidence interval of each point
Raup_means=read.csv("RAUP_temporal_10.csv") # there is 10 comparisons for each site ⨯ time point combination.
Raup_meansA1=Raup_means[Raup_means$plot=="A1",]
errorA1 <- qt(0.975,df=length(Raup_meansA1$index)-1)*sd(Raup_meansA1$index)/sqrt(length(Raup_meansA1$index))
leftA1 <- mean(Raup_meansA1$index)-errorA1
rightA1 <- mean(Raup_meansA1$index)+errorA1

Raup_meansB1=Raup_means[Raup_means$plot=="B1",]
errorB1 <- qt(0.975,df=length(Raup_meansB1$index)-1)*sd(Raup_meansB1$index)/sqrt(length(Raup_meansB1$index))
leftB1 <- mean(Raup_meansB1$index)-errorB1
rightB1 <- mean(Raup_meansB1$index)+errorB1

Raup_meansC1=Raup_means[Raup_means$plot=="C1",]
errorC1 <- qt(0.975,df=length(Raup_meansC1$index)-1)*sd(Raup_meansC1$index)/sqrt(length(Raup_meansC1$index))
leftC1 <- mean(Raup_meansC1$index)-errorC1
rightC1 <- mean(Raup_meansC1$index)+errorC1

Raup_meansD1=Raup_means[Raup_means$plot=="D1",]
errorD1 <- qt(0.975,df=length(Raup_meansD1$index)-1)*sd(Raup_meansD1$index)/sqrt(length(Raup_meansD1$index))
leftD1 <- mean(Raup_meansD1$index)-errorD1
rightD1 <- mean(Raup_meansD1$index)+errorD1



summary_plot=as.data.frame(c(1:4))
summary_plot$mean=c(mean(Raup_meansA1$index),mean(Raup_meansB1$index),mean(Raup_meansC1$index),mean(Raup_meansD1$index))
summary_plot$CILeft=c(leftA1,leftB1,leftC1,leftD1)
summary_plot$CIRight=c(rightA1,rightB1,rightC1,rightD1)
summary_plot$Age=as.numeric(c(12,15,20,25))
str(summary_plot)
colnames(summary_plot)=c("N","mean","CILeft","CIRight","Age")
summary_plot=as.data.frame(summary_plot)

### linear plot 

p2=ggplot(summary_plot, aes(Age, mean)) + geom_point() + geom_errorbar(aes(ymin=CILeft, ymax=CIRight), width=0.4) + theme_bw(base_size=12) +
  scale_y_continuous(limits = c(-1.1,1.25))+
  scale_x_continuous(breaks =c(12,15,20,25),
                     label=c("A1 - A2", "B1 - B2","C1 - C2","D1 - D2"))+
  labs(y = "Among successional plot dissimilarity (RCab)")+
  labs(x = "")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_hline(yintercept=1, linetype="dotted", color = "black")+
  geom_hline(yintercept=-1, linetype="dotted", color = "black")+
  annotate("text", x=18, y=1.1, label= "Deterministic (divergent) community assembly",size= 7) +
  annotate("text", x=18, y=-1.1, label= "Deterministic (convergent) community assembly",size= 7)+
  annotate("text", x=10, y=1.25, label= "(b)", size= 6)+
  geom_smooth(data=Raup_means, mapping = aes(x = age, y = index), method = "lm")

temporal=p2+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                   axis.text.x   = element_text(size=24, angle=90),axis.text.y   = element_text(size=20),axis.title=element_text(size=18))

#  regression ##
chorno_B=lmer(index~age+(1|plot)+(1|dist_group), RAUP)
plot_model(chorno_B,type='diag') # checking assumptions
qqnorm(resid(chorno_B))
summary(chorno_B)
chorno <- update(chorno_B, . ~ . -(1|dist_group))
summary(chorno)

AIC(chorno,chorno_B) # simpler model is the best one.

temporal=lmer(index~age+(1|trap), RAUP_temporal)
plot_model(temporal,type='diag') # checking assumptions
qqnorm(resid(temporal))
summary(temporal)

##################################
###### Figure 4 -Similarity index## 
####################################

CR_data <- read.csv("matrix_data.csv")
Richness_chao=SimilarityMult(CR_data,"abundance",q=0,nboot=200,"relative") 
Simi_Riche=as.data.frame(Richness_chao[["pairwise"]][["C02"]])# presence/absence - Sorensen similarity index
Simi_Riche_matrix=Richness_chao[["similarity.matrix"]][["C02"]]
rownames(Simi_Riche_matrix)=c("A1(LS)","B1","C1(LS)","D1","A2(LS)","B2","C2(LS)","D2","M(LS)")
colnames(Simi_Riche_matrix)=c("A1(LS)","B1","C1(LS)","D1","A2(LS)","B2","C2(LS)","D2","M(LS)")
corrplot(Simi_Riche_matrix, method ="circle", type = "upper",addCoef.col = "black", cl.lim=c(0,1),tl.col=c("red","red","red","red","blue","blue","blue","blue","black"),tl.cex = 1.2,tl.srt=45,col=brewer.pal(n=8, name="RdBu"))


Similarity_abundance= SimilarityMult(CR_data,"abundance",q=1,nboot=200,"relative") # 
Sim_Abun=as.data.frame(Similarity_abundance[["pairwise"]][["C12"]]) # weighted abundance Horn index
Sim_Abun_matrix=Similarity_abundance[["similarity.matrix"]][["C12"]]
rownames(Sim_Abun_matrix)=c("A1(LS)","B1","C1(LS)","D1","A2(LS)","B2","C2(LS)","D2","M")
colnames(Sim_Abun_matrix)=c("A1(LS)","B1","C1(LS)","D1","A2(LS)","B2","C2(LS)","D2","M")
corrplot(Sim_Abun_matrix, method ="circle", type = "upper",addCoef.col = "black", cl.lim=c(0,1),tl.col=c("red","red","red","red","blue","blue","blue","blue","black"),tl.cex = 1.2,tl.srt=45,col=brewer.pal(n=8, name="RdBu"))


Similarity_abundance_dom= SimilarityMult(CR_data,"abundance",q=2,nboot=200,"relative")
Sim_domi=as.data.frame(Similarity_abundance_dom[["pairwise"]][["C22"]]) # community Morisita-Horn similarity index 
Sim_domi_matrix=Similarity_abundance_dom[["similarity.matrix"]][["C22"]]
rownames(Sim_domi_matrix)=c("A1(LS)","B1","C1(LS)","D1","A2(LS)","B2","C2(LS)","D2","M")
colnames(Sim_domi_matrix)=c("A1(LS)","B1","C1(LS)","D1","A2(LS)","B2","C2(LS)","D2","M")
corrplot(Sim_domi_matrix, method ="circle", type = "upper",addCoef.col = "black", cl.lim=c(0,1),tl.col=c("red","red","red","red","blue","blue","blue","blue","black"),tl.cex = 1.2,tl.srt=45,col=brewer.pal(n=8, name="RdBu"))

#1. Similarity to mature forests ##

### between sites similarity  ----

Simi_Rich_M=Simi_Riche[c(1,2,3,9,10,16,27,28,29,31,32,34),] # the number represent the time by plot combination.
Simi_Rich_M$Timet=as.factor(c("T1","T1","T1","T1","T1","T1","T2","T2","T2","T2","T2","T2"))
Simi_Rich_M$Sites=as.factor(c("A(LS) - B","A(LS) - C(LS)","A(LS) - D","B - C(LS)","B - D","C(LS) - D","A(LS) - B","A(LS) - C(LS)","A(LS) - D","B - C(LS)","B - D","C(LS) - D"))
Simi_Rich_M$Index="Richness"
Data=Simi_Rich_M

display.brewer.pal(9, "Set1")
brewer.pal(9, "Set1")

col_man1=c("#E7298A","#7570B3","#1B9E77","#E6AB02","black","#A65628")
col_man2=c("#1B9E77","#E7298A","#E6AB02","black")
p1=ggplot(Data, aes(x=Timet, y=Estimate,colour=Sites)) + 
  scale_y_continuous(limits=c(0,1))+
  geom_line(aes(group = Sites),size=0.7) +
  geom_point(size = 3,shape=1)+
  ylab("Similarity Index")+
  xlab("")+
  ggtitle("(d) Presence/Absence")+
  theme_bw()+
  theme(plot.title = element_text(size = 8.5))+
  scale_color_manual(values = col_man1)+
  scale_x_discrete(labels=c("T1" = "1997-1999", "T2" = "2015-2017"))

p1.1=p1+theme(axis.title.y =element_text(size = 8),axis.text =element_text(size = 7),
              legend.text =element_text(size = 6),legend.title =element_text(size = 7))



##
Sim_Abun_M=Sim_Abun[c(1,2,3,9,10,16,27,28,29,31,32,34),]
Sim_Abun_M$Timet=as.factor(c("T1","T1","T1","T1","T1","T1","T2","T2","T2","T2","T2","T2"))
Sim_Abun_M$Sites=as.factor(c("A(LS) - B","A(LS) - C(LS)","A(LS) - D","B - C(LS)","B - D","C(LS) - D","A(LS) - B","A(LS) - C(LS)","A(LS) - D","B - C(LS)","B - D","C(LS) - D"))
Sim_Abun_M$Index="Abundance"
Data=Sim_Abun_M
p2=ggplot(Data, aes(x=Timet, y=Estimate,colour=Sites)) + 
  scale_y_continuous(limits=c(0,1))+
  geom_line(aes(group = Sites),size=0.7) +
  geom_point(size = 3,shape=1)+
  ylab("")+
  xlab("")+
  ggtitle("(e)Abundance weighted")+
  theme_bw()+
  theme(plot.title = element_text(size = 8.5))+
  scale_color_manual(values = col_man1)+
  scale_x_discrete(labels=c("T1" = "1997-1999", "T2" = "2015-2017"))

p2.1=p2+theme(axis.title.y =element_text(size = 8),axis.text =element_text(size = 7),
              legend.text =element_text(size = 6),legend.title =element_text(size = 7))

###
Sim_domi_M=Sim_domi[c(1,2,3,9,10,16,27,28,29,31,32,34),]
Sim_domi_M$Timet=as.factor(c("T1","T1","T1","T1","T1","T1","T2","T2","T2","T2","T2","T2"))
Sim_domi_M$Sites=as.factor(c("A(LS) - B","A(LS) - C(LS)","A(LS) - D","B - C(LS)","B - D","C(LS) - D","A(LS) - B","A(LS) - C(LS)","A(LS) - D","B - C(LS)","B - D","C(LS) - D"))
Sim_domi_M$Index="Dominance"

Data=Sim_domi_M
p3=ggplot(Data, aes(x=Timet, y=Estimate,colour=Sites)) + 
  scale_y_continuous(limits=c(0,1))+
  geom_line(aes(group = Sites),size=0.7,position=position_dodge(w=0.04)) +
  geom_point(size = 3,shape=1)+
  ylab("")+
  xlab("")+
  ggtitle("(f) Dominant species")+
  theme_bw()+
  theme(plot.title = element_text(size = 8.5))+
  scale_color_manual(values = col_man1)+
  scale_x_discrete(labels=c("T1" = "1997-1999", "T2" = "2015-2017"))

p3.1=p3+theme(axis.title.y =element_text(size = 8),axis.text =element_text(size = 7),
              legend.text =element_text(size = 6),legend.title =element_text(size = 7))


a=ggarrange(p1.1, p2.1, p3.1, ncol=3, nrow=1, common.legend = TRUE, legend="top")
tiff(file = "between.tiff", width = 6, height = 3.5,units = "in", res = 300)
plot(a)
dev.off()

stat.test=t_test(Data,Estimate ~ Time) %>% adjust_pvalue(method = "BH")

T_Rich=Simi_Rich_M
Ttest1=stat.test=t_test(T_Rich,Estimate ~ Timet,paired=TRUE)

T_Abun=Sim_Abun_M
Ttest2=stat.test=t_test(T_Abun,Estimate ~ Timet,paired=TRUE)

T_domi=Sim_domi_M
Ttest3=stat.test=t_test(T_domi,Estimate ~ Timet,paired=TRUE)

### similarity to mature forest-----

Simi_Riche$N=c(1:36)
Simi_Rich_M=Simi_Riche[c(8,15,21,26,30,33,35,36),]
Simi_Rich_M$Timet=as.factor(c("T1","T1","T1","T1","T2","T2","T2","T2"))
Simi_Rich_M$Sites=as.factor(c("A(LS)", "B","C(LS)","D","A(LS)", "B","C(LS)","D"))
Simi_Rich_M$Index="Richness"
Data=Simi_Rich_M
p4=ggplot(Data, aes(x=Timet, y=Estimate,colour=Sites)) + 
  scale_y_continuous(limits=c(0,1))+
  #geom_errorbar(aes(ymin=Data$`95%.LCL`, ymax=Data$`95%.UCL`)) +
  geom_line(aes(group = Sites),size=0.7) +
  geom_point(size = 3,shape=1)+
  ylab("Similarity Index")+
  xlab("")+
  theme_bw()+
  ggtitle("(a) Presence/Absence")+
  theme(plot.title = element_text(size = 8.5))+
  scale_color_manual(values = col_man2)+
  scale_x_discrete(labels=c("T1" = "1997-1999", "T2" = "2015-2017"))

p4.2=p4+theme(axis.title.y =element_text(size = 8),axis.text =element_text(size = 7),
              legend.text =element_text(size = 6),legend.title =element_text(size = 7))

Simi_Rich_M$ages=c(12,15,20,25,32,35,40,45)
Simi_Rich_M2=Simi_Rich_M[5:8,]
data=lm(Simi_Rich_M2$Estimate~Simi_Rich_M2$ages)
summary(data) # not significnat 

##
Sim_Abun_M=Sim_Abun[c(8,15,21,26,30,33,35,36),]
Sim_Abun_M$Timet=as.factor(c("T1","T1","T1","T1","T2","T2","T2","T2"))
Sim_Abun_M$Sites=as.factor(c("A(LS)", "B","C(LS)","D","A(LS)", "B","C(LS)","D"))
Sim_Abun_M$Index="Abundance"
Data=Sim_Abun_M
p5=ggplot(Data, aes(x=Timet, y=Estimate,colour=Sites)) + 
  scale_y_continuous(limits=c(0,1))+
  #geom_errorbar(aes(ymin=Data$`95%.LCL`, ymax=Data$`95%.UCL`)) +
  geom_line(aes(group = Sites),size=0.7) +
  geom_point(size = 3,shape=1)+
  ylab("")+
  xlab("")+
  ggtitle("(b) Abundance weighted")+
  theme_bw()+
  theme(plot.title = element_text(size = 8.5))+
  scale_color_manual(values = col_man2)+
  scale_x_discrete(labels=c("T1" = "1997-1999", "T2" = "2015-2017"))

p5.2=p5+theme(axis.title.y =element_text(size = 8),axis.text =element_text(size = 7),
              legend.text =element_text(size = 6),legend.title =element_text(size = 7))
###
Sim_domi=as.data.frame(Similarity_abundance_dom[["pairwise"]][["C22"]]) # community Morisita-Horn similarity index 

Sim_domi_M=Sim_domi[c(8,15,21,26,30,33,35,36),]
Sim_domi_M$Timet=as.factor(c("T1","T1","T1","T1","T2","T2","T2","T2"))
Sim_domi_M$Sites=as.factor(c("A(LS)", "B","C(LS)","D","A(LS)", "B","C(LS)","D"))
Sim_domi_M$Index="Dominance"
Data=Sim_domi_M
p6=ggplot(Data, aes(x=Timet, y=Estimate,colour=Sites)) + 
  scale_y_continuous(limits=c(0,1))+
  #geom_errorbar(aes(ymin=Data$`95%.LCL`, ymax=Data$`95%.UCL`)) +
  geom_line(aes(group = Sites),size=0.7) +
  geom_point(size = 2.5,shape=1)+
  ylab("")+
  xlab("")+
  ggtitle("(c) Dominant species")+
  theme_bw()+
  theme(plot.title = element_text(size = 8.5))+
  scale_color_manual(values = col_man2)+
  scale_x_discrete(labels=c("T1" = "1997-1999", "T2" = "2015-2017"))

p6.2=p6+theme(axis.title.y =element_text(size = 8),axis.text =element_text(size = 7),
              legend.text =element_text(size = 6),legend.title =element_text(size = 7))

b=ggarrange(p4.2,p5.2,p6.2, ncol=3, nrow=1, common.legend = TRUE, legend="top")

tiff(file = "to_mature.tiff", width = 6, height = 3.5,units = "in", res = 300)
plot(b)
dev.off()

Simi_Rich_M$Age=c(12,15,20,25,32,35,42,45)
Sim_Abun_M$Age=c(12,15,20,25,32,35,42,45)
Sim_domi_M$Age=c(12,15,20,25,32,35,42,45)

model1 <- lm(Estimate~Age+Site,data=Simi_Rich_M)
plot_model(model1,type='diag') # checking assumptions
summary(model1)
model2 <- lm(Estimate~Age+Site,data=Sim_Abun_M)
plot_model(model2,type='diag') # checking assumptions
summary(model2)
model3 <- lm(Estimate~Age+Site,data=Sim_domi_M)
plot_model(model2,type='diag') # checking assumptions
summary(model2)

#####################################################
### Figure 7 - species richeness by funtional groups  
####################################################

#1. classic rarefactions using vegan

S <- specnumber(CR_data) # observed number of species
raremax <- min(rowSums(CR_data))
Srare <- rarefy(CR_data, raremax)
min_n<- min(rowSums(CR_data)) # minimun number of species use to rarefaction limit
max_n<- max(rowSums(CR_data)) # max number of species use to rarefaction limit

lty <- c("dotted","dotted","dotted","dotted","longdash","longdash","longdash","longdash","solid")
col<- c("red","red","red","red","blue","blue","blue","blue","black")
rarecurve(CR_data, step = 20, col = col,sample = min_n, lty=lty, cex = 0.7,ylab="Number of Species",xlab="Number of Individuals")
lty2=c("solid","longdash","dotted")
legend("bottomright",bty = "n",legend = c("Mature Forest","Late successional","Mid successional"),
       col = "black",lwd = 0.8, lty= lty2)


#2. Rarefaction with INEXT package
CR_data2=as.data.frame(t(CR_data))

#rarefy to the end ###
out <- iNEXT(CR_data2, q=c(0),datatype="abundance",endpoint=max_n) # endpoint you can modify until what number of individuals you want to rarify 
ggiNEXT(out,se=F,color.var="site")

# rarefy to the minimun number of individuals by functional group
CR_data2=read.csv("~matrix_Large_small_Seeds.csv")
min_n<- min(rowSums(CR_data2)) # minimum number of species use to rarefaction limit
max_n<- max(rowSums(CR_data2)) # max number of species use to rarefaction limit

out <- iNEXT(CR_data2,datatype="abundance",endpoint=min_n) # endpoint you can modify until what number of individuals you want to rarify 
write.csv(out$AsyEst,"Large_and_small.csv") ##
data=as.data.frame(data2[["AsyEst"]])
data=data[data$Diversity=="Species richness",]
data.1=data[c(1:9),] #large
data.2=data[c(10:18),] # small
data.1$total=data.1$Estimator+data.2$Estimator
data.2$total=data.1$total
data.1$large=data.1$Estimator/data.1$total
data.2$small=data.2$Estimator/data.2$total
data.1$Seed_size="Large"
data.2$Seed_size="Small"

Large_vs_small=merge(data.1,data.2, by=Place) # rownames is plot name here
write.csv(Large_vs_small,"Large_vs_small.csv")

speciesproportion <- read_csv("Large_vs_small.csv")
colnames(speciesproportion)[4]="Seed_size"
colnames(speciesproportion)
speciesproportion$Plot=factor(speciesproportion$Plot,levels=c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2","M"))
colnames(speciesproportion)=c( "Plot","mean","se","(a) Seed size", "Years","LCL", "UCL", "Place"  )
speciesproportion2=speciesproportion[-which(speciesproportion$Plot=="M"),]
speciesproportion3=speciesproportion[which(speciesproportion$Plot=="M"),]
speciesproportion3$Years=54

model4 <- lm(mean~Years+Place,data=speciesproportion)
plot_model(model4,type='diag') # checking assumptions
summary(model4)


text_high <- textGrob("M", gp=gpar(fontsize=10)) #fontface="bold"


p<-  speciesproportion2 %>% ggplot(aes(Years,mean, group=`(a) Seed size`,color=`(a) Seed size`,fill=`(a) Seed size`)) +
  stat_smooth(method="lm",se=T)+
  geom_point(size=2) +
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(breaks =c(12,15,20,25,32,35,40,45))+
  ylab("Percentage of species")+
  xlab("Successional Age (y)")+
  expand_limits(x = 55)+
  geom_point(data = speciesproportion3,size=2)+
  annotate("text", x=51, y=0.98, label= paste0("R^2== 0.99"), size=3,parse=TRUE)+
  annotate("text", x=51, y=0.88, label= "p < 0.001", size=3)#+
N1=p+
  theme_bw()+
  theme(plot.margin=unit(c(20.5, 5.5, 15.5, 5.5),"points"),panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = NULL),axis.text.y = element_text(size=10.5),
        axis.text.x   = element_text(size=10),axis.title=element_text(size=10),strip.text.x = element_text(size = 10),legend.position = c(0.5,1.17),legend.direction = "horizontal",legend.spacing.y = unit(0.5, 'cm'),legend.text=element_text(size=8),legend.title=element_text(size=8))+ #face="bold" negrita+
  annotation_custom(text_high,xmin=54,xmax=54,ymin=-0.15,ymax=-0.15)+
  coord_cartesian(clip = "off")


### light tolerance #####
speciesproportion <- read_csv("ST_vs_LD.csv")
colnames(speciesproportion)
speciesproportion$Plot=factor(speciesproportion$Plot,levels=c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2","M"))

colnames(speciesproportion)=c( "Plot"  ,    "mean"   ,   "se"    , "(b) Seed type", "Years"   ,  "LCL"   ,    "UCL"   ,    "Place"  )

speciesproportion$`(b) Seed type`[speciesproportion$`(b) Seed type`=="Shade_tolerant"]="Shade tolerant"
speciesproportion$`(b) Seed type`[speciesproportion$`(b) Seed type`=="Light_demanding"]="Light demanding"


speciesproportion2=speciesproportion[-which(speciesproportion$Plot=="M"),]
speciesproportion3=speciesproportion[which(speciesproportion$Plot=="M"),]
speciesproportion3$Years=54

model5 <- lm(mean~Years+Place,data=speciesproportion)
plot_model(model5,type='diag') # checking assumptions
summary(model5)


text_high <- textGrob("M", gp=gpar(fontsize=10)) #fontface="bold"

p<-  speciesproportion2 %>% ggplot(aes(Years,mean, group=`(b) Seed type`,color=`(b) Seed type`,fill=`(b) Seed type`)) +
  stat_smooth(method="lm",se=T)+
  geom_point(size=2) +
  scale_y_continuous(breaks =c(0.00,0.25,0.50,0.75,1.00))+
  scale_x_continuous(breaks =c(12,15,20,25,32,35,40,45))+
  ylab("Percentage of species")+
  xlab("Successional Age (y)")+
  expand_limits(x = 55)+
  geom_point(data = speciesproportion3,size=2)+
  annotate("text", x=51, y=0.985, label= paste0("R^2== 0.78"), size=3,parse=TRUE)+
  annotate("text", x=51, y=0.87, label= "p = 0.09", size=3)

N2=p+
  theme_bw()+
  theme(plot.margin=unit(c(20.5, 5.5, 16.5, 5.5),"points"),panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = NULL),axis.text.y = element_text(size=10),
        axis.text.x   = element_text(size=10),axis.title=element_text(size=10),strip.text.x = element_text(size = 10),legend.position = c(0.5,1.16),legend.direction = "horizontal",legend.spacing.y = unit(0.5, 'cm'),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.box.background=element_rect(fill='transparent', colour='transparent'),)+ #face="bold" negrita+
  annotation_custom(text_high,xmin=54,xmax=54,ymin=-0.21,ymax=-0.21)+
  coord_cartesian(clip = "off")



## Animal_non amimal#####
speciesproportion <- read_csv("Animal_vs_non_Animal.csv")

colnames(speciesproportion)
speciesproportion$Plot=factor(speciesproportion$Plot,levels=c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2","M"))
colnames(speciesproportion)=c( "Plot"  ,    "mean"   ,   "se"    ,    "(c) Seed type", "Years"   ,  "LCL"   ,    "UCL"   ,    "Place"  )

speciesproportion$`(c) Seed type`[speciesproportion$`(c) Seed type`=="Animal_dispersed"]="Animal dispersed"
speciesproportion$`(c) Seed type`[speciesproportion$`(c) Seed type`=="Non_Animal_dispersed"]="Non Animal dispersed"

speciesproportion2=speciesproportion[-which(speciesproportion$Plot=="M"),]
speciesproportion3=speciesproportion[which(speciesproportion$Plot=="M"),]
speciesproportion3$Years=54

model6 <- lm(mean~Years+Place,data=speciesproportion)
plot_model(model6,type='diag') # checking assumptions
summary(model6)


text_high <- textGrob("M", gp=gpar(fontsize=10)) #fontface="bold"

p<-  speciesproportion2 %>% ggplot(aes(Years,mean, group=`(c) Seed type`,color=`(c) Seed type`,fill=`(c) Seed type`)) +
  stat_smooth(method="lm",se=T)+
  geom_point(size=2) +
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(breaks =c(12,15,20,25,32,35,40,45))+
  ylab("Percentage of species")+
  xlab("Successional Age (y)")+
  expand_limits(x = 55)+
  geom_point(data = speciesproportion3,size=2)+
  annotate("text", x=51, y=0.98, label= paste0("R^2== 0.65"), size=3,parse=TRUE)+
  annotate("text", x=51, y=0.88, label= "p = 0.28", size=3)
#annotate("text", x=12, y=1.1, label= "(c)", size=5)
N3=p+
  theme_bw()+
  theme(plot.margin=unit(c(20.5, 5.5, 5.5, 5.5),"points"),panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = NULL),axis.text.y = element_text(size=10),
        axis.text.x   = element_text(size=10),axis.title=element_text(size=10),strip.text.x = element_text(size = 10),legend.position = c(0.5,1.145),legend.direction = "horizontal",legend.spacing.y = unit(0.5, 'cm'),legend.text=element_text(size=8),legend.title=element_text(size=8))+ #face="bold" negrita+
  annotation_custom(text_high,xmin=54,xmax=54,ymin=-0.15,ymax=-0.15)+
  coord_cartesian(clip = "off") 

multiplot(N1,N2,N3, cols=1)

tiff(file = "Fig7.tiff", width = 3.9, height = 6.5,units = "in", res = 600)
multiplot(N1,N2,N3, cols=1)
dev.off()

#####################################################
### Figure 6 - number of stems and number of seeds  
#####################################################

Wholedata <- read_csv("Wholedata_F.csv") # all seeds data
Todo_CR= Wholedata
Wholedata= Todo_CR
table(Wholedata$AllPosible)
table(Wholedata$LifeForm)
just_trees_seeds= Wholedata[-which(Wholedata$LifeForm=="L"),] # since Liana was not included in the tree monitoring 
table(just_trees_seeds$uniqueplace)

Wholedataby_plot_seed=just_trees_seeds%>% group_by(Data,uniqueplace,Code2) %>% 
  summarize(NSeeds=sum(AllPosible))
seed_just_trees=Wholedataby_plot_seed[,2:4]
colnames(seed_just_trees)
seed_just_trees=spread(seed_just_trees, Code2, NSeeds, fill=0)
table(seed_just_trees$uniqueplace)
seed_just_trees=as.data.frame(seed_just_trees)

##trees ___
whole_trees_matrix=read_csv("whole_trees_matrix.csv")
whole_trees_matrix2=whole_trees_matrix[,2:291]
whole_trees_matrix2$uniqueplace=as.factor(c("H_D","N_D","H_A","N_A","N_E","H_C","N_C","H_B","N_B"))
whole_trees_matrix2=whole_trees_matrix2[,2:291]
whole_trees_matrix2=as.data.frame(whole_trees_matrix2)

trees_list <- melt(whole_trees_matrix2, id="uniqueplace")
colnames(trees_list)
seed_list <- melt(seed_just_trees, id="uniqueplace")
colnames(seed_list)=c("LocalityCode", "variable"  ,   "Seed_density" )
colnames(trees_list)=c("LocalityCode", "variable"  ,   "value" )
seed_list=as.data.frame(seed_list)

# double_check_zero
trees_list2=trees_list[-which(trees_list$value==0),]
seed_list2=seed_list[-which(seed_list$Seed_density==0),]
table(seed_list2$LocalityCode)
table(trees_list2$LocalityCode)

abundance_whole = merge(trees_list2, seed_list2,all = TRUE)

abundance_whole$value[is.na(abundance_whole$value)] <- 0 # make zero where there is not tree for seed
abundance_whole$Seed_density[is.na(abundance_whole$Seed_density)] <- 0 # make zero where there is not seed for tree

abundance_whole$log_tree <- decostand(abundance_whole$value, "log")
abundance_whole$log_seed <- decostand(abundance_whole$Seed_density, "log")

table(abundance_whole$LocalityCode)
str(abundance_whole)
table(abundance_whole$LocalityCode)
abundance_whole$Years=NULL
abundance_whole$Years[(abundance_whole$LocalityCode== "H_A" )] = 12
abundance_whole$Years[(abundance_whole$LocalityCode== "H_B" )] = 15
abundance_whole$Years[(abundance_whole$LocalityCode== "H_C" )] = 20
abundance_whole$Years[(abundance_whole$LocalityCode== "H_D" )] = 25
abundance_whole$Years[(abundance_whole$LocalityCode== "N_A" )] = 33
abundance_whole$Years[(abundance_whole$LocalityCode== "N_B" )] = 36
abundance_whole$Years[(abundance_whole$LocalityCode== "N_C" )] = 41
abundance_whole$Years[(abundance_whole$LocalityCode== "N_D" )] = 46
abundance_whole$Years[(abundance_whole$LocalityCode== "N_E" )] = 100


abundance_whole$LocalityCode2=NULL
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "H_A" )] = "A1(12)"
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "H_B" )] = "B1(15)"
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "H_C" )] = "C1(20)"
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "H_D" )] = "D1(25)"
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "N_A" )] = "A2(33)"
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "N_B" )] = "B2(36)"
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "N_C" )] = "C2(42)" 
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "N_D" )] = "D2(46)" 
abundance_whole$LocalityCode2[(abundance_whole$LocalityCode== "N_E" )] = "E"

abundance_whole$LocalityCode2=factor(abundance_whole$LocalityCode2,levels = c("A1(12)", "B1(15)",   "C1(20)",  "D1(25)",  "A2(33)",  "B2(36)",  "C2(42)",  "D2(46)",  "E"))
table(abundance_whole$LocalityCode2)
abundance_whole$data=sapply(strsplit(as.character(abundance_whole$LocalityCode), ""), tail, 1)

#lm function.
df=abundance_whole
regression=function(df){
  #setting the regression function. 
  reg_fun<-lm(formula=df$log_seed~df$log_tree) #regression function
  #getting the slope, intercept, R square and adjusted R squared of 
  #the regression function (with 3 decimals).
  slope<-round(coef(reg_fun)[2],3)  
  #intercept<-round(coef(reg_fun)[1],3) 
  #R2<-round(as.numeric(summary(reg_fun)[8]),3)
  R2.Adj<-round(as.numeric(summary(reg_fun)[9]),3)
  c(slope,R2.Adj)
  #c(slope,intercept,R2,R2.Adj)
}
library(plyr)
regressions_data<-ddply(df,"LocalityCode2",regression) # LocalityCode2 is the plot site
colnames(regressions_data)<-c ("LocalityCode2","slope","intercept","R2","R2.Adj") # assumptions ok when log# seed and trees but no in row numberd

#abundancebytree=glm.nb(log_seed~log_tree*Years2, data =abundance_whole) # did not change significantly simpler lm approach so keep lm aproach
#summary(abundancebytree)

p <- abundance_whole %>% ggplot(aes(log_tree,log_seed)) +
  facet_grid(cols = vars(LocalityCode2)) +
  stat_smooth(method='lm', colour=c("black"),aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)))+
  scale_y_continuous(breaks = c(0, 2, 4, 6,8))+
  scale_x_continuous(breaks = c(0, 2, 4, 6))+
  ylab("Number of Seeds (log)")+
  xlab("Number of Trees (log)")
p2=p+
  geom_point(aes(log_tree,log_seed), colour = factor(abundance_whole$Time2),position = "jitter")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = NULL),axis.text.y = element_text(size=13),
        axis.text.x   = element_text(size=13),axis.title=element_text(size=15),strip.text.x = element_text(size = 12))+ #face="bold" negrita+
  geom_text(data = ann_text,label=ann_text$label, size=4) #+
  #geom_label(data=regressions_data, inherit.aes=FALSE,aes(x = 1.2, y = 32, label=paste("slope=",slope,","," ","
  #intercept=",intercept,","," ","R^2=",R2,","," ","R^2.Adj=",R2.Adj))


tiff(file = "Fig.6.tiff", width = 7, height = 4,units = "in", res = 600)
plot(p2)
dev.off()

#####################################################
### Figure 5 - inmigrants seeds  
#####################################################

Inmi= abundance_whole[abundance_whole$value==0,] # this is proportion of seed with no tree # using here raw numbers
Inmi2=Inmi[!Inmi$Years==100,]
Inmi_control_group= abundance_whole[abundance_whole$value>0,] # all seeds from at least 1 tree or non_immigrants
Inmi_control_group2=Inmi_control_group[!Inmi_control_group$Years==100,]

# number of seeds proportion ###
group1= aggregate(Inmi2[, 6], list(Inmi2$LocalityCode), sum) # inmigrants
group2= aggregate(Inmi_control_group2[, 6], list(Inmi_control_group2$LocalityCode), sum)
proportion=as.data.frame(group1$V1/(group2$V1+group1$V1))
proportion$Age=c(12,15,20,25,32,35,40,45)
proportion$unit=c("A","B","C","D","A","B","C","D")
colnames(proportion)=c("Proportion", "Age", "unit")

proportion2=as.data.frame(group1$V1/(group2$V1+group1$V1)) # proportion of inmi
proportion2$Age=c(45,40,35,32,25,20,15,12)
colnames(proportion2)=c("Proportion", "Age")

R=lmer(proportion$Proportion~proportion$Age+(1|unit))
R=lm(proportion$Proportion~proportion$Age)

plot_model(R,type='diag') # checking assumptions
summary(R)

Pr=proportion %>% ggplot(aes(Age,Proportion)) +
  stat_smooth(method='lm', colour=c("black"))+
  geom_point(size=2.5)+
  scale_y_continuous(limits=c(0,0.8))+
  scale_x_continuous(breaks = c(12, 15, 20, 25, 32, 35,40,45))+
  ylab("Proportion of immigrant seeds")+
  xlab("\nSuccessional Age (y)")+
  annotate("text", x=43, y=0.7, label= paste0("R^2 == 0.84"), size=4.5,parse=TRUE)+
  annotate("text", x=43, y=0.65, label= "p = < 0.01", size=4.5)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.y = element_text(size=20),axis.text.x   = element_text(size=20),
        axis.title=element_text(size=20))

tiff(file = "Immigrant_proportion2.tiff", width = 6, height = 5,units = "in", res = 150)
plot(Pr)
dev.off()






