# 2023-02-21 pulling all stats for McLaskey et al. 2023 paper together


library(tidyverse)
library(ggfortify)
library(vegan)
library(Hmisc)
library(corrplot)
library(cowplot)
library(pairwiseAdonis)


allData <- read.csv("processed_data/QU39 2015 zoop POM FA SI Chl biomass 20230126.csv") # keeping 5/11 250 net samples in
str(allData)
allData$Date <- as.Date(allData$Date, "%Y-%m-%d")
allData <- allData %>% filter(Date < as.Date("2015-12-31"))

allData$Month <- format(as.Date(allData$Date), "%m")
allData$Month <- as.factor(allData$Month)
allData <- allData %>% mutate(season = ifelse( Month %in% c("03","04","05"), "spring" ,
                                                             ifelse(Month %in% c("06", "07", "08"), "summer", 
                                                                    "fall")))
allData$season <- factor(allData$season, levels = c("spring", "summer","fall"))

allData$Size.Fraction <- as.factor(allData$Size.Fraction)
levels(allData$Size.Fraction)  # shows the different factor levels
allData$Size.Fraction <- factor(allData$Size.Fraction, levels = c("2000", "1000","500", "250", "125","64", "POM"))



# Modify POM dates to match zoop collection dates ------------------------------------------------------

# To make comparisons between POM parameters and zooplankton, 
# I need to modify POM dates to match the closest zooplankton collection
# And I need to manually change some POM dates to be able to match up w zoop collections

# All are within two days, except the last sampling 
allData$Date <- as.character(allData$Date)

allData <- allData %>%  mutate(Date2 = ifelse(Date=="2015-04-14", "2015-04-15",
                                                 ifelse(Date=="2015-05-26", "2015-05-27",
                                                    ifelse(Date=="2015-07-07", "2015-07-06",
                                                        ifelse(Date=="2015-08-12", "2015-08-10",
                                                           ifelse(Date=="2015-09-09", "2015-09-10",
                                                              ifelse(Date=="2015-10-06", "2015-10-05",
                                                                  ifelse(Date=="2015-11-03", "2015-11-02",
                                                                     ifelse(Date=="2015-12-11", "2015-11-30", # this is the greatest different between collections 
                                                                         ifelse(Date=="2015-03-18", "2015-03-17",
                                                                            ifelse(Date=="2015-06-15", "2015-06-16",
                                                                               ifelse(Date=="2015-06-24", "2015-06-22",
                                                                                   ifelse(Date=="2015-05-13", "2015-05-11",
                                                                                      (Date))))))))))))))


allData$Date <- as.Date(allData$Date, "%Y-%m-%d")

# two different versions of the dataset I have used
allData$Prop.ID <- allData$SumFA / allData$Total.FA_G_AREA

fatty.acid.all <- allData[!is.na(allData$C16.0_PERCENT),]
QU39.2015.64 <- fatty.acid.all[fatty.acid.all$Size.Fraction!="POM",]






ggplot(allData, aes(x=Total.FA_G_AREA, y=SumFA)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0)


mean(QU39.2015.64$Prop.ID)
# 0.8553788
# Only for the zooplankton, POMFA haven't kept the AREAs 





# Fig 1 - Chl and Zoop Biomass ---------------------------------------------


chlSumm.props.small <- allData %>% select(Date, chla_prop_GF.F, chla_prop_20um, chla_prop_3um)

# Chl size class proportions
# chlSumm.props.small <- Chl.2015 %>% select(Date, chla_prop_GF.F, chla_prop_20um, chla_prop_3um)
chlSumm.props.small.long <- chlSumm.props.small %>% pivot_longer(-Date, names_to = "Size.class", values_to = "Prop")
chlSumm.props.small.long$Size.class <- as.factor(chlSumm.props.small.long$Size.class)
ord_size_class <- c( "chla_prop_GF.F",  "chla_prop_3um", "chla_prop_20um")
chlSumm.props.small.long<- chlSumm.props.small.long %>%  mutate(Size.class = factor(Size.class, 
                                                                                    levels = ord_size_class))

pChl.allProps <- ggplot(chlSumm.props.small.long, aes(Date, Prop, fill=Size.class)) + geom_area() + 
  theme_bw()+ 
  scale_fill_manual(values=c("#c2e699","#78c679", "#006837")) + 
  #theme(legend.position="none") +
  labs(y=element_blank(), fill="Size Class")   + 
  theme(legend.position = "none",
        panel.grid.minor=element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 0.3)) + 
  scale_x_date(limits=c(as.Date("2015-01-01"), as.Date("2015-12-31")),
               expand=c(0,0)) +  
  scale_y_continuous(expand=c(0,0), labels = c("0", "25%", "50%", "75%", ""))  + 
  geom_text(label="(d)", x=as.Date("2014-11-25"), y=0.98, color="black") + 
  coord_cartesian(clip = "off") 
pChl.allProps



# Chl size class concentrations
chlSumm.concs.small <- Chl.2015 %>% select(Date, chl_GF.F, chl_20um, chl_3um)
chlSumm.concs.small.long <- chlSumm.concs.small %>% pivot_longer(-Date, names_to = "Size.class", values_to = "Conc")
chlSumm.concs.small.long$Size.class <- as.factor(chlSumm.concs.small.long$Size.class)
ord_size_class <- c( "chl_GF.F",  "chl_3um", "chl_20um")
chlSumm.concs.small.long<- chlSumm.concs.small.long %>%  mutate(Size.class = factor(Size.class, 
                                                                                    levels = ord_size_class))

pChl.allconcs <- ggplot(chlSumm.concs.small.long, aes(Date, Conc, fill=Size.class)) + geom_area() + 
  theme_bw()+ 
  scale_fill_manual(values=c("#c2e699","#78c679", "#006837")) + 
  labs(x=element_blank(), y=expression(paste("Chl a (",mu,"g L"^-1*")")), fill="Size Class")  + 
  theme(panel.grid.minor=element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.3)) + 
  scale_x_date(limits=c(as.Date("2015-01-01"), as.Date("2015-12-31")),
               expand=c(0,0)) +  
  scale_y_continuous(expand=c(0,0)) + 
  geom_text(label="(c)", x=as.Date("2014-11-25"), y=19.5, color="black") + 
  coord_cartesian(clip = "off") 
pChl.allconcs



# plot for Chl legend 
pnew <- ggplot(chlSumm.concs.small.long, aes(Date, Conc, fill=Size.class)) + geom_area() + 
  scale_fill_manual(values=c("#c2e699","#78c679", "#006837"), 
                    labels = c("0.7", "3", "20")) + 
  labs(fill="Size")  
pnew

legend.Chl <- get_legend(
  # create some space to the left of the legend
  pnew + theme(legend.title = element_text(size=9),
               legend.text = element_text(size = 7),
               legend.position = c(0.45,0.6))
)



# Zooplankton Biomass

zoop.biomass2015 <- allData %>% select(Date, Size.Fraction, Biomass.mg.m3)
zoop.biomass2015 <- zoop.biomass2015[!is.na(zoop.biomass2015$Biomass.mg.m3),]
zoop.biomass2015$Size.Fraction <- factor(zoop.biomass2015$Size.Fraction, levels = c("64", "125","250", "500", "1000","2000" ))

pZoop.biomass <- ggplot(zoop.biomass2015, aes(x=Date, y=Biomass.mg.m3, fill=Size.Fraction)) + geom_area() +
  theme_bw() + 
  labs(y=expression(paste("Biomass (mg m"^-3*")")), 
       x=element_blank()) + 
  theme(panel.grid.minor=element_blank(),
        panel.border = element_rect(linewidth = 0.3) ) + 
  scale_fill_manual(values =  c("#47ebb4", "#66ccff","#9966ff",
                                "#db5764",  "#ff7433",  "#ffbf00")) + 
  scale_x_date(limits=c(as.Date("2015-01-01"), as.Date("2015-12-31")),
               expand=c(0,0)) +  
  scale_y_continuous(expand=c(0,0)) + 
  geom_text(label="(e)", x=as.Date("2014-11-25"), y=66, color="black") + 
  coord_cartesian(clip = "off")

pZoop.biomass


pZoop.props <- ggplot(zoop.biomass2015, aes(x=Date, y=Biomass.mg.m3, fill=Size.Fraction)) + geom_area(position = "fill") +
  labs(y=element_blank(), x=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.minor=element_blank(),
        legend.position = "none" ,
        panel.border = element_rect(linewidth = 0.3)) + 
  scale_fill_manual(values =  c("#47ebb4", "#66ccff","#9966ff",
                                "#db5764",  "#ff7433",  "#ffbf00")) + 
  scale_x_date(limits=c(as.Date("2015-01-01"), as.Date("2015-12-31")),
               expand=c(0,0)) +  
  scale_y_continuous(expand=c(0,0), labels = c("0", "25%", "50%", "75%", ""))  + 
  geom_text(label="(f)", x=as.Date("2014-11-25"), y=0.98, color="black") + 
  coord_cartesian(clip = "off") 

pZoop.props


# plot for legend 
pnew2 <- ggplot(zoop.biomass2015, aes(x=Date, y=Biomass.mg.m3, fill=Size.Fraction)) + geom_area(position = "fill") +
  theme_bw() +  
  scale_fill_manual(values =  c("#47ebb4", "#66ccff","#9966ff",
                                "#db5764",  "#ff7433",  "#ffbf00")) + 
  labs(fill= "Size") 
pnew2

legend.Zoop <- get_legend(
  # create some space to the left of the legend
  pnew2 + theme(legend.title = element_text(size=9),
                legend.text = element_text(size = 7),
                legend.position = c(0.45,0.6))
)


png("Biomass time series Fig 2 20221212.png", width=190, height=110, units="mm", res=300)
plot_grid(pChl.allconcs, pChl.allProps, legend.Chl, 
          pZoop.biomass, pZoop.props, legend.Zoop, 
          rel_widths = c(1.5, 1.5, 0.3), nrow=2, 
          rel_heights = c(0.9, 1))
dev.off()


tiff("Biomass time series Fig 2 20220113.tiff", width=190, height=110, units="mm", res=300)
plot_grid(pChl.allconcs, pChl.allProps, legend.Chl, 
          pZoop.biomass, pZoop.props, legend.Zoop, 
          rel_widths = c(1.5, 1.5, 0.3), nrow=2, 
          rel_heights = c(0.9, 1))
dev.off()




# I also need POM TFA time series the same size as the Chl panel 
# Plot the POM FA time series separately to overlay in Illustrator 

QU39.2015.POM <- fatty.acid.all[fatty.acid.all$Size.Fraction=="POM",]

p1 <- ggplot(allData, aes(x=Date, y=SumFA_ug.L, color="#734f22")) + 
  geom_point() + theme_bw() + 
  scale_x_date( limits = as.Date(c("2015-01-01", "2015-12-31")),
                expand = c(0,0)) + 
  scale_y_continuous(position="right",
                     expand = c(0,0), 
                     limits = c(0,71)) + 
  theme(legend.position = "none",
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text.x = element_blank(), 
        panel.border = element_blank(),
        axis.text.y = element_text(size = 14)) + 
  labs(y = "", x="") + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))

p1


ggsave(
  plot = p1,
  filename = "POM FA time series 3 no smooth.png",
  bg = "transparent"
)





# Trophic position --------------------------------------------------------


# Need to calculate 14 day running average of POM Del15N 
allData.sm <- allData %>% select(Date, Size.Fraction, delta15n)
allData.sm <- allData.sm[complete.cases(allData.sm$delta15n),]
allData.sm.agg =  aggregate(allData.sm[,3], 
                            by= list(allData.sm$Date, allData.sm$Size.Fraction), 
                            FUN = mean, na.omit=TRUE)

colnames(allData.sm.agg) <- c("Date", "Size.Fraction",  "delta15n")
allData.sm.wide <- allData.sm.agg %>% pivot_wider(values_from = delta15n, names_from = Size.Fraction)

# Going to use JDay to calculate window rather than dates 
allData.sm.wide$JDay <- format(as.Date(allData.sm.wide$Date), "%j")
allData.sm.wide$JDay <- as.numeric(allData.sm.wide$JDay)


# CAlculate 14 day moving average of POM DelN15
allData.sm.wide$meanPOM7 <- NA
allData.sm.wide$meanPOM14 <- NA
allData.sm.wide$meanPOM21 <- NA
allData.sm.wide$meanPOM28 <- NA
for(i in 1:nrow(allData.sm.wide)){
  allData.sm.wide$meanPOM7[i] <- mean(allData.sm.wide$POM[allData.sm.wide$JDay >= (allData.sm.wide$JDay[i]-7) & 
                                                            allData.sm.wide$JDay<=allData.sm.wide$JDay[i] ], na.rm = T)
  allData.sm.wide$meanPOM14[i] <- mean(allData.sm.wide$POM[allData.sm.wide$JDay >= (allData.sm.wide$JDay[i]-14) & 
                                                             allData.sm.wide$JDay<=allData.sm.wide$JDay[i] ], na.rm = T)
  allData.sm.wide$meanPOM21[i] <- mean(allData.sm.wide$POM[allData.sm.wide$JDay >= (allData.sm.wide$JDay[i]-21) & 
                                                             allData.sm.wide$JDay<=allData.sm.wide$JDay[i] ], na.rm = T)
  allData.sm.wide$meanPOM28[i] <- mean(allData.sm.wide$POM[allData.sm.wide$JDay >= (allData.sm.wide$JDay[i]-28) & 
                                                             allData.sm.wide$JDay<=allData.sm.wide$JDay[i] ], na.rm = T)
}

# Calculate trophic position after El-Sabaawi et al. 2009
allData.sm.wide$TP.2000 <- ((allData.sm.wide$`2000` - allData.sm.wide$meanPOM14) / 3.5) +1
allData.sm.wide$TP.1000 <- ((allData.sm.wide$`1000` - allData.sm.wide$meanPOM14)/ 3.5) +1
allData.sm.wide$TP.500 <- ((allData.sm.wide$`500` - allData.sm.wide$meanPOM14)/ 3.5) +1
allData.sm.wide$TP.250 <- ((allData.sm.wide$`250` - allData.sm.wide$meanPOM14)/ 3.5) +1
allData.sm.wide$TP.125 <- ((allData.sm.wide$`125` - allData.sm.wide$meanPOM14)/ 3.5) +1
allData.sm.wide$TP.64 <- ((allData.sm.wide$`64` - allData.sm.wide$meanPOM14)/ 3.5) +1


allData.sm.long <- allData.sm.wide %>% pivot_longer(TP.2000:TP.64, names_to = "Size", values_to = "TP")

allData.sm.long$Size <- factor(allData.sm.long$Size, levels = c("TP.2000", "TP.1000","TP.500", "TP.250", "TP.125","TP.64"))




allData.sm.wide.smaller <- allData.sm.wide %>% select(Date, TP.2000:TP.64)
allData.smaller.long <- allData.sm.wide.smaller %>% pivot_longer(TP.2000:TP.64, names_to = "Size", values_to = "TP")
allData.smaller.long$Size <- substring(allData.smaller.long$Size, 4,7)
colnames(allData.smaller.long)[2] <- "Size.Fraction"
QU39.2015.64 <- full_join(allData.smaller.long, QU39.2015.64)





# *Corr plot TP vs FATM ---------------------------------------------------------------

QU39.2015.64$SFA.PUFA <- QU39.2015.64$percent.SFA / QU39.2015.64$percent.PUFA

QU39.2015.64.smer <- QU39.2015.64 %>% select(Date, Size.Fraction, SFA.PUFA, Carn18.1N9_n7, DHA.EPA) 
QU39.2015.64.agg = aggregate(QU39.2015.64.smer[,c(3:5)],
                             by = list(QU39.2015.64.smer$Date, QU39.2015.64.smer$Size.Fraction),
                             FUN = mean, na.rm=TRUE)

colnames(QU39.2015.64.agg)[1:2] <- c("Date", "Size.Fraction")
QU39.2015.64.agg <- QU39.2015.64.agg[!is.nan(QU39.2015.64.agg$SFA.PUFA),]

QU39.2015.64.wide <- QU39.2015.64.agg %>% pivot_wider(names_from = Size.Fraction, values_from = c(SFA.PUFA, Carn18.1N9_n7, 
                                                                                                  DHA.EPA))
colnames(QU39.2015.64.wide)[2:7] <- paste0("FATM_", colnames(QU39.2015.64.wide)[2:7])

QU39.2015.64.smer2 <- QU39.2015.64 %>% select(Date, Size.Fraction, TP) 
QU39.2015.64.agg2 = aggregate(QU39.2015.64.smer2[,c(3)],
                              by = list(QU39.2015.64.smer2$Date, QU39.2015.64.smer2$Size.Fraction),
                              FUN = mean, na.rm=TRUE)
colnames(QU39.2015.64.agg2) <- c("Date", "Size.Fraction", "TP")

QU39.2015.64.agg2 <- QU39.2015.64.agg2[!is.nan(QU39.2015.64.agg2$TP),]
QU39.2015.64.wide2 <- QU39.2015.64.agg2 %>% pivot_wider(names_from = Size.Fraction, values_from = TP)
colnames(QU39.2015.64.wide2)[2:7] <- paste0("TP_", colnames(QU39.2015.64.wide2)[2:7])


QU39.2015.64.wide.all <- full_join(QU39.2015.64.wide, QU39.2015.64.wide2)


resTP <- rcorr(as.matrix(QU39.2015.64.wide.all[2:13]), type = "spearman")

corrplot(resTP$r, type="lower", method="number",
         p.mat = resTP$P, sig.level = 0.05, insig = "blank") 




# del13C correlation ------------------------------------------------------

# aggregate the multiple zooplankton nets (two 64 um nets done on same day) before I can pivot

allSIdatasm2 <- allData %>% select(Date2, Size.Fraction, delta13c) 
allSIdatasm2 <- allSIdatasm2[!is.na(allSIdatasm2$delta13c),]

SI.agg = aggregate(allSIdatasm2[,c(3)],
                   by = list(allSIdatasm2$Date2, allSIdatasm2$Size.Fraction),
                   FUN = mean, na.omit=TRUE)
colnames(SI.agg) <- c("Date", "Size.Fraction", "delta13c")

SI.agg.wide <- SI.agg %>% pivot_wider(names_from = Size.Fraction, values_from = delta13c)


res2<-rcorr(as.matrix(SI.agg.wide[2:8]), type = "pearson")
corrplot(res2$r, type="full", method="number",
         p.mat = res2$P, sig.level = 0.1,  tl.col = "black",
         title = "delta13c")





# FATM correlations -------------------------------------------------------

# All FATM to test
# C14.0_PERCENT,   C16.0_PERCENT,  C16.1n.7_PERCENT, C18.0_PERCENT,  C16.3n.4_PERCENT,  
# C18.1n.7_PERCENT,  C18.1n.9c_PERCENT,  C18.4n.3_PERCENT,  C20.5n.3_PERCENT
# C22.6n.3_PERCENT,  C16.2n.4_PERCENT,  Bacteria_15_17


fatty.acid.smer <- fatty.acid.all %>% select(Date2, Size.Fraction, Bacteria_15_17) 
# arcsine transform FA props
fatty.acid.smer[,3] <- asinTransform(fatty.acid.smer[,3])


fatty.acid.agg = aggregate(fatty.acid.smer[,c(3)],
                           by = list(fatty.acid.all$Date2, fatty.acid.all$Size.Fraction),
                           FUN = mean, na.rm=TRUE)
colnames(fatty.acid.agg) <- c("Date", "Size.Fraction", "FA")

fatty.acid.wide <- fatty.acid.agg %>% pivot_wider(names_from = Size.Fraction, values_from = FA)


res2 <- rcorr(as.matrix(fatty.acid.wide[2:8]), type = "spearman")

corrplot(res2$r, type="lower", method="number",
         p.mat = res2$P, sig.level = 0.1, insig = "blank") 




# Fig 4 formatted corr tables ----------------------------------------------------


TPcorrtable <- read.csv("processed_data/trophic position correlation table.csv", stringsAsFactors = FALSE)
str(TPcorrtable)

# columns with the R2 values to base color shading on
TPcorrtable_r <- TPcorrtable %>% select(FATM, r.64:r.2000)
TPcorrtable_r_long <- TPcorrtable_r %>% pivot_longer(r.64:r.2000, names_to = "Size", values_to = "R")

# columns with what I want printed in the cells
TPcorrtable_label <- TPcorrtable %>% select(FATM, z.64:z.2000)
# Make colnames match table of R values so I can bind them later
colnames(TPcorrtable_label)[2:7] <- colnames(TPcorrtable_r)[2:7]
TPcorrtable_label_long <- TPcorrtable_label %>% pivot_longer(r.64:r.2000, names_to = "Size", 
                                                             values_to = "Label",
                                                             values_transform = as.character)



TPcorrtable_all_long <- full_join(TPcorrtable_r_long, TPcorrtable_label_long)



# column names for variables I want on horizontal side
TPcorrtable_all_long$FATM
mylevels1 <- c("SFA/PUFA", "18:1n-9 / 18:1n-7", "DHA:EPA")
mylevels2 <- c("r.64", "r.125" , "r.250","r.500","r.1000", "r.2000")

TPcorrtable_all_long <- TPcorrtable_all_long %>% filter(FATM %in% mylevels1)
# reorder the variables to how I want them
TPcorrtable_all_long$FATM <- ordered(TPcorrtable_all_long$FATM, levels = mylevels1)
TPcorrtable_all_long$Size <- ordered(TPcorrtable_all_long$Size, levels = mylevels2)



# label names I actually want to use for variables on horizontal side
labels1 <- c(expression(paste("64 ",mu,"m")), 
             expression(paste("125 ",mu,"m")), 
             expression(paste("250 ",mu,"m")), 
             expression(paste("500 ",mu,"m")),
             expression(paste("1000 ",mu,"m")),
             expression(paste("2000 ",mu,"m")))

# label names I actually want to use for variables on horizontal side
labels3 <- c("DHA:EPA", expression(paste("18:1",omega,"9 / 18:1",omega,"7")),"SFA / PUFA" )
#labels2 <- c("SFA/PUFA" , expression(paste("18:1",omega,"9 / 18:1",omega,"7")),"DHA:EPA")


mycolors <- c( "#47ebb4", "#66ccff", "#9966ff",
               "#db5764","#ff7433", "#ffbf00")
# slightly modified
mycolors <- c( "#3bc496", "#51a5cf", "#9966ff",
               "#db5764","#ff7433", "#ffbf00")
mycolors <- c( "#3bc496", "#51a5cf", "#9966ff",
               "#db5764","#ff7433", "#e8aa00")


# Making the plot

pCorr2 <- TPcorrtable_all_long %>% 
  ggplot(aes(Size, FATM, fill=R, label=Label), linewidth=2) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation") + 
  scale_fill_gradient2(mid="#ffffff",low=paste("red"),high=paste("blue"), limits=c(-1,1), na.value="ivory2") +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0), position = "top", labels=(labels1)) +
  scale_y_discrete(expand=c(0,0),
                   labels=rev(labels3)) +
  theme(legend.position = "none", 
        axis.text.y = element_text(hjust = 0, size= 11),
        axis.text.x = element_text(colour = mycolors, size= 13)) + 
  geom_text(label="(b)", x=-0.8, y=3.69, color="black") + 
  coord_cartesian(clip = "off")



#png("TP correlations Fig 6 20221202.png", width=150, height=45, units="mm", res=300)
pCorr2
#dev.off()







# POM FATM correlation table


POMcorrtable <- read.csv("processed_data/POM FATM correlation table.csv", stringsAsFactors = FALSE)
str(POMcorrtable)
POMcorrtable$z.1000 <- as.character(POMcorrtable$z.1000)

# columns with the R2 values to base color shading on
POMcorrtable_r <- POMcorrtable %>% select(FATM, r.64:r.2000)
POMcorrtable_r_long <- POMcorrtable_r %>% pivot_longer(r.64:r.2000, names_to = "Size", values_to = "R")

# columns with what I want printed in the cells
POMcorrtable_label <- POMcorrtable %>% select(FATM, z.64:z.2000)
# Make colnames match table of R values so I can bind them later
colnames(POMcorrtable_label)[2:7] <- colnames(POMcorrtable_r)[2:7]
POMcorrtable_label_long <- POMcorrtable_label %>% pivot_longer(r.64:r.2000, names_to = "Size", 
                                                               values_to = "Label",
                                                               values_transform = as.character)



POMcorrtable_all_long <- full_join(POMcorrtable_r_long, POMcorrtable_label_long)



# column names for variables I want on horizontal side
POMcorrtable_all_long$FATM
mylevels1 <- c( "18:1w7", "DHA",  "EPA", "C16 PUFAs ", "SDA (18:4w4)", "ALA (18:3w3)","?13C ") 

mylevels2 <- c("r.64", "r.125" , "r.250","r.500","r.1000", "r.2000")

# POMcorrtable_all_long <- POMcorrtable_all_long %>% filter(FATM %in% mylevels1)
# reorder the variables to how I want them
POMcorrtable_all_long$FATM <- ordered(POMcorrtable_all_long$FATM, levels = mylevels1)
POMcorrtable_all_long$Size <- ordered(POMcorrtable_all_long$Size, levels = mylevels2)


# label names I actually want to use for variables on horizontal side
labels1 <- c(expression(paste("64 ",mu,"m")), 
             expression(paste("125 ",mu,"m")), 
             expression(paste("250 ",mu,"m")), 
             expression(paste("500 ",mu,"m")),
             expression(paste("1000 ",mu,"m")),
             expression(paste("2000 ",mu,"m")))


# label names I actually want to use for variables on horizontal side
labels2 <- c(expression(paste(delta^{13}, "C")), 
             expression(paste("ALA (18:3",omega,"3) ", italic("chlorophytes"))),
             expression(paste("SDA (18:3",omega,"3)", italic("cryptophytes"))),
             expression(paste("C16 PUFAs\n", italic("diatoms"))),
             expression(paste("EPA\n", italic("diatoms"))),
             expression(paste("DHA\n", italic("dinoflagelates"))),
             "18:1\u03C97\nchlorophytes")

labels2 <- c(expression(paste(delta^{13}, "C")), 
             "ALA (18:3\u03C93)\nchlorophytes",
             "SDA (18:3\u03C94)\ncryptophytes",
             "C16 PUFAs\ndiatoms",
             "EPA\ndiatoms",
             "DHA\ndinoflagelates",
             "18:1\u03C97\nchlorophytes")

# Making the plot

mycolors <- c( "#47ebb4", "#66ccff", "#9966ff",
               "#db5764","#ff7433", "#ffbf00")
# slightly modified
mycolors <- c( "#3bc496", "#51a5cf", "#9966ff",
               "#db5764","#ff7433", "#e8aa00")

pCorrPOM <- POMcorrtable_all_long %>% 
  ggplot(aes(Size, FATM, fill=R, label=Label), linewidth=2) +
  geom_tile() + 
  labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation") + 
  scale_fill_gradient2(mid="#ffffff",low=paste("red"),high=paste("blue"), limits=c(-1,1), na.value="ivory2") +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0), position = "top", labels=(labels1)) +
  scale_y_discrete(expand=c(0,0),
                   labels=rev(labels2), ) +
  theme(legend.position = "none", 
        axis.text.y = element_text(hjust = 0, vjust = 1,
                                   colour = "#089908", size= 12, face="bold"),
        axis.text.x = element_text(colour = mycolors, size= 13, face="bold")) + 
  geom_text(label="(a)", x=-0.78, y=7.8, color="black") + 
  coord_cartesian(clip = "off")



#png("POM correlations Fig 5 20221202.png", width=150, height=100, units="mm", res=300)
pCorrPOM
#dev.off()



# Make combined version w both as two panels, with a legend for fill color 

pCorrLegend <- POMcorrtable_all_long %>% 
  ggplot(aes(Size, FATM, fill=R, label=Label), linewidth=2) +
  scale_fill_gradient2(mid="#ffffff",low=paste("red"),high=paste("blue"), limits=c(-1,1), na.value="ivory2") +
  geom_tile() + 
  labs(fill="correlation\ncoefficient")  

# extract the legend from the montly PCA
corrLegend <- get_legend(
  # create some space to the left of the legend
  pCorrLegend + theme(legend.title = element_text(size=9),
                      legend.text = element_text(size = 7),
                      legend.position = c(0.5, 0.3)) 
)


png("Combined POM correlations Fig 4 20230201.png", width=180, height=130, units="mm", res=300)
plot_grid(pCorrPOM,corrLegend,pCorr2  ,
          rel_widths = c(1.5, 0.2), nrow=2, 
          rel_heights = c(2, 1),
          align = "v", 
          axis = "l")
dev.off()








# *Fig 2 ------------------------------------------------------------------



ptrophicP <- ggplot(allData.sm.long, aes(x=Date, y=TP, color=Size, fill=Size)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y="Trophic Position", x="", color="Size class" , fill="Size class") + 
  scale_x_date(expand = c(0,0),
               breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01", "2016-01-01")),
               labels = c("Jan", "Apr", "Jul", "Oct", "Jan"),
               limits = as.Date(c("2015-01-01", "2015-12-31"))) + 
  theme(legend.position = "none", 
        plot.margin = unit(c(0.1, 0.2, 0, 0.2), "cm"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=7), 
        axis.title = element_text(size=9)) +
  geom_text(label="(b)", x=as.Date("2014-11-13"), y=3, color="black") + 
  coord_cartesian(clip = "off")
ptrophicP


pDel13C <- ggplot(allData, aes(x=Date, y=delta13c, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste(delta^{13}, "C")), x="", color="" , fill="", size=7) + 
  scale_x_date(expand = c(0,0),
               limits = as.Date(c("2015-01-01", "2015-12-31"))) + 
  theme(legend.position = c(0.32,0.3), legend.text = element_text(size = 6), 
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        plot.margin = unit(c(0.1, 0, -0.3, 0.2), "cm"),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text = element_text(size=7), 
        axis.title = element_text(size=9), 
        legend.key.height = unit(0.35, 'cm')) +
  geom_text(label="(a)", x=as.Date("2014-11-13"), y=-18, color="black") + 
  coord_cartesian(clip = "off")

pDel13C



prow <- plot_grid(pDel13C, ptrophicP, 
                  align = 'v',
                  nrow = 2,
                  rel_heights = c(0.9,1))

prow




png("20230201 carbon and TP time series.png", width=90, height=125, units="mm", res=300)
prow
dev.off()





# abundant FAs in Zoops ------------------------------------------------------------


# Galloway and Winder 2015 arcsine-square root transformed the FA proportions
asinTransform <- function(p) { asin(sqrt(p)) }

# Make reduced data matrix of Zoop FAs
QU39.2015.64 <- QU39.2015.64[!is.na(QU39.2015.64$C22.6n.3_PERCENT),]
Q20.matrix <- QU39.2015.64 %>%  select(C14.0_PERCENT:C22.1n.11_PERCENT, Bacteria_15_17)



# What are the average contributions of each FA
contributions <- data.frame(matrix(nrow=ncol(Q20.matrix), ncol=3))
for(i in 1:ncol(Q20.matrix)){
  contributions[i,1] <- names(Q20.matrix[i])
  contributions[i,2] <- mean(Q20.matrix[,i])*100
  temp <- (Q20.matrix[i]==0)
  contributions[i,3] <- length(temp[temp==TRUE])
}
names(contributions) <- c("FA", "Mean.Percentage", "Num zeros")


# separate out peaks that are >1%
abundant.all <- contributions$FA[contributions$Mean.Percentage>=1]
abundant.all <- abundant.all[!is.na(abundant.all)]
(abundant.all <- abundant.all[-c(8)])


# remove FATM for zooplankton and 18:2n-6 because it coeluted with 16:4n-1
abundant <- abundant.all[-c(9,14,12)]
Q20.matrix.abund <- Q20.matrix[,abundant]

colnames(Q20.matrix.abund) <- paste(str_remove(colnames(Q20.matrix.abund), "_PERCENT"), "", sep = "")



# Arcsine square root Transform data
Q20.matrix.abund.transformed <- asinTransform(Q20.matrix.abund)




# FA summary table --------------------------------------------------------


QU39.2015.64.sm <- QU39.2015.64[!is.na(QU39.2015.64$C16.0_PERCENT),]
# QU39.2015.64.bySize <- QU39.2015.64.sm %>% select(Size.Fraction, C14.0_PERCENT:C22.1n.11_PERCENT, C12.0_PERCENT:Diatom.II_PERCENT)
QU39.2015.64.bySize <- QU39.2015.64.sm %>% select(Size.Fraction, abundant.all, Prop.ID)
QU39.2015.64.bySize <- QU39.2015.64.bySize[complete.cases(QU39.2015.64.bySize),]
QU39.2015.64.bySize <- QU39.2015.64.bySize %>% mutate(SumFA=rowSums(.[c(2:16)]))

QU39.2015.64.bySize.grouped <- group_by(QU39.2015.64.bySize, Size.Fraction)    # create an internal grouping structure

# Summary functions
summ.all.mean <- summarise_all(QU39.2015.64.bySize.grouped, mean)
summ.all.sd <- summarise_all(QU39.2015.64.bySize.grouped, sd)

write.csv(summ.all.sd, "20230206 FA contributions by size class sd.csv")
write.csv(QU39.2015.64.bySize, "20230206 FA table.csv")

# n for each size class
table(QU39.2015.64.sm$Size.Fraction)




# NMDS --------------------------------------------------------------------

Size <- QU39.2015.64$Size.Fraction
Months <- QU39.2015.64$Month
Seasons <- QU39.2015.64$season
taxa.names <- QU39.2015.64$Size.Fraction
taxa.names <- droplevels(taxa.names)

#create matrix for NMDS calculation:
species.matrix.sm <- as.matrix(Q20.matrix.abund)
#change diet dataframe into a matrix
class(species.matrix.sm) <- "numeric"
#make sure your numbers are treated as numbers
proportions_matrix <- decostand(species.matrix.sm, "total")
#calculations proportional biomass for each fish stomach
#total = 1 so it's expressed as a decimal. It is NOT total = 100 and a percentage.
transformed_matrix <- asin(sqrt(proportions_matrix))
#arc sine square root transformation of diet data


set.seed(505)
eco.nmds.bc <- metaMDS(transformed_matrix, distance="bray",labels=Size, trymax = 100, autotransform = FALSE)
# eco.nmds.bc <- metaMDS(transformed_matrix, distance="bray",labels=Size, trymax = 100, autotransform = FALSE, k=3)
eco.nmds.bc[1]

stressplot(eco.nmds.bc)

plot(eco.nmds.bc,  type = "p")


species.scores <- as.data.frame(scores(eco.nmds.bc, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- paste(str_remove(row.names(species.scores), "_PERCENT"), "", sep = "")

species.scores  #look at the data

# NMDS1: 22:1n-11, 15:0+17:0, 16:3n-4, 18:0
# NMDS2: 18:1n-9, 14:0, 16:3n-4, 16:2n-4


sites.scores <- as.data.frame(scores(eco.nmds.bc, "sites"))  
sites.scores <- cbind(sites.scores, QU39.2015.64$Size.Fraction, QU39.2015.64$Date)
colnames(sites.scores)[3:4] <- c("Size.Fraction", "Date")




# *PERMANOVA ---------------------------------------------------------------


permanova_eco.bc<-adonis2(transformed_matrix ~ taxa.names*Seasons, permutations = 999, method="bray")
permanova_eco.bc<-adonis2(transformed_matrix ~ taxa.names, permutations = 999, method="bray")
permanova_eco.bc #if significant, then plot it
# significant w and wo POM

pairwise.adonis(transformed_matrix, taxa.names)


eco.nmds.bc <- metaMDS(transformed_matrix, distance="bray",labels=Size, trymax = 100, autotransform = FALSE)

NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2],group=taxa.names, month= Months, season=Seasons)
#dataframe for plotting NMDS

eco.nmds.bc$species



# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Ellipses are standard deviation, no scaling of data (can use standard error, scaling, and confidence limit options)

ord.bc<-ordiellipse(eco.nmds.bc,taxa.names,display="sites",kind="sd", conf = 0.95, label=T)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


df_ell.bc <- data.frame()
for(g in levels(NMDS.bc$group)){
  df_ell.bc <- rbind(df_ell.bc, cbind(as.data.frame(with(NMDS.bc[NMDS.bc$group==g,],
                                                         veganCovEllipse(ord.bc[[g]]$cov,ord.bc[[g]]$center))),group=g))
}
#https://www.rpubs.com/RGrieger/545184


df_ell.bc$group <- as.factor(df_ell.bc$group) 
df_ell.bc$group <- factor(df_ell.bc$group, levels = c("2000", "1000","500", "250", "125","64"))


# Working on putting species (FAs) onto plot
species.dataframe <- data.frame(FA = rownames(eco.nmds.bc$species), 
                                NMDS.1 = eco.nmds.bc$species[,1],
                                NMDS.2 = eco.nmds.bc$species[,2])


# NMDS plot colored by size class
# With shape for season


# adding omega symbol
p <- ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+
  geom_path(data=df_ell.bc, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=2) +
  geom_point(stat = "identity", aes(fill=taxa.names,  shape=Seasons), size=2)+
  scale_fill_manual(values=c("#ffbf00",  "#ff7433","#db5764", "#9966ff", "#66ccff", "#47ebb4"), name="Size", guide="legend") +
  guides(fill= guide_legend(override.aes = list(shape=21)))+
  scale_color_manual(values=c("#ffbf00",  "#ff7433","#db5764","#9966ff","#66ccff", "#47ebb4"), guide=FALSE) + 
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(angle=90,size=9),
        axis.text.y=element_text(size=8),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"),
        panel.border = element_rect(linewidth = 0.3)) + 
  geom_segment(data = species.scores, aes(x = 0, xend=NMDS1*0.9, y=0, yend=NMDS2*0.9), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.4) + #add vector arrows of significant species
  # geom_text(data = species.dataframe, aes(NMDS.1, NMDS.2, label=FA)) + 
  #  16:0
  annotate("text", x=species.dataframe[c(2),2]*1.2, y=(species.dataframe[c(2),3]*1.5), 
           label=c("16:0"),  size=3.5) +
  # 16:1n-7
  annotate("text", x=species.dataframe[c(3),2]*1.5, y=species.dataframe[c(3),3]*1.2, 
           label=expression(paste("16:1",omega,"7")), fontface =2, size=3.5) + 
  # 14:0 
  annotate("text", x=species.dataframe[c(1),2], y=species.dataframe[c(1),3.5]*1.1, 
           label="14:0",  size=3.5) + 
  # 16:3n-4
  annotate("text", x=species.dataframe[c(5),2]*0.8, y=(species.dataframe[c(5),3]*1.05), 
           label=expression(paste("16:3",omega,"4")), fontface =2, size=3.5) +
  # 18:1n-7
  annotate("text", x=species.dataframe[c(6),2], y=(species.dataframe[c(6),3])*1.6, 
           label=expression(paste("18:1",omega,"7")), fontface =2, size=3.5) + 
  # 18:1n-9
  annotate("text", x=species.dataframe[c(7),2], y=(species.dataframe[c(7),3]*1.05), 
           label=expression(paste("18:1",omega,"9")), fontface =2, size=3.5)  + 
  # 18:4n-3
  annotate("text", x=species.dataframe[c(8),2]*0.4, y=(species.dataframe[c(8),3]*1.8), 
           label=expression(paste("18:4",omega,"3")), fontface =2, size=3.5) +
  # EPA
  annotate("text", x=species.dataframe[c(9),2]*1.65, y=(species.dataframe[c(9),3]*1.3), 
           label=expression(paste("20:5",omega,"3")), fontface =2, size=3.5) + 
  # DHA
  annotate("text", x=species.dataframe[c(10),2]*2.3, y=(species.dataframe[c(10),3]*1.2), 
           label=expression(paste("22:6",omega,"3")), fontface =2, size=3.5) + 
  # 16:2n-4
  annotate("text", x=species.dataframe[c(11),2]*0.88, y=(species.dataframe[c(11),3]*1.15), 
           label=expression(paste("16:2",omega,"4")), fontface =2, size=3.5) +
  # BFA and 18:0
  annotate("text", x=species.dataframe[c(12,4),2]*1.15, y=(species.dataframe[c(12,4),3]), 
           label=c("BAFA", "18:0"), size=3.5) +
  labs(x=element_blank(),
       y="NMDS2") + 
  scale_shape_manual(values = c(22,24, 21), name = "Season") + 
  geom_text(label="(b)", x=-0.285, y=0.22, color="black") + 
  coord_cartesian(clip = "off")

p 




# Colored by months

p2 <- ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+
  geom_path(data=df_ell.bc, aes(x=NMDS1, y=NMDS2, colour=group), size=1, linetype=2) +
  geom_point(stat = "identity", aes(fill=month, shape = season), size=2)+
  #to change color: fill = other_ids
  scale_fill_manual(values = c("#66C2A5", "#7FBC41", 
                               "#A6D96A", "#FEE08B", "#FDAE61", "#F46D43", 
                               "#D53E4F", "#9E0142", "#67001F", "white"), name="Month", guide="legend") +
  guides(fill= guide_legend(override.aes = list(shape=21)))+
  scale_color_manual(values=c("#ffbf00",  "#ff7433","#db5764","#9966ff","#66ccff", "#47ebb4"), guide=FALSE) + 
  theme_bw()+
  theme(axis.text.x=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(angle=90,size=9),
        axis.text.y=element_text(size=8),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"),
        panel.border = element_rect(linewidth = 0.3)) + 
  geom_segment(data = species.scores, aes(x = 0, xend=NMDS1*0.9, y=0, yend=NMDS2*0.9), arrow = arrow(length = unit(0.2, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  #  16:0
  annotate("text", x=species.dataframe[c(2),2]*1.2, y=(species.dataframe[c(2),3]*1.5), 
           label=c("16:0"),  size=3.5) +
  # 16:1n-7
  annotate("text", x=species.dataframe[c(3),2]*1.5, y=species.dataframe[c(3),3]*1.2, 
           label=expression(paste("16:1",omega,"7")), fontface =2, size=3.5) + 
  # 14:0 
  annotate("text", x=species.dataframe[c(1),2], y=species.dataframe[c(1),3.5]*1.1, 
           label="14:0",  size=3.5) + 
  # 16:3n-4
  annotate("text", x=species.dataframe[c(5),2]*0.8, y=(species.dataframe[c(5),3]*1.05), 
           label=expression(paste("16:3",omega,"4")), fontface =2, size=3.5) +
  # 18:1n-7
  annotate("text", x=species.dataframe[c(6),2], y=(species.dataframe[c(6),3])*1.6, 
           label=expression(paste("18:1",omega,"7")), fontface =2, size=3.5) + 
  # 18:1n-9
  annotate("text", x=species.dataframe[c(7),2], y=(species.dataframe[c(7),3]*1.05), 
           label=expression(paste("18:1",omega,"9")), fontface =2, size=3.5)  + 
  # 18:4n-3
  annotate("text", x=species.dataframe[c(8),2]*0.4, y=(species.dataframe[c(8),3]*1.8), 
           label=expression(paste("18:4",omega,"3")), fontface =2, size=3.5) +
  # # 20:1n-9
  # annotate("text", x=species.dataframe[c(9),2]*1.15, y=(species.dataframe[c(9),3]*1.15), 
  #          label=expression(paste("20:1",omega,"9")), fontface =2, size=3.5) + 
  # EPA
  annotate("text", x=species.dataframe[c(9),2]*1.65, y=(species.dataframe[c(9),3]*1.3), 
           label=expression(paste("20:5",omega,"3")), fontface =2, size=3.5) + 
  # DHA
  annotate("text", x=species.dataframe[c(10),2]*2.3, y=(species.dataframe[c(10),3]*1.2), 
           label=expression(paste("22:6",omega,"3")), fontface =2, size=3.5) + 
  # # 24:1n-9
  # annotate("text", x=species.dataframe[c(11),2]*1.2, y=(species.dataframe[c(11),3]*1.2), 
  #          label=expression(paste("24:1",omega,"9")), fontface =2, size=3.5) +
  # 16:2n-4
  annotate("text", x=species.dataframe[c(11),2]*0.88, y=(species.dataframe[c(11),3]*1.15), 
           label=expression(paste("16:2",omega,"4")), fontface =2, size=3.5) +
  # # 22:1n-11
  # annotate("text", x=species.dataframe[c(14),2]*0.86, y=(species.dataframe[c(14),3]*1.05), 
  #          label=expression(paste("22:1",omega,"11")), fontface =2, size=3.5) +
  # BFA and 18:0
  annotate("text", x=species.dataframe[c(12,4),2]*1.15, y=(species.dataframe[c(12,4),3]), 
           label=c("BAFA", "18:0"), size=3.5) +
  labs(x="NMDS1",
       y="NMDS2") + 
  scale_shape_manual(values = c(22,24, 21), name = "Season") + 
  geom_text(label="(d)", x=-0.285, y=0.22, color="black") + 
  coord_cartesian(clip = "off")

p2 








# SI biplot --------------------------------------------------------------

allData$Size.Fraction <- as.factor(allData$Size.Fraction)

allSIdata <- allData %>% select(Date, Size.Fraction, delta13c, delta15n, C_N, Month, season)
allSIdata.sm <- allSIdata %>% filter(Size.Fraction != "POM")
allSIdata.sm <- allSIdata.sm[complete.cases(allSIdata.sm),]
taxa.names.SI <- allSIdata.sm$Size.Fraction
months.SI <- allSIdata.sm$Month
seasons.SI  <- allSIdata.sm$season

allSIdata.sm <- allSIdata.sm %>% select(delta13c, delta15n)
allSIdata.sm$delta13c <- allSIdata.sm$delta13c * -1 

permanova_SI<-adonis2(allSIdata.sm ~ taxa.names.SI, permutations = 999, method="bray")
permanova_SI<-adonis2(allSIdata.sm ~ taxa.names.SI*seasons.SI, permutations = 999, method="bray")
permanova_SI #if significant, then plot it


#dataframe for plotting NMDS
NMDS.bc.SI<-data.frame(NMDS1.bc=allSIdata.sm[,1],NMDS2.bc=allSIdata.sm[,2],group=taxa.names.SI, month=months.SI, season=seasons.SI)

ord.bc.SI <-ordiellipse(allSIdata.sm, taxa.names.SI, display="sites",kind="sd", conf = 0.95, label=T)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.bc.SI <- data.frame()

for(g in levels(NMDS.bc.SI$group)){
  df_ell.bc.SI <- rbind(df_ell.bc.SI, cbind(as.data.frame(with(NMDS.bc.SI[NMDS.bc.SI$group==g,],
                                                               veganCovEllipse(ord.bc.SI[[g]]$cov,ord.bc.SI[[g]]$center))),group=g))
}

#https://www.rpubs.com/RGrieger/545184


# Pairwise comparison
pairwise.adonis(allSIdata.sm, taxa.names.SI)




pSIsizes <- ggplot(NMDS.bc.SI, aes(-NMDS1.bc, NMDS2.bc))+
  geom_path(data=df_ell.bc.SI, aes(x=-delta13c, y=delta15n, colour=group), size=1, linetype=2) +
  geom_point(stat = "identity", aes(fill=group,  shape=season), size=2)+
  #to change color: fill = other_ids
  scale_fill_manual(values=c("#ffbf00",  "#ff7433","#db5764", "#9966ff", "#66ccff", "#47ebb4", "#00b300"), name="Size", guide="legend") +
  guides(fill= guide_legend(override.aes = list(shape=21)))+
  scale_color_manual(values=c("#ffbf00",  "#ff7433","#db5764","#9966ff","#66ccff", "#47ebb4", "#00b300"), guide=FALSE) + 
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(angle=90,size=10),
        axis.text.y=element_text(size=8),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"),
        panel.border = element_rect(linewidth = 0.3)) + 
  labs(x=element_blank(),
       y=expression(paste(delta^{15}, "N"))) +  
  scale_shape_manual(values = c(22,24,21)) + 
  geom_text(label="(a)", x=-25.6, y=11.75, color="black") + 
  coord_cartesian(clip = "off")

pSIsizes
#   coord_fixed() + 


df_ell.bc.SI$group <- factor(df_ell.bc.SI$group, levels = c("2000", "1000","500", "250", "125","64", "POM"))

# Colored by months
pSImonths <- ggplot(NMDS.bc.SI, aes(-NMDS1.bc, NMDS2.bc))+
  geom_path(data=df_ell.bc.SI, aes(x=-delta13c, y=delta15n, colour=group), size=1, linetype=2) +
  geom_point(stat = "identity", aes(fill= month, shape=season), size=2)+
  #to change color: fill = other_ids
  scale_fill_manual(values = c( "#66C2A5", "#7FBC41", 
                                "#A6D96A", "#FEE08B", "#FDAE61", "#F46D43", 
                                "#D53E4F", "#9E0142", "#67001F", "#40004B"), name="Month", guide="legend") +
  guides(fill= guide_legend(override.aes = list(shape=c(21))))+
  scale_color_manual(values=c("#ffbf00",  "#ff7433","#db5764","#9966ff","#66ccff", "#47ebb4", "#00b300"), guide=FALSE) + 
  theme_bw()+
  theme(axis.text.x=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(angle=90,size=10),
        axis.text.y=element_text(size=8),
        panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
        panel.border = element_rect(linewidth = 0.3)) + 
  labs(x=expression(paste(delta^{13}, "C")),
       y=expression(paste(delta^{15}, "N"))) +  
  scale_shape_manual(values = c(22,24, 21)) + 
  geom_text(label="(c)", x=-25.6, y=11.75, color="black") + 
  coord_cartesian(clip = "off")

pSImonths



# plot for legend 
pnew2 <- ggplot(QU39.2015.64, aes(x=Date, y=C14.0_PERCENT)) + 
  geom_point(aes(color= Size.Fraction)) +
  theme_bw() +  
  scale_color_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                 "#9966ff",  "#66ccff", "#47ebb4")) + 
  labs(color = "Size") 
pnew2

legend.Zoop <- get_legend(
  # create some space to the left of the legend
  pnew2 + theme(legend.title = element_text(size=9),
                legend.text = element_text(size = 7),
                legend.position = c(0.5,0.6))
)

# plot for legend 
pnew3 <- ggplot(QU39.2015.64, aes(x=Date, y=C14.0_PERCENT)) + 
  geom_point(aes(shape=season, color= Month)) +
  theme_bw() +  
  scale_color_manual(values =  c("#66C2A5", "#7FBC41", 
                                 "#A6D96A", "#FEE08B", "#FDAE61", "#F46D43", 
                                 "#D53E4F", "#9E0142", "#67001F", "#40004B")) + 
  labs(color = "Month",
       shape = "") + 
  scale_shape_manual(values = c(22,24, 21)) 

pnew3

legend.months <- get_legend(
  # create some space to the left of the legend
  pnew3 + theme(legend.title = element_text(size=9),
                legend.text = element_text(size = 7),
                legend.position = c(0.5,0.7))
)


# Fig. 3 NMDS + Biplot ----------------------------------------------------

# Putting together in a grid 

prow <- plot_grid(pSIsizes, p, legend.Zoop, 
                  pSImonths, p2, legend.months,
                  rel_widths = c(0.5, 0.52, 0.12 ),
                  rel_heights = c(0.9, 1),
                  nrow = 2)
prow



png("20230106 NMDS and biplot Fig 4.png", width=180, height=150, units="mm", res=300)
prow
dev.off()







# IndVal analysis ------------------------------------------------------------------

# If I want to cut the tree at different heights then I need to specify h instead of k 

require(labdsv)

# Always use the raw data - not transformed 
taxa.names <- QU39.2015.64$Size.Fraction
(taxa.names <- droplevels(taxa.names))

# Use 
Q20.matrix.abund.all <- Q20.matrix[,abundant.all]

# Indicator species
(iva <- indval(Q20.matrix.abund.all, taxa.names, numitr = 10000))


# Correct the p-values for multiple testing:
pval.adj <- p.adjust(iva$pval)

gr <- iva$maxcls[pval.adj <= 0.05]
iv <- iva$indcls[pval.adj <= 0.05]
pv <- iva$pval[pval.adj <= 0.05]
fr <- apply(Q20.matrix.abund > 0, 2, sum)[pval.adj <= 0.05]
fidg <- data.frame(
  group = gr,
  indval = iv,
  pvalue = pv,
  freq = fr
)
fidg <- fidg[order(fidg$group, -fidg$indval), ]
fidg

# 2000 : C18.1n.9c, C18.1n.7, C16:0
# 1000 : 22:1n-11, 20:1n-9
# 500  : 16:1n-7, 20:5n-3 
# 250  : none
# 125  : 18:0, 14:0
# 64   : 15:0+17:0, 16:2n-4 
taxa.names







# time series plots (supplemental) ----------------------------------------



# time series of copepod FATM (Fig. S5)
ptimeseries <- ggplot(fatty.acid.all, aes(x=Date, y=I(C20MUFAs*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% C20 MUFAs")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")



ptimeseries
png("20230110 C20 MUFA time series.png", width=90, height=80, units="mm", res=300)
ptimeseries
dev.off()


# 4 panel diatom FATMs plot (Fig. S1)
fatty.acid.all$EPA.DHA <- fatty.acid.all$C20.5n.3_PERCENT / fatty.acid.all$C22.6n.3_PERCENT

pEPA.DHA <- ggplot(fatty.acid.all, aes(x=Date, y=EPA.DHA, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("EPA:DHA")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pEPA.DHA

pRatio16.1 <- ggplot(fatty.acid.all, aes(x=Date, y=Ratio16.1, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("16:1",omega,"7 / 16:0")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pRatio16.1


pC16PUFA <- ggplot(fatty.acid.all, aes(x=Date, y=I(C16PUFA*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% C16 PUFAs")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pC16PUFA

pD.F <- ggplot(fatty.acid.all, aes(x=Date, y=Diatom.Flag, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("D / F")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pD.F


prow <- plot_grid(pEPA.DHA, pRatio16.1, 
                  pC16PUFA, pD.F, 
                  align = 'v',
                  nrow = 2,
                  rel_heights = c(1,1))

png("20230110 diatom FATM time series.png", width=180, height=170, units="mm", res=300)
prow
dev.off()


# 4 panel flagellate FATMs plot (Fig. S3)
pDHA <- ggplot(fatty.acid.all, aes(x=Date, y=I(C22.6n.3_PERCENT*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% DHA")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pDHA

pSDA <- ggplot(fatty.acid.all, aes(x=Date, y=I(C18.4n.3_PERCENT*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% SDA (18:4",omega,"3)")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pSDA


pALA <- ggplot(fatty.acid.all, aes(x=Date, y=I(C18.3n.3_PERCENT*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% ALA (18:3",omega,"3)")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pALA

pC18.1n.7 <- ggplot(fatty.acid.all, aes(x=Date, y=I(C18.1n.7_PERCENT*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% 18:1",omega,"7")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pC18.1n.7


prow <- plot_grid(pDHA, pSDA, 
                  pALA, pC18.1n.7, 
                  align = 'v',
                  nrow = 2,
                  rel_heights = c(1,1))

png("20230110 flagellate FATM time series.png", width=180, height=170, units="mm", res=300)
prow
dev.off()



# 4 panel carnivory FATMs plot (Fig. S4)
pDHA.EPA <- ggplot(fatty.acid.all, aes(x=Date, y=DHA.EPA, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("DHA:EPA")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pDHA.EPA


pcarn18.1 <- ggplot(fatty.acid.all, aes(x=Date, y=Carn18.1N9_n7, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("18:1",omega,"9 / 18:1",omega,"7")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pcarn18.1

fatty.acid.all$SFA.PUFA
pSFA.PUFA <- ggplot(fatty.acid.all, aes(x=Date, y=SFA.PUFA, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("Sat. FA / PUFA")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pSFA.PUFA

pBAFA <- ggplot(fatty.acid.all, aes(x=Date, y=I(Bacteria_15_17*100), color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("% BAFA (15:0 + 17:0)")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31")),
               expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")
pBAFA


prow <- plot_grid(pDHA.EPA, pcarn18.1, 
                  pSFA.PUFA, pBAFA, 
                  align = 'v',
                  nrow = 2,
                  rel_heights = c(1,1))

png("20230110 carnivory FATM time series.png", width=180, height=170, units="mm", res=300)
prow
dev.off()




# Time series of total FA concentrations (Fig. S6) 
pSumFA <- ggplot(fatty.acid.all, aes(x=Date, y=SumFA_mg.g, color=Size.Fraction, fill=Size.Fraction)) + 
  geom_point() + theme_bw() + geom_smooth(se=F)  + 
  scale_color_manual(values = c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300")) + 
  scale_fill_manual(values =  c("#ffbf00",  "#ff7433","#db5764",
                                "#9966ff",  "#66ccff", "#47ebb4", "#00b300"))  +
  labs(y=expression(paste("Total FA (",mu,"g mg"^-1*")")), 
       x="", color="Size class" , fill="Size class") + 
  scale_x_date(breaks = as.Date(c("2015-01-01", "2015-04-01", "2015-07-01",
                                  "2015-10-01")),
               labels = c("Jan", "Apr", "Jul", "Oct"),
               limits = as.Date(c("2015-01-01", "2015-12-31"))) + 
  theme(legend.position = "none")

pSumFA
png("20221006 SumFA time series.png", width=120, height=100, units="mm", res=300)
pSumFA
dev.off()






