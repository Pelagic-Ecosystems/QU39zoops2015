
# 2022-04-04
# Building data table for 2015 QU39 Zooplankton analyses

library(tidyverse)
library(ggfortify)



# FA data -----------------------------------------------------------------

# QU39 Old Runs are samples analyzed prior to Anna joining PEL in 2019
# This file is made within the PEL fatty acid OneDrive folder
QU39oldruns <- read.csv("raw_data/QU39 Old Runs compiled 12.13.21.csv", stringsAsFactors = FALSE, na="")
str(QU39oldruns)
QU39oldruns$Date <- as.Date(QU39oldruns$Date, "%Y-%m-%d")

# Subset samples from 2015
QU39.2015 <- QU39oldruns %>% filter(Date < as.Date("2015-12-31"))

# Make variable for net mesh size
QU39.2015 <- QU39.2015 %>%  mutate(mesh.size = case_when(grepl("_64_", Hakai.Net.ID) ~ "64 um",
                                                         grepl("_250_", Hakai.Net.ID) ~ "250 um"))

# Make two separate data sets for the two mesh sizes
QU39.2015.64 <- QU39.2015 %>% filter(mesh.size == "64 um")
QU39.2015.250 <- QU39.2015 %>% filter(mesh.size == "250 um")



# Calculate Markers -------------------------------------------------------
# These have already been done for the POM FA data package

# 16:1n7/16:0
QU39.2015.64$Ratio16.1 <- (QU39.2015.64$C16.1n.7_PERCENT/ QU39.2015.64$C16.0_PERCENT)

# Diatom FAs / Flagellate FAs
QU39.2015.64$Diatom.Flag <- ((QU39.2015.64$C16.1n.7_PERCENT + QU39.2015.64$C20.5n.3_PERCENT +  QU39.2015.64$Diatom.III_PERCENT + QU39.2015.64$C18.1n.12_PERCENT) 
                             / (QU39.2015.64$C22.6n.3_PERCENT + QU39.2015.64$C18.3n.3_PERCENT + QU39.2015.64$C18.3n.6_PERCENT + QU39.2015.64$C18.4n.3_PERCENT))

QU39.2015.64$C16PUFA <- rowSums(QU39.2015.64[,c("Diatom.III_PERCENT", "C18.1n.12_PERCENT")], na.rm=TRUE)

# DHA/EPA : Dinos/Diatoms
QU39.2015.64$DHA.EPA <- (QU39.2015.64$C22.6n.3_PERCENT/ QU39.2015.64$C20.5n.3_PERCENT)

# 18:2(n-6)+ 18:3(n-3) > 2.5 Terrestrial signal
QU39.2015.64$Terr_18.2n6_18.3n3 <- (QU39.2015.64$C18.2n.6c_PERCENT + QU39.2015.64$C18.3n.3_PERCENT)

# 15:0 and 17:0 : Bacteria
QU39.2015.64$Bacteria_15_17 <- (QU39.2015.64$C15.0_PERCENT + QU39.2015.64$C17.0_PERCENT)

# 18:1n-9/18:1n-7 Carnivory or omnivory
QU39.2015.64$Carn18.1N9_n7 <- (QU39.2015.64$C18.1n.9c_PERCENT/ QU39.2015.64$C18.1n.7_PERCENT)

# Make 18:3n-3_18:4n-4 (pico-Chl)
QU39.2015.64$C18.3andC18.4 <- QU39.2015.64$C18.3n.3_PERCENT + QU39.2015.64$C18.4n.3_PERCENT

QU39.2015.64$C20MUFAs <- QU39.2015.64$C20.1n.9_PERCENT + QU39.2015.64$C22.1n.11_PERCENT + QU39.2015.64$C20.1n.11_PERCENT + QU39.2015.64$C22.1n.9_PERCENT


# Calculate proportions for each main class

percent.SFA <- vector(length=nrow(QU39.2015.64))
for(i in 1:nrow(QU39.2015.64)){
  percent.SFA[i]=sum(
    QU39.2015.64$C10.0_PERCENT[i], 
    QU39.2015.64$C11.0_PERCENT[i], 
    QU39.2015.64$C12.0_PERCENT[i], 
    QU39.2015.64$C13.0_PERCENT[i], 
    QU39.2015.64$C14.0_PERCENT[i], 
    QU39.2015.64$C15.0_PERCENT[i], 
    QU39.2015.64$C16.0_PERCENT[i], 
    QU39.2015.64$C17.0_PERCENT[i], 
    QU39.2015.64$C18.0_PERCENT[i], 
    QU39.2015.64$C20.0_PERCENT[i],
    QU39.2015.64$C22.0_PERCENT[i],
    QU39.2015.64$C23.0_PERCENT[i],
    QU39.2015.64$C24.0_PERCENT[i],
    na.rm=TRUE)
}
QU39.2015.64$percent.SFA <- percent.SFA

percent.MUFA <- vector(length=nrow(QU39.2015.64))
for(i in 1:nrow(QU39.2015.64)){
  percent.MUFA[i]=sum(QU39.2015.64$C14.1_PERCENT[i], 
                      QU39.2015.64$C16.1n.7_PERCENT[i], 
                      QU39.2015.64$C18.1n.7_PERCENT[i], 
                      QU39.2015.64$C18.1n.9c_PERCENT[i], 
                      QU39.2015.64$C20.1n.9_PERCENT[i],
                      QU39.2015.64$C22.1n.9_PERCENT[i],
                      QU39.2015.64$C24.1n.9_PERCENT[i],na.rm=TRUE)
}
QU39.2015.64$percent.MUFA <- percent.MUFA

percent.PUFA <- vector(length=nrow(QU39.2015.64))
for(i in 1:nrow(QU39.2015.64)){
  percent.PUFA[i]=sum(QU39.2015.64$C18.2n.6c_PERCENT[i], 
                      QU39.2015.64$C18.3n.3_PERCENT[i], 
                      QU39.2015.64$C18.3n.6_PERCENT[i], 
                      QU39.2015.64$C18.4n.3_PERCENT[i],
                      QU39.2015.64$Diatom.III_PERCENT[i],
                      QU39.2015.64$C18.1n.12_PERCENT[i], 
                      QU39.2015.64$C20.3n.3_PERCENT[i], 
                      QU39.2015.64$C20.4n.6_PERCENT[i], 
                      QU39.2015.64$C20.5n.3_PERCENT[i],
                      QU39.2015.64$C22.2n.6_PERCENT[i], 
                      QU39.2015.64$C22.4n.6_PERCENT[i], 
                      QU39.2015.64$C22.5n.3_PERCENT[i], 
                      QU39.2015.64$C22.5n.6._PERCENT[i],
                      QU39.2015.64$C22.6n.3_PERCENT[i], na.rm=TRUE)
} 
QU39.2015.64$percent.PUFA <- percent.PUFA

percent.n3_PUFA <- vector(length=nrow(QU39.2015.64))
for(i in 1:nrow(QU39.2015.64)){
  percent.n3_PUFA[i]=sum(QU39.2015.64$C18.3n.3_PERCENT[i], 
                         QU39.2015.64$C18.4n.3_PERCENT[i],
                         QU39.2015.64$C20.3n.3_PERCENT[i], 
                         QU39.2015.64$C20.5n.3_PERCENT[i],
                         QU39.2015.64$C22.5n.3_PERCENT[i], 
                         QU39.2015.64$C22.6n.3_PERCENT[i], na.rm=TRUE)
} 
QU39.2015.64$percent.n3_PUFA <- percent.n3_PUFA

percent.n6_PUFA <- vector(length=nrow(QU39.2015.64))
for(i in 1:nrow(QU39.2015.64)){
  percent.n6_PUFA[i]=sum(QU39.2015.64$C18.2n.6_PERCENT[i], 
                         QU39.2015.64$C18.3n.6_PERCENT[i], 
                         QU39.2015.64$C20.4n.6_PERCENT[i], 
                         QU39.2015.64$C22.2n.6_PERCENT[i], 
                         QU39.2015.64$C22.4n.6_PERCENT[i], 
                         QU39.2015.64$C22.5n.6_PERCENT[i], na.rm=TRUE)
} 
QU39.2015.64$percent.n6_PUFA <- percent.n6_PUFA


# PUFA/SFA carnivory index
QU39.2015.64$PUFA.SFA <- QU39.2015.64$percent.PUFA/QU39.2015.64$percent.SFA

# n-3 : n-6 ratio (too many samples w zero n-6 in POMs)
QU39.2015.64 <- QU39.2015.64 %>% mutate(n3.n6 = percent.n3_PUFA/percent.n6_PUFA)

# EFA = 18C PUFA and longer; HUFA = 20C PUFA and longer
QU39.2015.64 <- QU39.2015.64 %>% mutate(EFA.mg.g = ((C18.3n.6_PERCENT + C18.3n.3_PERCENT + C18.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT + C22.5n.6_PERCENT + C22.4n.6_PERCENT + C22.5n.3_PERCENT + C20.3n.3_PERCENT)* SumFA_mg.g))
QU39.2015.64 <- QU39.2015.64 %>% mutate(HUFA.mg.g = ((C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT + C22.5n.6_PERCENT + C22.4n.6_PERCENT + C22.5n.3_PERCENT + C20.3n.3_PERCENT)* SumFA_mg.g))

QU39.2015.64 <- QU39.2015.64 %>% mutate(Prop.EFA = (C18.3n.6_PERCENT + C18.3n.3_PERCENT + C18.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT + C22.5n.6_PERCENT + C22.4n.6_PERCENT + C22.5n.3_PERCENT + C20.3n.3_PERCENT))
QU39.2015.64 <- QU39.2015.64 %>% mutate(Prop.HUFA = (C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT + C22.5n.6_PERCENT + C22.4n.6_PERCENT + C22.5n.3_PERCENT + C20.3n.3_PERCENT))



# POMFA data --------------------------------------------------------------

# Published data package from McLaskey et al. 2021 
# However there are more Chl and SI data than in this data package
POMFAdata <- read.csv("raw_data/POMFA data package Final 11.5.21.csv", stringsAsFactors = FALSE)
POMFAdata$Date <- as.Date(POMFAdata$Date, "%m/%d/%y")

# Filter for 2015 samples only
POMFA.2015 <- POMFAdata %>% filter(Date < as.Date("2015-12-31"))

POMFA.2015 <- POMFA.2015 %>% select(Date, Volume.filtered_mL, C14.0_PERCENT:Prop.HUFA) 

colnames(POMFA.2015)[16] <- "C18.2n.6c_PERCENT"
colnames(QU39.2015.64)[c(101, 71)] <- c("C16.2n.4_PERCENT", "C16.3n.4_PERCENT")

# Change the units/scale on POM TFA ug/L so it will plot w zooplankton ug/mg
POMFA.2015$SumFA_mg. <- (POMFA.2015$TotalFA_ug / POMFA.2015$Volume.filtered_mL)*2000
POMFA.2015$SumFA_ug.L <- (POMFA.2015$TotalFA_ug / POMFA.2015$Volume.filtered_mL)*1000
POMFA.2015$Size.Fraction..um. <- "POM"

POMFA.2015$Size.Fraction..um. <- as.factor(POMFA.2015$Size.Fraction..um.)
QU39.2015.64$Size.Fraction..um. <- as.factor(QU39.2015.64$Size.Fraction..um.)

# Missing FATM
POMFA.2015$C16PUFA <- rowSums(POMFA.2015[,c("C16.2n.4_PERCENT", "C16.3n.4_PERCENT")], na.rm=TRUE)



# Join POM FA  w ZoopFA ----------------------------------------------------------

fatty.acid.2015all  <- full_join(POMFA.2015, QU39.2015.64)

fatty.acid.2015all$Size.Fraction..um. <- factor(fatty.acid.2015all$Size.Fraction..um., levels = c("2000", "1000","500", "250", "125","64", "POM"))




# POM SI ------------------------------------------------------------------

# THere are more SI samples than are in the POMFA data package 

POMsi <- read.csv("raw_data/POM SI QU39 all dates surface 5m.csv", stringsAsFactors = FALSE)
 
POMsi$Date <- as.Date(POMsi$Date)

# Only the year of the time series
POMsi.2015 <- POMsi %>% filter(Date < as.Date("2015-12-31"))
colnames(POMsi.2015)[5:6] <- c("delta15n", "delta13c")
 
# Calculate correction for first two non-acidified carbon samples 
POMsi.2015$Del13c.diff <- POMsi.2015$delta13c.notacid -POMsi.2015$delta13c
mean.offset <- mean(POMsi.2015$Del13c.diff[8:17])
POMsi.2015$delta13c[1:6] <- POMsi.2015$delta13c.notacid[1:6] - mean.offset

ggplot(POMsi.2015, aes(x=Date, y=delta13c)) + geom_point() + theme_bw() + 
   geom_smooth()
 
POMsi.2015.sm <- POMsi.2015 %>% select(Date, delta13c, delta15n)
POMsi.2015.sm$Size.Fraction..um. <- "POM"




# Zoop SI -----------------------------------------------------------------

# Load the data file from the Portal
zoopSI <- read.csv("raw_data/2022-07-14_095901_HakaiData_zooplankton_isotope.csv", na="")
zoopSI$Date <- as.Date(zoopSI$Date, "%Y-%m-%d")
colnames(zoopSI)
# Make smaller dataset
zoopSI.sm <- zoopSI %>% select(Hakai.Net.ID, Date, Site.ID, Size.Fraction..um., Hakai.ID, 
                               Analyzing.Lab,  corr_delta15n, corr_delta13c, ug.C, ug.N, C.Flag, N.Flag, Sample.Status) 

zoopSI.sm <- zoopSI.sm %>% filter(Sample.Status != "Not Available")

# Add column for net mesh size 
zoopSI.sm["Net.Mesh.Size"] <- NA 
for(i in 1:nrow(zoopSI.sm)){
  if(grepl("64", zoopSI.sm$Hakai.Net.ID[i])){
    zoopSI.sm$Net.Mesh.Size[i] <- "64 um"
  }   else if(grepl("250", zoopSI.sm$Hakai.Net.ID[i]))
  { zoopSI.sm$Net.Mesh.Size[i] <- "250 um"}}

# This is the exact same thing I have done with the biomass data
# Compare the number of unique net IDs to the number of mesh sizes
netsperday <- zoopSI.sm %>% select(Date, Hakai.Net.ID)
meshesperday <- zoopSI.sm %>% select(Date, Net.Mesh.Size)

netsperday.nodups <- netsperday[duplicated(netsperday)==FALSE,]
netsperday.nodups <- netsperday.nodups[complete.cases(netsperday.nodups)==TRUE,]
meshesperday.nodups <- meshesperday[duplicated(meshesperday)==FALSE,]
meshesperday.nodups <- meshesperday.nodups[complete.cases(meshesperday.nodups)==TRUE,]

numNetsDay <- data.frame(table(netsperday.nodups$Date))
colnames(numNetsDay)[2] <- "Num.Nets"
numMeshSizeDay <- data.frame(table(meshesperday.nodups$Date))
colnames(numMeshSizeDay)[2] <- "Num.Meshes"

netSummary <- full_join(numNetsDay, numMeshSizeDay)
netSummary$difference <- netSummary$Num.Nets - netSummary$Num.Meshes
table(netSummary$difference)
# five dates have more nets cast than mesh sizes, i.e. more than one cast of the same mesh
# 6/6/2017      64 um net is called 250
# 6/13/2017     64 um net is called 250
# 7/25/2017     64 um net is called 250
# 12/4/2017     the FA samples from this day are listed as isotope
# 5/29/2019     64 um net is called 250

# Manually change the mesh size on these nets
zoopSI.sm$Net.Mesh.Size[zoopSI.sm$Hakai.Net.ID=="2017-06-06_QU39_4_250_3"] <- "64 um"
zoopSI.sm$Net.Mesh.Size[zoopSI.sm$Hakai.Net.ID=="2017-06-13_QU39_4_250_1"] <- "64 um"
zoopSI.sm$Net.Mesh.Size[zoopSI.sm$Hakai.Net.ID=="2017-07-25_QU39_4_250_1"] <- "64 um"
zoopSI.sm$Net.Mesh.Size[zoopSI.sm$Hakai.Net.ID=="2019-05-29_QU39_6_250_3"] <- "64 um"
# Remove the FA samples that are mislabeled as isotopes
zoopSI.sm <- zoopSI.sm[zoopSI.sm$Hakai.Net.ID!="2017-12-04_QU39_3_250_2",]  

# Check quality flags
table(zoopSI.sm$C.Flag)
table(zoopSI.sm$N.Flag)

# Create unique ID for each size fraction within each net
zoopSI.sm["Net_sizefrac"] <- NA 
for(i in 1:nrow(zoopSI.sm)){
  zoopSI.sm$Net_sizefrac[i] <- paste(zoopSI.sm$Hakai.Net.ID[i], zoopSI.sm$Size.Fraction[i], sep="_")
}

# Calculate C:N
zoopSI.sm <- zoopSI.sm %>% mutate(C_N =  ((ug.C*12.01) / (ug.N*14.001)) )

# Do a lipid correction on Del13C following regression for SoG copepods
# δ13C_corrected = δ13C_bulk + (0.38 * C:N_bulk) - 1.85  (El-Sabaawi et al. 2009)
zoopSI.sm <- zoopSI.sm %>% mutate(delta13c_corrected =  corr_delta13c + (0.38 * C_N) - 1.85 )
colnames(zoopSI.sm)
colnames(zoopSI.sm)[c(4,7:8,17)] <- c("Size.Fraction", "delta15n", "delta13c_uncorrected", "delta13c")
zoopSI.sm$Size.Fraction <- as.factor(zoopSI.sm$Size.Fraction)

# I haven't filtered any flags yet so C or N doesn't matter 
zoopSI.sm.complete <- zoopSI.sm[complete.cases(zoopSI.sm$delta13c),]

ggplot(zoopSI.sm.complete, aes(x=Date, y=delta13c, color=Size.Fraction)) + geom_point() + theme_bw()

# Check for filters that were run separate, but actually same sample 
zoopSI.sm.complete[duplicated(zoopSI.sm.complete$Net_sizefrac),]
zoopSI.sm.complete[duplicated(zoopSI.sm.complete$Net_sizefrac, fromLast = T),]
# 2016-07-31_QU39_2_250_2_500: QF2238 and QF2236
# 2018-05-25_QU39_1_64_3_4000: QF4526 and QF4629
# Neither will be in this dataset 

# Filter out data flags
zoopSI.no.flags <- zoopSI.sm.complete
zoopSI.no.flags$delta13c[!is.na(zoopSI.no.flags$C.Flag)] <- NA
zoopSI.no.flags$delta15n[!is.na(zoopSI.no.flags$N.Flag)] <- NA



# Filter to 2015 only
# There are only a very few 64 and 125 fractions prior to March, drop them 
# There are a very few samples from the 250 net on 5/11/15 - which I am keeping in for completeness
# zoopSI.2015 <- zoopSI.no.flags %>% filter(Date < as.Date("2015-12-31") & 
#                                    Date > as.Date("2015-03-01"))

zoopSI.2015 <- zoopSI.no.flags %>% select(Date, Site.ID, Size.Fraction, delta15n, delta13c, 
                                      ug.C, ug.N, C_N)

colnames(POMsi.2015.sm)[4] <- "Size.Fraction"


# Join ZoopSI and POMSI
all.2015.SI <- full_join(zoopSI.2015, POMsi.2015.sm)


# Join AllSI and AllFA
fatty.acid.2015all <- fatty.acid.2015all %>% rename(Size.Fraction = Size.Fraction..um.)

full.2015.data <- full_join(fatty.acid.2015all, all.2015.SI)





# Chlorophyll -------------------------------------------------------------

chl.phaeo_all <- read.csv("raw_data/QU39 Chlorophyll data.csv", stringsAsFactors = FALSE)
str(chl.phaeo_all)
chl.phaeo_all$Date <- as.Date(chl.phaeo_all$Date, "%Y-%m-%d")


Chl.2015 <- chl.phaeo_all %>% filter(Date < as.Date("2015-12-31") & Date > as.Date("2015-01-01")) %>% 
  select(Date, chl_20um, chl_3um, chl_GF.F, chl_Bulk.GF.F, SumChl)

str(Chl.2015)

# Chlorophyll proportions
Chl.2015 <- Chl.2015 %>%   mutate(chla_prop_GF.F = chl_GF.F / (chl_20um + chl_GF.F + chl_3um))  %>% 
  mutate(chla_prop_20um = chl_20um / (chl_20um + chl_GF.F + chl_3um))  %>% 
  mutate(chla_prop_3um = chl_3um / (chl_20um + chl_GF.F + chl_3um))  


Chl.2015$Size.Fraction <- "POM"

full.2015.wChl <- full_join(full.2015.data, Chl.2015)



# Biomass -----------------------------------------------------------------

data64nets <- read.csv("raw_data/Biomass_64net_20221215.csv")
str(data64nets)
data64nets$Date <- as.Date(data64nets$Date, "%Y-%m-%d")

data64nets.2015 <- data64nets %>% filter(Date < as.Date("2015-12-31") & Date > as.Date("2015-01-01"))
data64nets.2015 <- data64nets.2015[!is.na(data64nets.2015$screen_size),]


data64nets.2015.sm <- data64nets.2015 %>% select(Date, screen_size, Biomass.mg.m3)
colnames(data64nets.2015.sm)[2] <- "Size.Fraction"
data64nets.2015.sm$Size.Fraction <- as.factor(data64nets.2015.sm$Size.Fraction)

data64nets.2015.sm <- data64nets.2015.sm[data64nets.2015.sm$Size.Fraction!= "4000",]

full.2015.wChl.biomass <- full_join(full.2015.wChl, data64nets.2015.sm)



# write file --------------------------------------------------------------

# Sort data to match original data file (only influences ordination direction)
full.2015.wChl.biomass.new <- full.2015.wChl.biomass %>% arrange(Size.Fraction)
full.2015.wChl.biomass.new <- full.2015.wChl.biomass.new %>% arrange(Date)


write.csv(full.2015.wChl.biomass.new, "processed_data/QU39 2015 zoop POM FA SI Chl biomass 20231102.csv", row.names = F)

full.2015.wChl.biomass.new <- read.csv("processed_data/QU39 2015 zoop POM FA SI Chl biomass 20231102.csv")

colnames(full.2015.wChl.biomass.new)

full.2015.wChl.biomass.new$EPA.DHA <- full.2015.wChl.biomass.new$C20.5n.3_PERCENT / full.2015.wChl.biomass.new$C22.6n.3_PERCENT

full.2015.wChl.biomass.new2 <- full.2015.wChl.biomass.new %>% select(Date, Size.Fraction, 
                                                                     SumFA_ug.L, SumFA_mg.g, SumFA_mg.g_WetWt, 
                                                                     C14.0_PERCENT:C18.1n.7_PERCENT,
                                                                     C18.1n.9c_PERCENT:C18.2n.6c_PERCENT, 
                                                                     C18.3n.3_PERCENT:C22.1n.11_PERCENT,
                                                                    Ratio16.1:DHA.EPA, Bacteria_15_17:n3.n6, 
                                                                    Prop.EFA, Prop.HUFA, C16PUFA, EPA.DHA, C20MUFAs, 
                                                                    delta15n, delta13c, 
                                                                    chl_20um:chl_GF.F, SumChl:Biomass.mg.m3)


write.csv(full.2015.wChl.biomass.new2, "processed_data/QU39 2015 zoop POM FA SI Chl biomass 20231102.csv", row.names = F)
