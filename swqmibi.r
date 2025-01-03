####################################################################
# 
# MACROINVERTEBRATE METRIC and TEXAS IBI CALCULATION SCRIPTS
#
# 
#
# Ely Kosnicki 2017
####################################################################
setwd("~")
library(Hmisc); 
library(doBy); 
library(plyr);
bugs<-read.csv("RBP2021.csv", header=T)
bugs$Date = as.Date(bugs$Date, "%d-%b-%y") #Reformat date
bugs<-ddply(bugs, .(Site,	Location,	Date,	Season,	Class, Order,	Family,
                      FinalID,	num,	TolTX,	FFGTX,	FFGTX2),
              summarize, num = sum(num)) #make sure totals per taxon 
total<-ddply(bugs, .(Site, Location, Season, Date), summarize, total=sum(num))
rich<-ddply(bugs, .(Site, Location, Season, Date), summarize, Richness=length(num))

EPTprep<-bugs
EPTprep["EPT"]<-ifelse(bugs$Order=="Trichoptera"|bugs$Order=="Ephemeroptera"|
                        bugs$Order=="Plecoptera",1,0)
EPT<-ddply(EPTprep, .(Site, Location, Season, Date), summarize, EPT=sum(EPT))

richprep<-bugs
richprep["Noninsects"]<-ifelse(richprep$Class != "Insecta",1,0)
Noinrich<-ddply(richprep, .(Site, Location, Season, Date), summarize, Noninsects=sum(Noninsects))

##Calculate %Dominant
maxprep<-bugs
maxes <- function(x){max=max(x)}
dom1<-summaryBy(num~Site+Location+Season+Date, data=maxprep, FUN=maxes)
dom1["Dominant"]<-100*dom1$num.maxes/total$total
dom<-dom1[,-5]

##CALCULATE BIOTIC INDEX
biprep<-bugs
bired.1<-biprep[biprep$TolTX>-1,] #remove artifact of BI values <0
bired1<-bired.1[complete.cases(bired.1[,9]),] #remove na.s
bitotal<-ddply(bired1, .(Site, Location, Season, Date), summarize, total=sum(num))
bired2<-merge(bired1,bitotal)
bired2["BI"]<-bired2$TolTX*bired2$num/bired2$total
bi<-aggregate(BI~Site+Location+Season+Date, bired2, sum)

bired1["intol"]<-ifelse(bired1$TolTX<6,bired1$num,0)
intoltax<-ddply(bired1, .(Site, Location, Season, Date), summarize, intol=sum(intol))

bired1["tolerant"]<-ifelse(bired1$TolTX>=6,bired1$num,0)
toltax<-ddply(bired1, .(Site, Location, Season, Date), summarize, tolerant=sum(tolerant))

#Must put in a "Divider" which is really a multiplier at the end of the FFG
#and other composition calculations
ffgprep<-bugs
ffgprep["Divider"]<-ifelse(ffgprep$FFGTX=="" & ffgprep$FFGTX2=="",0,
                           ifelse(ffgprep$FFGTX2=="",1, 0.5))
ffgprep["chiro"]<-ifelse(ffgprep$Family=="Chironomidae",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
chiro<-ddply(ffgprep, .(Site, Location, Season, Date), summarize, chiro=sum(chiro))

ffgprep["pred"]<-ifelse(ffgprep$FFGTX=="Predator"|ffgprep$FFGTX2=="Predator",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
pred<-ddply(ffgprep, .(Site, Location, Season, Date, Divider), summarize, pred=sum(pred))

ffgprep["gat"]<-ifelse(ffgprep$FFGTX=="Gather/Collector"|ffgprep$FFGTX2=="Gather/Collector",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
gat<-ddply(ffgprep, .(Site, Location, Season, Date, Divider), summarize, gat=sum(gat))

ffgprep["shred"]<-ifelse(ffgprep$FFGTX=="Shredder"|ffgprep$FFGTX2=="Shredder",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
shred<-ddply(ffgprep, .(Site, Location, Season, Date, Divider), summarize, shred=sum(shred))

ffgprep["scrape"]<-ifelse(ffgprep$FFGTX=="Scraper"|ffgprep$FFGTX2=="Scraper",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
scrape<-ddply(ffgprep, .(Site, Location, Season, Date, Divider), summarize, scrape=sum(scrape))

ffgprep["filt"]<-ifelse(ffgprep$FFGTX=="Filterer/Collector"|ffgprep$FFGTX2=="Filterer/Collector",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
filt<-ddply(ffgprep, .(Site, Location, Season, Date, Divider), summarize, filt=sum(filt))

ffgprep["elm"]<-ifelse(ffgprep$Family=="Elmidae",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
elm<-ddply(ffgprep, .(Site, Location, Season, Date), summarize, elm=sum(elm))

ffgprep["trichop"]<-ifelse(ffgprep$Order=="Trichoptera",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
trichop<-ddply(ffgprep, .(Site, Location, Season, Date), summarize, trichop=sum(trichop))

ffgprep["hydrop"]<-ifelse(ffgprep$Family=="Hydropsychidae",ffgprep$num,0)
ffgprep[is.na(ffgprep)]<-0
hydrop<-ddply(ffgprep, .(Site, Location, Season, Date), summarize, hydrop=sum(hydrop))

##COMBINE TOTAL WITH FFG, HABIT, AND TAXA TOTALS TO CALCULATE %FFGs
ffgtots<-Reduce(function(x,y) merge(x,y, all=T), list(total, intoltax, toltax, 
    chiro, hydrop, trichop, elm))
ffgtots["Chironomidae"]<-100*chiro$chiro/total$total
ffgtots["Elmidae"]<-100*elm$elm/total$total
ffgtots["Hydropsychidae_Trichoptera"]<-100*hydrop$hydrop/trichop$trichop
ffgtots["Intolerant_Tolerant"]<-intoltax$intol/toltax$tolerant
ffgtots[is.na(ffgtots)]<-0

ffgs<-Reduce(function(x,y) merge(x,y, all=T), list(gat, pred, filt, shred, scrape))
ffgs["pgat"]<-100*gat$gat*gat$Divider
GGat<-ddply(ffgs, .(Site, Location, Season, Date), summarize, pgat=sum(pgat))
ffgs["ppred"]<-100*pred$pred*pred$Divider
PPred<-ddply(ffgs, .(Site, Location, Season, Date), summarize, ppred=sum(ppred))
ffgs["pfilt"]<-100*filt$filt*ffgs$Divider
FFilt<-ddply(ffgs, .(Site, Location, Season, Date), summarize, pfilt=sum(pfilt))
ffgs["pshred"]<-100*shred$shred*ffgs$Divider
SShred<-ddply(ffgs, .(Site, Location, Season, Date), summarize, pshred=sum(pshred))
ffgs["pscrape"]<-100*scrape$scrape*ffgs$Divider
SScrape<-ddply(ffgs, .(Site, Location, Season, Date), summarize, pscrape=sum(pscrape))

ffgtots["Gatherers"]<-GGat$pgat/total$total
ffgtots["Predators"]<-PPred$ppred/total$total
ffgtots["Filterers"]<-FFilt$pfilt/total$total
ffgtots["Shredders"]<-SShred$pshred/total$total
ffgtots["Scrapers"]<-SScrape$pscrape/total$total

##Calculate %Dominant FFG
ffgdom2<-ffgtots[,c("Gatherers","Filterers", "Predators", "Shredders", "Scrapers")]
ffgdom2$Guild<-apply(ffgdom2, 1, max)
ffgdom2$label<-row.names(ffgdom2)
total1<-total
total1$label<-row.names(total1)
ffgdom1<-merge(total1, ffgdom2, by="label")
ffgdom<-ffgdom1[,c(2:5,12)]

##Calculate SWQM MMI
##Metric list: rich, ept, bi, pchiro, pdom1, ffgdom, ppred, intol_tol, phydrop_tri,
##Noinrich, pgat, pelmidae
metrics<-Reduce(function(x,y) merge(x,y, all=T), list(rich, EPT, bi, dom,
                                                      ffgtots, ffgdom, Noinrich))

mmiprep<-metrics[,c("Site", "Location", "Season", "Date", "Richness", "EPT", "BI", "trichop",
            "Chironomidae", "Dominant", "Guild", "Predators", "Intolerant_Tolerant",
            "Hydropsychidae_Trichoptera", "Noninsects", "Gatherers", "Elmidae")]

mmiprep["sc_rich"]<-ifelse(mmiprep$Richness>21,4, ifelse(mmiprep$Richness>=15, 3,
                    ifelse(mmiprep$Richness>=8, 2, 1)))

mmiprep["sc_ept"]<-ifelse(mmiprep$EPT>9,4, ifelse(mmiprep$EPT>=7, 3,
                    ifelse(mmiprep$EPT>=4, 2, 1)))

mmiprep["sc_bi"]<-ifelse(mmiprep$BI>5.27,1, ifelse(mmiprep$BI>=4.53, 2,
                    ifelse(mmiprep$BI>=3.77, 3, 4)))

mmiprep["sc_pchiro"]<-ifelse(mmiprep$Chironomidae<0.79|mmiprep$Chironomidae>16.19,1, 
        ifelse(mmiprep$Chironomidae<4.11, 4, ifelse(mmiprep$Chironomidae<9.49, 3, 2)))

mmiprep["sc_pdom1"]<-ifelse(mmiprep$Dominant>39.88,1, ifelse(mmiprep$Dominant>=31.02, 2,
                    ifelse(mmiprep$Dominant>=22.15, 3, 4)))

mmiprep["sc_ffgdom"]<-ifelse(mmiprep$Guild>54.12,1, ifelse(mmiprep$Guild>=45.31, 2,
                    ifelse(mmiprep$Guild>=36.5, 3, 4)))

mmiprep["sc_ppred"]<-ifelse(mmiprep$Predators<4.73|mmiprep$Predators>36.14,1, 
        ifelse(mmiprep$Predators<15.21, 4, ifelse(mmiprep$Predators<25.68, 3, 2)))

mmiprep["sc_intol_tol"]<-ifelse(mmiprep$Intolerant_Tolerant>4.79,4, ifelse(mmiprep$Intolerant_Tolerant>=3.21, 3,
                    ifelse(mmiprep$Intolerant_Tolerant>=1.63, 2, 1)))

mmiprep["sc_phydrop_tri"]<-ifelse(mmiprep$trichop==0|mmiprep$Hydropsychidae_Trichoptera>75.5,1, 
        ifelse(mmiprep$Hydropsychidae_Trichoptera<25.5, 4, ifelse(mmiprep$Hydropsychidae_Trichoptera<50.51, 3, 2)))

mmiprep["sc_Noinrich"]<-ifelse(mmiprep$Noninsects>5,4, ifelse(mmiprep$Noninsects>=4, 3,
                    ifelse(mmiprep$Noninsects>=2, 2, 1)))

mmiprep["sc_pgat"]<-ifelse(mmiprep$Gatherers<8|mmiprep$Gatherers>41.68,1, 
        ifelse(mmiprep$Gatherers<19.24, 4, ifelse(mmiprep$Gatherers<30.47, 3, 2)))

mmiprep["sc_pelmidae"]<-ifelse(mmiprep$Elmidae<0.88|mmiprep$Elmidae>30.12,1, 
        ifelse(mmiprep$Elmidae<10.05, 4, ifelse(mmiprep$Elmidae<20.09, 3, 2)))

mmiprep["MMI"]<-apply(mmiprep[,18:29], 1, sum)    ##Calculate the score and
                                             ##Define the Aquatic-life-use point-score
str(mmiprep)

mmiprep["Aquatic-life-use"]<-ifelse(mmiprep$MMI>36,"Exceptional",
        ifelse(mmiprep$MMI>28, "High", ifelse(mmiprep$MMI>21, "Intermediate", "Limited")))

##Take out mmi metric scores
##Add MMI and Aquatic-life-use to metrics:
mmi<-mmiprep[,c(1:7,9:17,30:31)]
write.csv(mmiprep, file="MMI.scores.and.values.csv")
write.csv(mmi, file="MMI_Scores.csv")
