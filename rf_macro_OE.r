# R code to make predictions of O/E for a new set of sites,;
#   based on a Random Forest predictive model:
# Modified from:
# J. Van Sickle, 02/25/10;
# Version 4.1 - Given same upgrade as model.predict.4.1
# Version 4.2 - Added option for out-of-bag predictions on CAL data. 12/30/10;
#
# Modified by: E. Kosnicki, 2013-Sep-05
# NOTE: this is mostly in a state of exploratory analysis
# Main focus was to create a cap for the O scores (max = to E) so that you 
# will never have an O/E > 1

## IMPORT AND LOAD PACKAGES
library(cluster); library(Hmisc);
library(randomForest);
## STEP 1 read in data:
predall<-read.csv("mydirtyFinal2.csv",row.names="Sample",header=T);
bugall<-read.csv("mytaxa.csv",row.names="SAMPLE",header=T)

# candvar<-c("lat","long","fieldtemp","ddcumlog","sandHUC_OS","highest_point__m_","yr1_precip","stream_power")
# ,"DECIDUOUS_"
candvar<-c("sandHUC_OS","sandHUC_A","sandHUC_CF","sandHUC_PD","drainage_area__m2_","basin_length__m_",
          "stream_length","highest_point__m_","basin_relief_ratio","entire_stream_gradient","sinuosity",
          "lat","long","fieldtemp","Q","wetWDave","wk2_precip","mon6_precip","ddcumlog","stream_power",
	"liveave","woodave","stream_order__Strahler_","GRASSLAND","DECIDUOUS_","EVERGREEN_","MIXED_FORE",
	"wetWave","yr1_temp","yr1_precip")

bugall.pa<-bugall;
bugall.pa[bugall.pa>0]<-1;  #This is the presence absence metrix

bug.log<-log10(bugall+1);    #Calculate the log of abundance data

row.names(bug.log)==row.names(predall); #to see if predicor and taxa data are aligned
#  bugall<-bugall[row.names(predall),];   #if they need to be aligned

# OK A BIG DECISION HERE IS TO KEEP THE PROGRAM AS IS AND RERUN FOR DIFFERENT REFERENCE LEVELS
# OR JUST ADD TO THIS PROGRAM ALL THE WAY THROUGH; MIGHT BE EASIER TO RUN SEPARATELY
# BUT WE STILL NEED TO ADD ADDITIONAL VARIABLES FOR TESTING AND SENSITIVITY PURPOSES IN STEP 7
# NOTE TO PROGRAM IN THE oRIG FORM DOES NOT NEED TO BE MODIFIED FOR PRIMARY MODEL

predcal<-predall[predall[,'Mdl704']=='C',];  #predictor data -calibration sites;
predcal<-predall[predall$Mdl704=='ref2'|predall$Mdl704=='C',]; #Pulling out more than one site;

pred.vld<-predall[substr(as.character(predall[,'Mdl704']),1,1)=='V',];  #predictor data - validation sites;

bugcal<-bugall[predall$Mdl704=='ref2'|predall$Mdl704=='C',]; #Bug Abundance matrix, calibration sites;
bugcal.pa<-bugall.pa[predall$Mdl704=='ref2'|predall$Mdl704=='C',]; #Bug presence/absence matrix, calibration sites;
bugcal.log<-bug.log[predall$Mdl704=='ref2'|predall$Mdl704=='C',]   #Bug log transformed matrix, calibartion sites


#bug.vld.pa<-bugall.pa[predall[,'Mdl704']=='V',]
bug.vld.pa<-bugall.pa[substr(as.character(predall[,'Mdl704']),1,1)=='V',];

bug.vld.log<-bug.log[predall[,'Mdl704']=='V',]

bug.misc.log<-bug.log[predall[,'Mdl704']=='V',]

## NOTE THE BELOW CODE EXTRACTS SUBSETS BASED ON ANY 'v-' IN THE VARIABLE EXPRESSION.
## bug.vld.pa<-bugall.pa[substr(as.character(predall[,'Mdl704']),1,1)=='V',]; #Bug presence/absence matrix, validation sites;

## STEP 2 -- CLUSTER ANALYSIS
psite.occ<-colSums(bugcal.pa)/dim(bugcal.pa)[[1]];

## nonrare.taxa<-names(psite.occ)[psite.occ>0.1]
nonrare.taxa<-names(psite.occ)[psite.occ>.1&psite.occ<.9]

bugcal.nonrare<- bugcal[,nonrare.taxa];
bugcal.pa.nonrare<- bugcal.pa[,nonrare.taxa];
bugcal.log.nonrare<- bugcal.log[,nonrare.taxa];

source("dapply.r");

# Option 1 -- Jaccard dissimilarity;
# function computes jaccard P/A dissimilarity between one site pair,  siti and sitj;
#input can be P/A or abundance data;

jaccfun<-function(siti,sitj) {
  shared<-sum((siti>0)&(sitj>0));
  uniquei<-sum((siti>0)&(sitj==0));
  uniquej<-sum((siti==0)&(sitj>0));
  1-(shared/(shared+uniquei+uniquej)); #return Jaccard dissimilarity;
} #end of function;

#compute Jaccard dissimilarities among all calibration sites;
dissim.jac<-dapply(bugcal.pa,1,bugcal.pa,1,jaccfun); #presence absence
dissim.jaca<-dapply(bugcal,1,bugcal,1,jaccfun); #abundance / density
dissim.jac2<-dapply(bugcal.pa.nonrare,1,bugcal.pa.nonrare,1,jaccfun); #log transformed
#Proceed to clustering;

sornfun<-function(siti,sitj) {
  shared<-sum((siti>0)&(sitj>0));
  uniquei<-sum((siti>0)&(sitj==0));
  uniquej<-sum((siti==0)&(sitj>0));
  1-(2*shared/(2*shared+uniquei+uniquej)); #return Sorenson dissimilarity;
} #end of function;

#compute Sorenson dissimilarities among all calibration sites;
dissim.sor<-dapply(bugcal.pa,1,bugcal.pa,1,sornfun);

# Clustering of calibration sites;
# Use flexible-Beta method, with Beta=-0.6;
#See R documentation on agnes() in "cluster" package;
# When using agnes() with Flexible Beta strategy, set Beta=(1-2*Alpha) in Lance-Williams formula;
#A single value for par.method value specifies alpha, so alpha=0.8 gives Beta=-0.6;

#######################################################################################
#########YOU WILL HAVE TO CHOOSE ONE OF THESE MODELS AT A TIME TO PROCEED FROM HERE

clus1<-agnes(x=dissim.jac,diss=T,method="flexible", par.method=0.8,keep.diss=F,keep.data=F);
clusa<-agnes(x=dissim.jac,diss=T,method="ward",keep.diss=F,keep.data=F);
clus2<-agnes(x=dissim.jac2,diss=T,method="flexible", par.method=0.8,keep.diss=F,keep.data=F);
clusab<-agnes(x=dissim.jac2,diss=T,method="ward",keep.diss=F,keep.data=F);
clusa<-agnes(x=dissim.jaca,diss=T,method="flexible", par.method=0.8,keep.diss=F,keep.data=F);
clus.sor<-agnes(x=dissim.sor,diss=T,method="flexible", par.method=0.4,keep.diss=F,keep.data=F);
```

## Various plots of cluster outcome. Leaf labels are row numbers of dissim matrix;
#that is, the order of sites in the calibration data set;
#plot the dendrogram;
plot(clus1,which.plots=2,labels=row.names(bugcal.pa),cex=.5);
plot(clusa,which.plots=2,labels=row.names(bugcal),cex=.5);
plot(clusab,which.plots=2,labels=row.names(bugcal.pa.nonrare),cex=.5);
plot(clus2,which.plots=2,labels=row.names(bugcal.pa.nonrare),cex=.5);
#######################;
#Pruning the dendrogram to create a small number of groups;
# level pruning can be done by specifying the number of groups (k parameter);
#Also can prune at a specified height. See cutree help;
#result is a vector of site group assignments;
#can repeat this process to generate several candidate groupings from a single dendrogram;

grps.6<-cutree(clus1,k=4); #vector of group assignments is in the order of sites in the clustered data;
table(grps.6); #count number of sites in each group;
cbind(row.names(predcal),grps.6); #list calibration sites and their group assignments;
#candidate site groups complete;
#Cluster analysis is complete;
###########;
#as a check, plot the geographic  locations of the site clusters;
# Can also use the "maps" package to draw a state or regional boundary;
plot(predcal$X_coord[grps.6==1], predcal$Y_coord[grps.6==1],col='black',
     type='p',xlim=c(-125,-118),ylim=c(42,46));     ## need to change the settings here
points(predcal$X_coord[grps.6==2], predcal$Y_coord[grps.6==2],col='red')
points(predcal$X_coord[grps.6==3], predcal$Y_coord[grps.6==3],col='green')
points(predcal$X_coord[grps.6==4], predcal$Y_coord[grps.6==4],col='orange')
points(predcal$X_coord[grps.6==5], predcal$Y_coord[grps.6==5],col='blue')
points(predcal$X_coord[grps.6==5], predcal$Y_coord[grps.6==5],col='purple')

#####################################;
# STEP 3 . BUILD RANDOM fOREST MODEL TO PREDICT GROUP MEMBERSHIP;
# First, put the group IDs and predictor variables into a single data frame;
rfdat<-data.frame(predcal[,candvar],Groups=grps.6);
rfdat$Groups<-factor(rfdat$Groups);

#Build RF model;
## NOTE !!!-- I have just used default settings for key fitting parameters, like "mtry". ;
# These choices could use some research, for the RIVPACS context! ;
rf.mod<-randomForest(Groups ~ ., data=rfdat, ntree=500, importance=TRUE, norm.votes=TRUE, keep.forest=TRUE);
print(rf.mod);   #Displays out-of-bag (OOB) error matrix. Columns are predicted class, rows are true class;
#various diagnostics;
varImpPlot(rf.mod);  #plots 2 measures of predictor importance;
#Note that all candidate predictors are used in a RF model.;
# Suggest re-running RF with only a few of the most important predictors;
# May give almost-equal performance;
######################################################################################################
#STEP 4 - Save the final RF predictive model, for future use;
# To specify the entire, final RF model, you need to store 4 things as an R object;
# 4.1) The site-by-taxa matrix of observed presence/absence at calibration sites (bugcal.pa, already available);
# 4.2) Specify the vector of final group membership assignments at calibration sites(grps.final);
grps.final<-grps.6;
# 4.3) Specify the final predictor variables ;
preds.final<-candvar;
# 4.4) The final RF model (rfmod);

# Save the model components together in a single .Rdata file.;
# Any R user can load this file, along with model.predict.RF.r, to make predictions from the model;
save(bugcal.pa, predcal, grps.final, preds.final, rf.mod, file='PrimaryRef.PA.redpred.Rdata');
#NOTE - Predcal is not needed to define the model, but is included so that users see the 
#       required format for predictor data;
#################################;

#Step 5 - Further checks on performance of the final, chosen model;

#Option 5.1 - Make predictions of E and O/E for calibration (reference) sites. Examine O/E statistics and plots;
# To do this, run the model.predict.RanFor.4.2 function, using the calibration data as the 'new' data;
# See Step 7 below, for more info on making predictions;
# Also see internal documentation of model.predict.Ran.For.4.2;
source("model.predict.RanFor.4.2.r");
# Two options, for calibration data;
# Option 1 - Set Cal.OOB=TRUE, for out-of-bag predictions (see Cutler et.al). Gives more realistic (larger) SD(O/E), appropriate for new data; 
# OE.assess.cal<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.000000000005,Cal.OOB=TRUE);
# Option 2 - Set Cal.OOB=FALSE, for in-bag predictions. Gives optomistically small SD(O/E), because RF models are tightly tuned to in-bag calibration data; 
OE.assess.cal<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.00000005,Cal.OOB=FALSE);
OE.assess.cal50<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.5,Cal.OOB=FALSE);

OE.assess.callog<-model.predict.RanFor.4.2(bugcal.log,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.00000005,Cal.OOB=FALSE);
OE.assess.calog50<-model.predict.RanFor.4.2(bugcal.log,grps.final,preds.final, ranfor.mod=rf.mod,prednew=predcal,bugnew=bugcal.pa,Pc=0.5,Cal.OOB=FALSE);

#look at other prediction results, for calibration sites;
names(OE.assess.cal);   #names of 2 components of the prediction results list;
head(OE.assess.cal$OE.scores); #data frame of O/E scores, 1st 5 rows;
head(OE.assess.cal$Capture.Probs); #predicted capture probabilties, 1st 5 rows;
head(OE.assess.cal$Group.Occurrence.Probs); #predicted group occurrence probabilities, 1st 5 rows;
#check distribution of Calibration-site O/E scores. Is it Normal?;
#plot a histogram and a Normal q-q plot;
par(mfrow=c(2,1));
hist(OE.assess.cal$OE.scores$OoverE,xlab="O/E");
qqnorm(OE.assess.cal$OE.scores$OoverE);

#scatterplot of O (on y-axis) vs E (on x-axis). See Pineiro et al. Ecol. Modelling 2008, 316-322, for this choice of axes;
par(mfrow=c(1,1));
plot(OE.assess.cal$OE.scores[,c('E','O')],xlab='Expected richness',ylab='Observed richness');
abline(0,1); #add a 1-1 line;
#########;

#Option 5.1.5 - Calculate replicate sampling SD of O/E, as a "perfect model" lower bound for SD(O/E) on calibration data;
# reference is Van Sickle et al. null model paper;
#first, compile the following function;
rep.sam.sd<-function(occprb.cal,Pc) {
  cal.cut<-(occprb.cal>=Pc); #TRUE/FALSE matrix denoting which taxa are above predicted P cutoff;
  #use occurrence probabilities only of taxa above the cutoff (cal.cut='TRUE');
  E.cal<-apply(occprb.cal*cal.cut,1,sum); #vector of predicted E ;
  # numerator of site-specific replicate sampling var. Result is a site vector
  RS.cal<-apply(occprb.cal*cal.cut,1,function(x)sum(x*(1-x)));
  SDRS<-sqrt(mean(RS.cal/(E.cal^2))); #replicate sampling SD is sqrt(mean(site-specific replicate sampling variances));
  print(' ',quote=F)
  print(' Replicate sampling SD of O/E: ',quote=F)
  print(SDRS,digits=4);  }; #end of function;
#Then execute the above function, using either in-bag or OOB predicted occurrence probs ('Capture probs') for the calibration data;
rep.sam.sd(occprb.cal=OE.assess.cal$Capture.Probs,Pc=0.5);
###########;

### Option 5.2 - Repeat Step 5.1, but this time use validation data. Check especially for model bias (mean(O/E) differs from 1.0);
OE.assess.vld<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=rf.mod,prednew=pred.vld,bugnew=bug.vld.pa,Pc=0.5,Cal.OOB=FALSE)  ;
OE.assess.vld$OE.scores;
## RIGHT NOW I CANNOT GET THE ABOVE FUNCTION TO RUN... SUBSCRIPT OUT OF BOUNDS ERROR, BUT VARIABLES SEEM
## TO BE CORRECT LENGTH

source("model.predict.RanFor.4.2.r");
load('SecondaryRef.PA.redpred.Rdata');

# User must supply a sample-by-taxa matrix of taxon abundance or else presence/absence (coded as 1 or 0), for all new samples;
# User must also supply a corresponding file of predictor data for those same samples;
# These 2 files should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
# Notes on format --
#   A) The sample ID column in both files should be read into R as a row name (see Step 1 examples).
#   B) Predictor data set -- Must include columns with the same names, units, etc.,
#        as the model's predictor variables. All other columns will be ignored;
#        Column order does not matter;
#        Predictions and calculations of O/E will be made only for those samples that have;
#        complete data for all model predictors.;
#   C)  Sample-by-taxa matrix. Can contain abundance or presence/absence (1 or 0). Missing or empty cells now allowed;
#       Sample ID's (row names) must match those of predictor data.
#       Any names for new taxa (column names) are acceptable, in any order;
#       HOWEVER - Only those new-data taxa names that match the names in the
#            calibration data can be use to calculate observed richness;
#            All other taxa (columns) in the new-data bug matrix are ignored;
#        To see a list of the calibration-taxa names, do:
names(bugcal.pa)[colSums(bugcal.pa)>0];
##########;
# Example predictions: For nonreference sites in the Oregon DEQ data set that are labeled "N_lc" (see Step 1);

pred.test<-predall[as.character(predall[,'Mdl704'])=='T',];  #predictor data - test sites;
bug.test.pa<-bugall.pa[as.character(predall[,'Mdl704'])=='T',]; #Bug presence/absence matrix, test sites;

#Alternatively, use the following code to just calculate values from all sample sites. At the moment you will
#have to export the Capture.Probs with the Obs and calculated ceiling values etc... in SAS...

pred.test<-predall
bug.test.pa<-bugall.pa
bug.test.log<-bug.log

#Drop all samples/sites that do not not have complete data for the model predictors;
pred.test<-pred.test[complete.cases(pred.test[,preds.final]),];
bug.test.pa<-bug.test.pa[row.names(pred.test),];
bug.test.log<-bug.test.log[row.names(pred.test),];

#makes predictions for test data;
OE.assess.test<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final,ranfor.mod=rf.mod, prednew=pred.test,bugnew=bug.test.pa,Pc=0.000000005,Cal.OOB=FALSE);
OE.assess.testlog<-model.predict.RanFor.4.2(bugcal.log,grps.final,preds.final,ranfor.mod=rf.mod, prednew=pred.test,bugnew=bug.test.log,Pc=0.0000005,Cal.OOB=FALSE);

# look at O/E scores, for all samples;
OE.assess.test$OE.scores;
head(OE.assess.test$Capture.Probs)
