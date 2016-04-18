rm(list=ls())
setwd("/Users/jzhou/Documents/Research/MaxH/MaxH/")
options(na.action=na.exclude)

data.folder <- "./datasets/"

## this is original datasets ### 
pheno10000 <- read.table(paste(data.folder,"Final10000_dataset_02nov13.txt",sep=""),
                         header=T,sep="\t",na.strings = "-9")

race.index <-  which(names(pheno10000)=="race")
pheno10000.nhw <- pheno10000[which(pheno10000[,race.index]==1),]
pheno10000.aa <- pheno10000[which(pheno10000[,race.index]==2),]

## this is the dataset including PCs ##
pheno.nhw.plink <- read.table(paste(data.folder,"CG10kNHWPhenoPCAs_4plink.txt",sep=""),header=T,na.strings = "-9")
gender <- pheno.nhw.plink$gender
age <- pheno.nhw.plink$Age_Enroll

pcfile <- data.frame(FID=pheno.nhw.plink$FID,
                     pc1=pheno.nhw.plink$PC1,
                     pc2=pheno.nhw.plink$PC2,
                     pc3=pheno.nhw.plink$PC3,
                     pc4=pheno.nhw.plink$PC4,
                     pc5=pheno.nhw.plink$PC5)

subtypes.pheno.names <- 
  c("sid",
  "gender",
  "Age_Enroll",
  "Height_CM",
  "ATS_PackYears",
  "FEV1pp_utah",
  "FVCpp_utah",
  "FEV1_FVC_utah",
  "TLCpp_race_adjusted",
  "UpperThird_LowerThird_Slicer",
  "pctEmph_Slicer",
  "pctGasTrap_Slicer",
  "pctEmph_UpperThird_Slicer",
  "pctEmph_LowerThird_Slicer",
  "Pi10_SRWA",
  "WallAreaPct_seg",
  "BDR_pct_FEV1",
  "BDR_pct_FVC")

index <- which(names(pheno10000.nhw)%in%subtypes.pheno.names)
subtypes.phenotypes <- pheno10000.nhw[,index]

subtypes.phenotypes.out <- data.frame(FID=subtypes.phenotypes$sid,
                                      gender=subtypes.phenotypes$gender,
                                      age=subtypes.phenotypes$Age_Enroll,
                                      height=subtypes.phenotypes$Height_CM,
                                      packyears=subtypes.phenotypes$ATS_PackYears,
                                      fev=subtypes.phenotypes$FEV1pp_utah,
                                      FVCpp=subtypes.phenotypes$FVCpp_utah,
                                      fev_fvc=subtypes.phenotypes$FEV1_FVC_utah,
                                      TLCpp=subtypes.phenotypes$TLCpp_race_adjusted,
                                      logpctEmph=log(subtypes.phenotypes$pctEmph_Slicer),
                                      logpctGasTrap=log(subtypes.phenotypes$pctGasTrap_Slicer),
                                      logpctEmph_UT=log(subtypes.phenotypes$pctEmph_UpperThird_Slicer),
                                      logpctEmph_LT=log(subtypes.phenotypes$pctEmph_LowerThird_Slicer),
                                      logUT_LT=log(subtypes.phenotypes$UpperThird_LowerThird_Slicer),
                                      logPi10_SRWA=log(subtypes.phenotypes$Pi10_SRWA),
                                      WallAreaPctS=subtypes.phenotypes$WallAreaPct_seg,
                                      BDR_pct_FEV1=subtypes.phenotypes$BDR_pct_FEV1,
                                      BDR_pct_FVC=subtypes.phenotypes$BDR_pct_FVC)

all <- merge(subtypes.phenotypes.out,pcfile,by="FID")

rm(list= ls()[!(ls() %in% c('all'))])
attach(all)

#write.table(subtypes.phenotypes.out,
#            file=paste(data.folder,"CG10kNHWPheno4Subtyping.txt",sep=""),quote=F,row.names=F,na = "-9")
                                
### calculate standardized residuals ###
FID=as.vector(FID)
IID=FID
agesquare=age^2
heightsquare=height^2
packyearsquare=packyears^2

pcs3<-paste(c("pc1","pc2","pc3"),collapse="+")
bas<-paste("gender","age","height","packyears",sep="+")

fev.fmla<-as.formula(paste("fev ~ ",paste(bas,pcs3,sep="+")))
fev.lm<-lm(fev.fmla)
fev_std_res=rstandard(fev.lm)

FVCpp.fmla<-as.formula(paste("FVCpp ~ ",paste(bas,pcs3,sep="+")))
FVCpp.lm<-lm(FVCpp.fmla)
FVCpp_std_res=rstandard(FVCpp.lm)

TLCpp.fmla<-as.formula(paste("TLCpp ~ ",paste(bas,pcs3,sep="+")))
TLCpp.lm<-lm(TLCpp.fmla)
TLCpp_std_res=rstandard(TLCpp.lm)

fev_fvc.fmla<-as.formula(paste("fev_fvc ~ ",paste(bas,pcs3,sep="+")))
fev_fvc.lm<-lm(fev_fvc.fmla)
fev_fvc_std_res=rstandard(fev_fvc.lm)

logpctEmph.fmla<-as.formula(paste("logpctEmph ~ ",paste(bas,pcs3,sep="+")))
logpctEmph.lm<-lm(logpctEmph.fmla)
logpctEmph_std_res=rstandard(logpctEmph.lm)

logpctEmph_UT.fmla<-as.formula(paste("logpctEmph_UT ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_UT.lm<-lm(logpctEmph_UT.fmla)
logpctEmph_UT_std_res=rstandard(logpctEmph_UT.lm)

logpctEmph_LT.fmla<-as.formula(paste("logpctEmph_LT ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_LT.lm<-lm(logpctEmph_LT.fmla)
logpctEmph_LT_std_res=rstandard(logpctEmph_LT.lm)

logUT_LT.fmla <- as.formula(paste("logUT_LT ~ ",paste(bas,pcs3,sep="+")))
logUT_LT.lm <- lm(logUT_LT.fmla)
logUT_LT_std_res=rstandard(logUT_LT.lm)

logpctGasTrap.fmla<-as.formula(paste("logpctGasTrap ~ ",paste(bas,pcs3,sep="+")))
logpctGasTrap.lm<-lm(logpctGasTrap.fmla)
logpctGasTrap_std_res=rstandard(logpctGasTrap.lm)

logPi10_SRWA.fmla<-as.formula(paste("logPi10_SRWA ~ ",paste(bas,pcs3,sep="+")))
logPi10_SRWA.lm<-lm(logPi10_SRWA.fmla)
logPi10_SRWA_std_res=rstandard(logPi10_SRWA.lm)

WallAreaPctS.fmla<-as.formula(paste("WallAreaPctS ~ ",paste(bas,pcs3,sep="+")))
WallAreaPctS.lm<-lm(WallAreaPctS.fmla)
WallAreaPctS_std_res=rstandard(WallAreaPctS.lm)

BDR_pct_FEV1.fmla<-as.formula(paste("BDR_pct_FEV1 ~ ",paste(bas,pcs3,sep="+")))
BDR_pct_FEV1.lm<-lm(BDR_pct_FEV1.fmla)
BDR_pct_FEV1_std_res=rstandard(BDR_pct_FEV1.lm)

BDR_pct_FVC.fmla<-as.formula(paste("BDR_pct_FVC ~ ",paste(bas,pcs3,sep="+")))
BDR_pct_FVC.lm<-lm(BDR_pct_FVC.fmla)
BDR_pct_FVC_std_res=rstandard(BDR_pct_FVC.lm)

subtypes.res.out <- data.frame(FID=all$FID,
                               IID=FID,
                               fev_std_res=fev_std_res,
                               FVCpp_std_res=FVCpp_std_res,
                               TLCpp_std_res=TLCpp_std_res,
                               fev_fvc_std_res=fev_fvc_std_res,
                               logpctEmph_std_res=logpctEmph_std_res,
                               logpctEmph_UT_std_res=logpctEmph_UT_std_res,
                               logpctEmph_LT_std_res=logpctEmph_LT_std_res,
                               logUT_LT_std_res = logUT_LT_std_res,
                               logpctGasTrap_std_res=logpctGasTrap_std_res,
                               logPi10_SRWA_std_res=logPi10_SRWA_std_res,
                               WallAreaPctS_std_res=WallAreaPctS_std_res,
                               BDR_pct_FEV1_std_res=BDR_pct_FEV1_std_res,
                               BDR_pct_FVC_std_res=BDR_pct_FVC_std_res)


fam <- read.table("./datasets/CG10kNhwHg19Clean_v2_Mar2013.fam")
names(fam) <- c("FID","IID","MOM","DAD","SEX","CC")

subtypes.res.out <- merge(subtypes.res.out,fam,by="FID")
names(subtypes.res.out)[2] <- "IID"

write.table(subtypes.res.out[,seq(1,15)],
            file="./datasets/CG10kNHWRes4Subtyping.txt",
            quote=F,row.names=F,na = "-9")


