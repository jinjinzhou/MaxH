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
pheno.nhw.plink <-  read.table(paste(data.folder,"CG10kNHWPhenoPCAs_4plink.txt",sep=""),header=T,na.strings = "-9")
## output continous covariates ## 
##"gender","age","agesquare","age*gender","agesquare*gender"

gender<-pheno.nhw.plink$gender
age <- pheno.nhw.plink$Age_Enroll


#qcovar <- data.frame(pheno.nhw.plink$FID,pheno.nhw.plink$IID,age,age^2,gender*age,gender*age*age,
#                pheno.nhw.plink$PC1,pheno.nhw.plink$PC2,pheno.nhw.plink$PC3,pheno.nhw.plink$PC4,
#                pheno.nhw.plink$PC5,pheno.nhw.plink$PC6,pheno.nhw.plink$PC7,pheno.nhw.plink$PC8,
#                pheno.nhw.plink$PC9,pheno.nhw.plink$PC10)
qcovar <- data.frame(pheno.nhw.plink$FID,pheno.nhw.plink$IID,age,gender,age^2,gender*age,gender*age*age,
                     pheno.nhw.plink$PC1,pheno.nhw.plink$PC2,pheno.nhw.plink$PC3)
covar <- data.frame(pheno.nhw.plink$FID,pheno.nhw.plink$IID,gender)

#write.table(qcovar,file="./datasets/qcovarNHW.txt",quote=F,col.names=F,row.names=F)
#write.table(covar,file="./datasets/covarNHW.txt",quote=F,col.names=F,row.names=F)

subtypes.pheno.names <- c("pctEmph_UpperLobes",
                    "pctEmph_LowerLobes",
                    "pctEmph_UL_LL_ratio", #strange distribution, removed from estimation
                    "pctEmph_Slicer",
                    "Slicer_15pctIn_Total",
                    "pctGasTrap_Slicer",
                    "pctEmph_UpperThird_Slicer",
                    "pctEmph_LowerThird_Slicer",
                    "Pi10_SRWA",
                    "WallAreaPct_seg",
                    "TLCpp_race_adjusted",
                    "FRCpp_race_adjusted",
                    "FEV1pp_utah",
                    "FVCpp_utah",
                    "FEV1_FVC_utah",
                    "BDR_pct_FEV1",
                    "BDR_pct_FVC")


index <- which(names(pheno.nhw.plink)%in%subtypes.pheno.names)


subtypes.phenotypes <- pheno.nhw.plink[,index]

subtypes.phenotypes$Slicer_15pctIn_Total <- -1*subtypes.phenotypes$Slicer_15pctIn_Total 

subtypes.phenotypes.out <- cbind(log(subtypes.phenotypes$pctEmph_UpperLobes),
                                 log(subtypes.phenotypes$pctEmph_LowerLobes),
                                 log(subtypes.phenotypes$pctEmph_UL_LL_ratio),
                                 log(subtypes.phenotypes$pctEmph_Slicer),
                                 log(subtypes.phenotypes$Slicer_15pctIn_Total),
                                 log(subtypes.phenotypes$pctGasTrap_Slicer),
                                 log(subtypes.phenotypes$pctEmph_UpperThird_Slicer),
                                 log(subtypes.phenotypes$pctEmph_LowerThird_Slicer),
                                 log(subtypes.phenotypes$Pi10_SRWA),
                                 subtypes.phenotypes$WallAreaPct_seg,
                                 subtypes.phenotypes$TLCpp_race_adjusted,
                                 log(subtypes.phenotypes$FRCpp_race_adjusted),
                                 subtypes.phenotypes$FEV1pp_utah,
                                 subtypes.phenotypes$FVCpp_utah,
                                 subtypes.phenotypes$FEV1_FVC_utah,
                                 subtypes.phenotypes$BDR_pct_FEV1,
                                 subtypes.phenotypes$BDR_pct_FVC)

subtypes.phenotypes.out <- data.frame(FID=pheno.nhw.plink$FID,IID=pheno.nhw.plink$IID,subtypes.phenotypes.out)

names(subtypes.phenotypes.out) <- c("FID","IID",
                                    "logpctEmph_UpperLobes",
                                    "logpctEmph_LowerLobes",
                                    "logpctEmph_UL_LL_ratio", #strange distribution, removed from estimation
                                    "logpctEmph_Slicer",
                                    "logSlicer_15pctIn_Total",
                                    "logpctGasTrap_Slicer",
                                    "logpctEmph_UpperThird_Slicer",
                                    "logpctEmph_LowerThird_Slicer",
                                    "logPi10_SRWA",
                                    "WallAreaPct_seg",
                                    "TLCpp_race_adjusted",
                                    "logFRCpp_race_adjusted",
                                    "FEV1pp_utah",
                                    "FVCpp_utah",
                                    "FEV1_FVC_utah",
                                    "BDR_pct_FEV1",
                                    "BDR_pct_FVC")

#write.table(subtypes.phenotypes.out,
#            file=paste(data.folder,"CG10kNHWPheno4Subtyping.txt",sep=""),quote=F,row.names=F,na = "-9")
                                
                                 
### calculate standardized residuals ###
nhw <- pheno.nhw.plink
FID=as.vector(nhw$FID)
IID=FID
gender=nhw$gender
age=nhw$Age_Enroll
agesquare=(nhw$Age_Enroll)^2
height=nhw$Height_CM
heightsquare=height^2

fev=nhw$FEV1pp_utah
FVCpp=nhw$FVCpp_utah
TLCpp<-nhw$TLCpp_race_adjusted
logFRCpp <- log(nhw$FRCpp_race_adjusted)
fev_fvc=nhw$FEV1_FVC_utah
logpctEmph_UL=log(nhw$pctEmph_UpperLobes)
logpctEmph_LL=log(nhw$pctEmph_LowerLobes)
logpctEmph_UL_LL_ratio = log(nhw$pctEmph_UL_LL_ratio)
logpctEmph_slicer=log(nhw$pctEmph_Slicer)
logpctEmph_UT=log(nhw$pctEmph_UpperThird_Slicer)
logpctEmph_LT=log(nhw$pctEmph_LowerThird_Slicer)
logpctGasTrap=log(nhw$pctGasTrap_Slicer)
logSlicer_15pctIn_Total=log(-1*nhw$Slicer_15pctIn_Total)
logPi10_SRWA=log(nhw$Pi10_SRWA)
WallAreaPct_seg = nhw$WallAreaPct_seg
BDR_pct_FEV1 = nhw$BDR_pct_FEV1
BDR_pct_FVC = nhw$BDR_pct_FVC


packyears=nhw$ATS_PackYears
packyearsquare=packyears^2
BMI=nhw$BMI

pc1<-nhw$PC1
pc2<-nhw$PC2
pc3<-nhw$PC3
pc4<-nhw$PC4
pc5<-nhw$PC5

pcs3<-paste(c("pc1","pc2","pc3"),collapse="+")
bas<-paste("gender","age","agesquare","age*gender","agesquare*gender",sep="+")


fev.fmla<-as.formula(paste("fev ~ ",paste(bas,pcs3,sep="+")))
fev.lm<-lm(fev.fmla)
fev_std_res=rstandard(fev.lm)

FVCpp.fmla<-as.formula(paste("FVCpp ~ ",paste(bas,pcs3,sep="+")))
FVCpp.lm<-lm(FVCpp.fmla)
FVCpp_std_res=rstandard(FVCpp.lm)

TLCpp.fmla<-as.formula(paste("TLCpp ~ ",paste(bas,pcs3,sep="+")))
TLCpp.lm<-lm(TLCpp.fmla)
TLCpp_std_res=rstandard(TLCpp.lm)

logFRCpp.fmla<-as.formula(paste("logFRCpp ~ ",paste(bas,pcs3,sep="+")))
logFRCpp.lm<-lm(logFRCpp.fmla)
logFRCpp_std_res=rstandard(logFRCpp.lm)

fev_fvc.fmla<-as.formula(paste("fev_fvc ~ ",paste(bas,pcs3,sep="+")))
fev_fvc.lm<-lm(fev_fvc.fmla)
fev_fvc_std_res=rstandard(fev_fvc.lm)

logpctEmph_UL.fmla<-as.formula(paste("logpctEmph_UL ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_UL.lm<-lm(logpctEmph_UL.fmla)
logpctEmph_UL_std_res=rstandard(logpctEmph_UL.lm)

logpctEmph_LL.fmla<-as.formula(paste("logpctEmph_LL ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_LL.lm<-lm(logpctEmph_LL.fmla)
logpctEmph_LL_std_res=rstandard(logpctEmph_LL.lm)

logpctEmph_UL_LL_ratio.fmla<-as.formula(paste("logpctEmph_UL_LL_ratio ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_UL_LL_ratio.lm<-lm(logpctEmph_UL_LL_ratio.fmla)
logpctEmph_UL_LL_ratio_std_res=rstandard(logpctEmph_UL_LL_ratio.lm)

logpctEmph_slicer.fmla<-as.formula(paste("logpctEmph_slicer ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_slicer.lm<-lm(logpctEmph_slicer.fmla)
logpctEmph_slicer_std_res=rstandard(logpctEmph_slicer.lm)

logpctEmph_UT.fmla<-as.formula(paste("logpctEmph_UT ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_UT.lm<-lm(logpctEmph_UT.fmla)
logpctEmph_UT_std_res=rstandard(logpctEmph_UT.lm)

logpctEmph_LT.fmla<-as.formula(paste("logpctEmph_LT ~ ",paste(bas,pcs3,sep="+")))
logpctEmph_LT.lm<-lm(logpctEmph_LT.fmla)
logpctEmph_LT_std_res=rstandard(logpctEmph_LT.lm)

logpctGasTrap.fmla<-as.formula(paste("logpctGasTrap ~ ",paste(bas,pcs3,sep="+")))
logpctGasTrap.lm<-lm(logpctGasTrap.fmla)
logpctGasTrap_std_res=rstandard(logpctGasTrap.lm)

logSlicer_15pctIn_Total.fmla<-as.formula(paste("logSlicer_15pctIn_Total ~ ",paste(bas,pcs3,sep="+")))
logSlicer_15pctIn_Total.lm<-lm(logSlicer_15pctIn_Total.fmla)
logSlicer_15pctIn_Total_std_res=rstandard(logSlicer_15pctIn_Total.lm)

logPi10_SRWA.fmla<-as.formula(paste("logPi10_SRWA ~ ",paste(bas,pcs3,sep="+")))
logPi10_SRWA.lm<-lm(logPi10_SRWA.fmla)
logPi10_SRWA_std_res=rstandard(logPi10_SRWA.lm)

WallAreaPct_seg.fmla<-as.formula(paste("WallAreaPct_seg ~ ",paste(bas,pcs3,sep="+")))
WallAreaPct_seg.lm<-lm(WallAreaPct_seg.fmla)
WallAreaPct_seg_std_res=rstandard(WallAreaPct_seg.lm)

BDR_pct_FEV1.fmla<-as.formula(paste("BDR_pct_FEV1 ~ ",paste(bas,pcs3,sep="+")))
BDR_pct_FEV1.lm<-lm(BDR_pct_FEV1.fmla)
BDR_pct_FEV1_std_res=rstandard(BDR_pct_FEV1.lm)


BDR_pct_FVC.fmla<-as.formula(paste("BDR_pct_FVC ~ ",paste(bas,pcs3,sep="+")))
BDR_pct_FVC.lm<-lm(BDR_pct_FVC.fmla)
BDR_pct_FVC_std_res=rstandard(BDR_pct_FVC.lm)

                                 
subtypes.res.out <- data.frame(FID=pheno.nhw.plink$FID,IID=pheno.nhw.plink$IID,
                               fev_std_res=fev_std_res,
                               FVCpp_std_res=FVCpp_std_res,
                               TLCpp_std_res=TLCpp_std_res,
                               logFRCpp_std_res=logFRCpp_std_res,
                               fev_fvc_std_res=fev_fvc_std_res,
                               logpctEmph_UL_std_res=logpctEmph_UL_std_res,
                               logpctEmph_LL_std_res=logpctEmph_LL_std_res,
                               logpctEmph_UL_LL_ratio_std_res=logpctEmph_UL_LL_ratio_std_res,
                               logpctEmph_slicer_std_res=logpctEmph_slicer_std_res,
                               logpctEmph_UT_std_res=logpctEmph_UT_std_res,
                               logpctEmph_LT_std_res=logpctEmph_LT_std_res,
                               logpctGasTrap_std_res=logpctGasTrap_std_res,
                               logSlicer_15pctIn_Total_std_res=logSlicer_15pctIn_Total_std_res,
                               logPi10_SRWA_std_res=logPi10_SRWA_std_res,
                               WallAreaPct_seg_std_res=WallAreaPct_seg_std_res,
                               BDR_pct_FEV1_std_res=BDR_pct_FEV1_std_res,
                               BDR_pct_FVC_std_res=BDR_pct_FVC_std_res)
write.table(subtypes.res.out,
            file=paste(data.folder,"CG10kNHWRes4Subtyping.txt",sep=""),
            quote=F,row.names=F,na = "-9")
