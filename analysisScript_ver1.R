library(tidyverse)
library(ggplot2)
library(gridExtra)
rm(list = ls())
source('D:/毕业论文temp/lmeFunction.R')
dataGbbArea <- read.csv('D:/毕业论文temp/gbbRawDataZip/gbbDataBrainArea.csv')[,-c(1:2)]
dataGbbThickness <- read.csv('D:/毕业论文temp/gbbRawDataZip/gbbDataBrainThickness.csv')[,-c(1:2)]
dataGbbVolume <- read.csv('D:/毕业论文temp/gbbRawDataZip/gbbDataBrainVolume.csv')[,-c(1:2)]
dataGbbMeanCurve <- read.csv('D:/毕业论文temp/gbbRawDataZip/gbbDataBrainMeanCurve.csv')[,-c(1:2)]
dataGbbGausCurve <- read.csv('D:/毕业论文temp/gbbRawDataZip/gbbDataBrainGausCurve.csv')[,-c(1:2)]
dataGbbArea$site <- 2
dataGbbThickness$site <- 2
dataGbbVolume$site <- 2
dataGbbMeanCurve$site <- 2
dataGbbGausCurve$site <- 2

dataSlimArea <- read.csv('D:/毕业论文temp/rawDataZipSlim/slimDataBrainArea.csv')[,-c(1,2,420)]
dataSlimThickness <- read.csv('D:/毕业论文temp/rawDataZipSlim/slimDataBrainThickness.csv')[,-c(1,2,420)]
dataSlimVolume <- read.csv('D:/毕业论文temp/rawDataZipSlim/slimDataBrainVolume.csv')[,-c(1,2,420)]
dataSlimMeanCurve <- read.csv('D:/毕业论文temp/rawDataZipSlim/slimDataBrainMeanCurve.csv')[,-c(1,2,420)]
dataSlimGausCurve <- read.csv('D:/毕业论文temp/rawDataZipSlim/slimDataBrainGausCurve.csv')[,-c(1,2,420)]
dataSlimArea$site <- 1
dataSlimThickness$site <- 1
dataSlimVolume$site <- 1
dataSlimMeanCurve$site <- 1
dataSlimGausCurve$site <- 1

calDataArea <- rbind(dataGbbArea, dataSlimArea)
calDataArea <- calDataArea[order(calDataArea$ID, calDataArea$label),]

calDataThickness <- rbind(dataGbbThickness, dataSlimThickness)
calDataThickness <- calDataThickness[order(calDataThickness$ID, calDataThickness$label),]

calDataVolume <- rbind(dataGbbVolume, dataSlimVolume)
calDataVolume <- calDataVolume[order(calDataVolume$ID, calDataVolume$label),]

calDataMeanCurve <- rbind(dataGbbMeanCurve, dataSlimMeanCurve)
calDataMeanCurve <- calDataMeanCurve[order(calDataMeanCurve$ID, calDataMeanCurve$label),]

calDataGausCurve <- rbind(dataGbbGausCurve, dataSlimGausCurve)
calDataGausCurve <- calDataGausCurve[order(calDataGausCurve$ID, calDataGausCurve$label),]

subcSlim <- read.csv('D:/毕业论文temp/rawDataZipSlim/subCorticalSlim.csv')
subcGbb <- read.csv('D:/毕业论文temp/gbbRawDataZip/subCorticalGbb.csv')
subcGbb$site <- 2
subcSlim$site <- 1
calDataSubc <- rbind(subcGbb, subcSlim)[,-1]
nrow(calDataArea[calDataArea$Sex == 1,])
nrow(calDataArea[calDataArea$Sex == 2,])
# male: 1 female: 2
printBaseInfo <- function(x){
  meanTemp <- x[x$seq!=0,]
  data.frame(info=c(nrow(x), length(unique(x$ID)), mean(meanTemp$seq), sd(meanTemp$seq), nrow(meanTemp[meanTemp$Sex==1,]), nrow(meanTemp[meanTemp$Sex==2,])),
             row.names = c("obs", "subjNum", "meanSeq", "stdSeq",'n_male','n_female'))
}
splitSubj <- function(x, timeLow, timeHigh){
  # x <- calDataArea
  # timeLow <- 1
  # timeHigh <- 1.2
  
  tempIndex <- rep(0, nrow(x))
  for (i in 1:length(tempIndex)) {
    if (x[i, "seq"]<=timeHigh&x[i, "seq"]>timeLow) {
      tempIndex[i] = i
      tempIndex[i-1] = i-1
      }
  }
  tempIndex <- tempIndex[tempIndex!=0]
  x <- x[tempIndex,]
  
  return(x)
}
splitSubj2 <- function(x, timeHigh){
  # x <- calDataArea
  # timeHigh <- 0.5
  
  tempIndex <- rep(0, nrow(x))
  for (i in 1:length(tempIndex)) {
    if (x[i, "seq"]<=timeHigh&x[i, "seq"]>0) {
      tempIndex[i] = i
      tempIndex[i-1] = i-1
    }
  }
  tempIndex <- tempIndex[tempIndex!=0]
  x <- x[tempIndex,]
  
  return(data.frame(x))
}

# basic informations 
tempInfo <- calDataArea[calDataArea$seq!=0,]
longInfo <- dataSlimArea[dataSlimArea$seq!=0,]

mean(longInfo$seq)
sd(longInfo$seq)
# mean and sd seq
mean(tempInfo$seq)
sd(tempInfo$seq)


calDataMSTemp <- calMs(calDataThickness, calDataGausCurve, calDataVolume, calDataMeanCurve, calDataArea)
calDataMS <- calDataMSTemp[["regionalMs"]]
calDataMSBehav <- calDataArea[, -c(50:409)]
calDataMS <- inner_join(calDataMS, calDataMSBehav, by = c("ID", "label"))
calDataMS <- select(calDataMS, colnames(calDataArea)[1:49], everything())
colnames(calDataMS)[50:409] <- paste0(colnames(calDataMS)[50:409], "_ROI_MS")
max(calDataArea$seq)
min(calDataArea$seq)
seq <- c(0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.9, 2.2,
         2.5, 99)
seqSlide <- data.frame(s=c(0, 0.3, 0.6, 0.9, 1.2, 1.5),
                       e=c(0.6, 0.9, 1.2, 1.5, 1.8, 99))
seqSlide2 <- data.frame(s=c(0, 0.4, 0.8, 1.2, 1.6, 2.0),
                       e=c(0.6, 1.0, 1.4, 1.8, 2.2, 99))
seqWindow <- c(0.8, 1, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 99)
seqR <- c(0, 0.9, 1.8, 99)
seqL <- nrow(seqSlide) ## length(seq)-1 or nrow(seqSlide) or nrow(seqSlide2)
measures <- c("Thickness", "Volume", "Area", "MS", 'GausCurve', 'MeanCurve')
measures1 <- c("Thickness", "Volume", "Area", "MS", "Subc", 'GausCurve', 'MeanCurve')
for (r in 1:seqL) {

  for (m in 1:length(measures1)) {
    M <- measures1[m]
    print(paste0('now calculating ', "reDf", r, M, '...', m+(r-1)*7, '/', length(measures1)*seqL))
    cmd <- paste0("reDf", r, M, " <- splitSubj(calData", M, ",", seqSlide[r,1], ",", seqSlide[r,2], ")")
    cmd2 <- paste0("reDf", r, M, " <- calLongAndCrossCov(", "reDf", r, M, ")")
    # cmd3 <- paste0("reDfValidation", r, M, " <- splitSubj(calData", M, ",", seqR[r], ",", seqR[r+1], ")")
    # cmd4 <- paste0("reDfValidation", r, M, " <- calLongAndCrossCov(", "reDfValidation", r, M, ")")
    eval(parse(text = c(cmd, cmd2)))
  }
  
}
## 0.5年为起始点，每经过0.2年计算一次
for (r in 1:length(seqWindow)) {
  
  for (m in 1:length(measures1)) {
    M <- measures1[m]
    print(paste0('now calculating ', m+(r-1)*7, '/', length(seqWindow)*length(measures1)))
    cmd <- paste0("reDfWindow.reDf", r, M, " <- splitSubj2(calData", M, ",", seqWindow[r], ")")
    cmd2 <- paste0("reDfWindow.reDf", r, M, " <- calLongAndCrossCov(", "reDfWindow.reDf", r, M, ")")
    eval(parse(text = c(cmd, cmd2)))
  }
  
}

measuresWindow <- ls()[grep('reDf.MS', ls())]
for (w in 1:length(measuresWindow)) {
  tmp <- get(measuresWindow[w])
  print(paste0('now plotting ', measuresWindow[w], '...', w, '/', length(measuresWindow)))
  demPlot(tmp)
  ggsave(paste0('D:/毕业论文temp/resultsPlt/timeWindow/demographInfo/', measuresWindow[w], '.pdf'), width = 4, height = 5)
}
demPlot(reDf1Area)
printBaseInfo(reDf1MS)
printBaseInfo(reDf2MS)
printBaseInfo(reDf3MS)
printBaseInfo(reDf4MS)
printBaseInfo(reDf5MS)
printBaseInfo(reDf6MS)

d1 <- demPlot(reDf1MS)
d2 <- demPlot(reDf2MS)
d3 <- demPlot(reDf3MS)
d4 <- demPlot(reDf4MS)
d5 <- demPlot(reDf5MS)
d6 <- demPlot(reDf6MS)
dall <- grid.arrange(d1, d2, d3, d4, d5, d6, ncol = 3)
ggsave('D:/毕业论文temp/resultsPlt/timeWindow/demographInfo/all.pdf', dall)

printBaseInfo(reDfValidation1Area)

demPlot(reDf2Area)
printBaseInfo(reDf2Area)
printBaseInfo(reDfValidation2Area)

demPlot(reDf6Area)
printBaseInfo(reDf5Area)
printBaseInfo(reDfValidation3Area)

demPlot(calDataArea)
printBaseInfo(calDataArea)
outputPath1 <- 'D:/毕业论文temp/outputPath/no_overlap/'
outputPath2 <- 'D:/毕业论文temp/outputPath/'
outputPathVal <- 'D:/毕业论文temp/outputPath/validation/'
# 不能直接套用！！！！要重新计算age_mean和centertime！！！！！！！！！！callong那个函数就是重新计算
# 计算皮层区域
for (l in 1:seqL) {
  
  for (m in 1:length(measures)) {
    M <- measures[m]
    cmd <- paste0("resSex", l, M, " <- lmePipelineSex(", "reDf", l, M, ", ", "outputPath1, ", l, ")")
    # cmd2 <- paste0("resSexVal", l, M, " <- lmePipelineSex(", "reDfValidation", l, M, ", ", "outputPathVal, ", l, ")")
    eval(parse(text = c(cmd)))
  }
  
}

for (l in 1:length(seqWindow)) {
  options (warn = -1)
  for (m in 1:length(measures)) {
    M <- measures[m]
    cmd <- paste0("windresSex", l, M, " <- lmePipelineSex(", "reDfWindow.reDf", l, M, ", ", "outputPath2, ", l, ")")
    eval(parse(text = cmd))
  }
  
}
# 计算皮层下区域
fitStoreTime <- list()
fitStoreTimeBySex <- list()
dataTmpAll <- c('reDf1Subc', 'reDf2Subc', 'reDf3Subc', 'calDataSubc')
for (l in 1:length(dataTmpAll)) {
  dataTmp <- get(dataTmpAll[l])
  cmd <- paste0('resDfSubc',l,' <- data.frame(array(0,c(14,8)))')
  eval(parse(text = cmd))
  for (r in 1:14) {
    
    if(l<4){
      print(paste0(dataTmpAll[l],'  region: ',colnames(dataTmp)[56+r]))
      l1 <- as.formula(paste0(colnames(dataTmp)[56+r],'~Sex+age_mean+Sex*centertime+site+(1|ID)'))
      l2 <- as.formula(paste0(colnames(dataTmp)[56+r],'~Sex+age_mean+centertime+site+(1|ID)'))
    }
    else{
      print(paste0(dataTmpAll[l],'  region: ',colnames(dataTmp)[58+r]))
      l1 <- as.formula(paste0(colnames(dataTmp)[58+r],'~Sex+age_mean+Sex*centertime+site+(1|ID)'))
      l2 <- as.formula(paste0(colnames(dataTmp)[58+r],'~Sex+age_mean+centertime+site+(1|ID)'))
    }

    fitStoreTimeBySex[[r]] <- lmer(l1, data = dataTmp)
    fitStoreTime[[r]] <- lmer(l2, data = dataTmp)
  }
  pvectTimeBySex <- do.call(cbind,lapply(fitStoreTimeBySex, function(x) summary(x)$coefficients['Sex:centertime',]))
  pvectTime <- do.call(cbind,lapply(fitStoreTime, function(x) summary(x)$coefficients['centertime',]))
  resTmp <- get(paste0('resDfSubc',l))
  rownames(resTmp) <- colnames(dataTmp)[57:70]
  if(l==4) rownames(resTmp) <- colnames(dataTmp)[59:72]
  colnames(resTmp) <- c('pval.timebysex','tval.timebysex','pval.time','tval.time',
                           'pval.timebysex.corr','tval.timebysex.corr','pval.time.corr','tval.time.corr')
  resTmp[,'pval.timebysex'] <- pvectTimeBySex[5,]
  resTmp[,'tval.timebysex'] <- pvectTimeBySex[4,]
  resTmp[,'pval.time'] <- pvectTime[5,]
  resTmp[,'tval.time'] <- pvectTime[4,]
  resTmp[,'pval.timebysex.corr'] <- p.adjust(pvectTimeBySex[5,], method = 'fdr')
  resTmp[,'pval.time.corr'] <- p.adjust(pvectTime[5,], method = 'fdr')
  for (c in 1:14) {
    if(resTmp[,'tval.timebysex'][c] < 0) {
      resTmp[,'tval.timebysex.corr'][c] <- -qt(resTmp[,'pval.timebysex.corr'][c]/2, pvectTimeBySex[3,][c], lower.tail = F)
    }else{
      resTmp[,'tval.timebysex.corr'][c] <- qt(resTmp[,'pval.timebysex.corr'][c]/2, pvectTimeBySex[3,][c], lower.tail = F)
    }
    if(resTmp[,'tval.time'][c] < 0) {
      resTmp[,'tval.time.corr'][c] <- -qt(resTmp[,'pval.time.corr'][c]/2, pvectTime[3,][c], lower.tail = F)
    }else{
      resTmp[,'tval.time.corr'][c] <- qt(resTmp[,'pval.time.corr'][c]/2, pvectTime[3,][c], lower.tail = F)
    }
  }
  cmd2 <- paste0('resDfSubc',l,' <- resTmp')
  eval(parse(text = cmd2))
}

# 计算全时段发育
measuresAll <- c('calDataArea', 'calDataVolume', 'calDataThickness', 'calDataMS', 'calDataMeanCurve', 'calDataGausCurve')
for (m in 1:length(measuresAll)) {
  M <- measuresAll[m]
  cmd <- paste0("resSex4", M, " <- lmePipelineSex(", M, ", outputPath2,  4)")
  eval(parse(text = cmd))
}



# check split results
corValTest <- matrix(0, length(measures)*3, 2)
colnames(corValTest) <- c("corr", "pval")
rownames(corValTest) <- paste0(rep(measures, seqL), "-", rep(measures, seqL))
for (l in 1:seqL) {
  
  for (m in 1:length(measures)) {
    M <- measures[m]
    cmd <- paste0("cor.test(resSex", l, M, "$Sex.tvalue, resSexVal", l, M, "$Sex.tvalue)")
    rTemp <- eval(parse(text = cmd))
    r <- rTemp[["estimate"]][["cor"]]
    p <- rTemp[["p.value"]]
    corValTest[m+(l-1)*length(measures), 1] <- r
    corValTest[m+(l-1)*length(measures), 2] <- p
  }
  
}

# time or Sex

for (p in 1:seqL) {
  
  for (m in 1:length(measures)) {
    
    M <- measures[m]
    res <- paste0("resSex", p, M)
    resVal <- paste0("resSexVal", p, M)
    cmd <- paste0(res, "$rank <- 1:360")
    cmd2 <- paste0("sexRes", p, M, "Plt <- ", res, "[", res, "$Sex.pvalue<=0.05,]")
    # cmd3 <- paste0(resVal, "$rank <- 1:360")
    # cmd4 <- paste0("sexResVal", p, M, "Plt <- ", resVal, "[", resVal, "$time.pvalue<=0.05,]")
    eval(parse(text = c(cmd, cmd2)))
    
  }
}


for (p in 1:length(seqWindow)) {
  
  for (m in 1:length(measures)) {
    
    M <- measures[m]
    res <- paste0("windresSex", p, M)
    cmd <- paste0(res, "$rank <- 1:360")
    cmd2 <- paste0("windresTimeThreshold", p, M, "Plt <- ", res, "[", res, "$time.pvalue<=0.05,]")
    eval(parse(text = c(cmd, cmd2)))
    
  }
}


# plot test
library(ggseg)
library(ggseg3d)
library(ggsegExtra)
library(ggsegGlasser)

# 皮层图
## timeWindow plot ## 这里画图有要手动更改的地方
for (m in 1:length(measures)) {
  
  for (p in 1:length(seqWindow)) {
    
    measure <- measures[m]
    index <- paste0("windresTimeThreshold", p, measure, "Plt")
    pltIndex <- get(index)
    SigRegionTvaluePlot <- ggseg (pltIndex,atlas = glasser, 
                                  mapping = aes(fill = time.tvalue.Corr),  #change: time.tvalue.Corr or Sex.tvalue.Corr
                                  show.legend = T, 
                                  position = "stacked") +
      scale_fill_gradient2(midpoint = 0,mid = "grey",low = "purple",high = "yellow")+
      # scale_fill_gradient(low="purple",high="yellow")+
      # scale_fill_gradient(low="purple",high="yellow",limits=c(-9,3))+
      theme_void()+
      theme(text=element_text(size=16,  family="serif"),
            legend.position = 'none')
    
    plot(SigRegionTvaluePlot)
    if(dir.exists(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure))==F) dir.create(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure))
    ggsave(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure, "/", index, ".pdf"), SigRegionTvaluePlot)
    
    eval(parse(text=paste0('plt', p, ' <- SigRegionTvaluePlot')))
  }
  check <- grid.arrange(plt2, plt3, plt4, plt5, plt6, plt7, plt8, plt9, 
                        plt10, plt11, plt12, plt13, ncol = 4) ## 这里要根据seqWindow手动更改下,还是代码功力不足啊
  ggsave(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure, "/", "All.pdf"), check)

}

## 原来的plot
for (m in 1:length(measures)) {
  
  for (p in 1:seqL) {
    
    measure <- measures[m]
    index <- paste0("sexRes", p, measure, "Plt")
    pltIndex <- get(index)
    SigRegionTvaluePlot <- ggseg (pltIndex,atlas = glasser, 
                                  mapping = aes(fill = Sex.tvalue.Corr),  #change: time.tvalue.Corr or Sex.tvalue.Corr
                                  show.legend = T, 
                                  position = "stacked") +
      scale_fill_gradient2(midpoint = 0,mid = "grey",low = "purple",high = "yellow")+
      # scale_fill_gradient(low="purple",high="yellow")+
      # scale_fill_gradient(low="purple",high="yellow",limits=c(-9,3))+
      labs(title = index, fill = "t-value")+
      theme_void()+
      theme(text=element_text(size=16,  family="serif"),
            legend.position = 'none')
    
    plot(SigRegionTvaluePlot)
    if(dir.exists(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure))==F) dir.create(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure))
    ggsave(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure, "/", index, "timebysex.pdf"), SigRegionTvaluePlot)
    
    eval(parse(text=paste0('plt', p, ' <- SigRegionTvaluePlot')))
  }
  check <- grid.arrange(plt1, plt2, plt3, plt4, plt5, plt6, ncol = 3)
  ggsave(paste0("D:/毕业论文temp/resultsPlt/timeWindow/", measure, "/", "timebysexAll.pdf"), check)
}

# 皮层下区域图
indexTmp <- c('tval.time.corr', 'tval.timebysex.corr')
options (warn = -1)
for (l in 1:4) {
  dfTmp <- get(paste0('resDfSubc',l))
  # if (l == 4 ) dfTmp <- resDfSubc4
  dfTmp$label <- rownames(dfTmp)
  # 运行两次，因为有些区域两个点，为什么不能一次替换完呢？
  dfTmp$label <- str_replace(dfTmp$label, '\\.', '-')
  dfTmp$label <- str_replace(dfTmp$label, '\\.', '-')
  dfTmp$label <- str_replace(dfTmp$label, 'Thalamus', 'Thalamus-Proper')
  dfTmp <- dfTmp[dfTmp$pval.timebysex<0.05,]  ##change!!: pval.time.corr or pval.timebysex.corr
  plt <- ggseg(dfTmp, atlas = aseg,
        mapping = aes(fill = tval.time.corr), show.legend = T, view = 'coronal')+
    scale_fill_gradient2(midpoint = 0,mid = "grey",low = "purple",high = "yellow")
    
  
  if(dir.exists(paste0("D:/毕业论文temp/resultsPlt/subCortical"))==F) dir.create(paste0("D:/毕业论文temp/resultsPlt/subCortical"))
  ggsave(paste0("D:/毕业论文temp/resultsPlt/subCortical/subcortical", l, "timebysex.pdf"), plt) ##change!!:time or timebysex
}

##全时段发育图
resAll <- ls()[grep('Sex4cal', ls())]
for (i in 1:length(resAll)) {
  options (warn = -1)
  print(paste0('now plotting ',resAll[i]))
  dfTmp <- get(resAll[i])
  dfTmp <- dfTmp[dfTmp$Sex.pvalue<0.05,]
  plt <- ggseg(dfTmp, atlas = glasser,
               mapping = aes(fill = Sex.tvalue.Corr), show.legend = T, position = "stacked", color = 'black')+
    scale_fill_gradient2(midpoint = 0,mid = "white",low = "purple",high = "yellow")+
    labs(title = resAll[i], fill = "t-value")+
    theme_void()+
    theme(text=element_text(size=12))
  if (nrow(dfTmp) <= 1) {
    plt <- ggseg(dfTmp, atlas = glasser,
                 mapping = aes(fill = Sex.tvalue.Corr), show.legend = T, position = "stacked", color = 'black')+
      labs(title = resAll[i], fill = "t-value")+
      theme_void()+
      theme(text=element_text(size=12))
  }
  
  if(dir.exists(paste0("D:/毕业论文temp/resultsPlt/AllTime"))==F) dir.create(paste0("D:/毕业论文temp/resultsPlt/AllTime"))
  ggsave(paste0("D:/毕业论文temp/resultsPlt/AllTime/", resAll[i], "_timebysex.pdf"), plt)
}

## mediation effect plot
library(lmerTest)
library(ggplot2)
medDf <- calDataMS[, c(1,3,419,7,411,232,119)]
medDf[,4:7] <- scale(medDf[,4:7])
fit <- lmer(rh_R_V6_ROI_MS~age_mean+centertime*Sex+site+(1|ID), data = medDf)
fit2 <- lmer(lh_L_8BL_ROI_MS~age_mean+centertime*Sex+site+(1|ID), data = medDf)
coef <- summary(fit)$coefficients[,1]
coef2 <- summary(fit2)$coefficients[,1]
lineMV6 <- 1*coef['(Intercept)']+medDf$centertime*coef['centertime']+1*coef['Sex']+medDf$centertime*1*coef['centertime:Sex']+coef['site']
lineFV6 <- 1*coef['(Intercept)']+medDf$centertime*coef['centertime']+2*coef['Sex']+medDf$centertime*2*coef['centertime:Sex']+coef['site']
lineM8BL <- 1*coef2['(Intercept)']+medDf$centertime*coef2['centertime']+1*coef2['Sex']+medDf$centertime*1*coef2['centertime:Sex']+coef2['site']
lineF8BL <- 1*coef2['(Intercept)']+medDf$centertime*coef2['centertime']+2*coef2['Sex']+medDf$centertime*2*coef2['centertime:Sex']+coef2['site']

pltDf <- medDf
pltDf$Sex <- factor(pltDf$Sex, levels = c('2','1'), labels = c('Female', 'Male'))
pltV6 <- ggplot(pltDf,aes(centertime, rh_R_V6_ROI_MS))+
  geom_point(aes(color = Sex), alpha=0.5)+
  geom_line(aes(group = ID, color = Sex), alpha=0.5)+
  geom_line(aes(y=lineFV6),color = '#F8766D',size = 2)+
  geom_line(aes(y=lineMV6),color = '#00BFC4',size = 2)+
  xlab('Scan time') +
  ylab('Rh_V6_MS')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12))

plt8BL <- ggplot(pltDf,aes(centertime, lh_L_8BL_ROI_MS))+
  geom_point(aes(color = Sex), alpha=0.5)+
  geom_line(aes(group = ID, color = Sex), alpha=0.5)+
  geom_line(aes(y=lineF8BL),color = '#F8766D',size = 2)+
  geom_line(aes(y=lineM8BL),color = '#00BFC4',size = 2)+
  xlab('Scan time') +
  ylab('Lh_8BL_MS')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12))
ggsave('D:/毕业论文temp/resultsPlt/medPltV6MS.pdf',pltV6,height = 3, width = 4)
ggsave('D:/毕业论文temp/resultsPlt/medPlt8BLMS.pdf',plt8BL,height = 3, width = 4)

# 3d plot

ggseg3d(sexRes1AreaPlt, atlas = glasser_3d, colour = 'time.tvalue.Corr')%>%
  pan_camera("right lateral")

##  Cell Type–Specific Enrichment Analysis (CSEA) correct ##
library(gdata)
library(pSI)
library(pSI.data)
plsRes <- read.csv('D:/毕业论文temp/PLS/PLS1_geneWeights_ZMSSexbyTime.csv')
chromosome <- read.csv('D:/毕业论文temp/PLS/chromosomeInfo.csv')[,-1]
plsRes <- inner_join(plsRes, chromosome, by='Gene')
plsRes$ZtoP <- 0
plsRes$p.corr <- 0
Zcorr <- function(x){
  for (p in 1:nrow(x)) {
    if(x[p, 'Z']<=0) x[p, 'ZtoP'] <- 2*pnorm(x[p, 'Z'])
    if(x[p, 'Z']>0) x[p, 'ZtoP'] <- 2*pnorm(x[p, 'Z'], lower.tail = F)
    
    x$p.corr <- p.adjust(x$ZtoP, method = 'fdr')
  }
  ret <- x[x$p.corr<=0.01,]
  return(ret)
}

plsResCorr <- Zcorr(plsRes)
write.csv(plsResCorr, 'D:/毕业论文temp/PLS/pls_chromosome_corr.csv')

## GO and KEGG



## chromosome median test
# x <- '22'
chr_median <- sapply(as.character(unique(chromosome$Cluster)),function(x) median(rank(plsRes$Z)[which(plsRes$Cluster==x)]))-median(rank(plsRes$Z))
chr_se <- sapply(as.character(unique(chromosome$Cluster)),function(x) sd(rank(plsRes$Z)[which(plsRes$Cluster==x)])/sqrt(length(which(plsRes$Cluster==x))))
chr_n_genes <- sapply(as.character(unique(chromosome$Cluster)),function(x) length(which(plsRes$Cluster==x)))
names(chr_n_genes)
chr_perm <- data.frame(array(0, c(1001, 24)))
colnames(chr_perm) <- names(chr_n_genes)
set.seed(1)
for (c in 1:length(names(chr_n_genes))) {
  for (p in 1:1000) {
    print(paste0('chr: ', names(chr_n_genes)[c], ' perm: ', p))
    chr_perm[p+1, c] <- median(sample(rank(plsRes$Z),chr_n_genes[c],replace = F))-median(rank(plsRes$Z))
  }
}

for (i in 1:ncol(chr_perm)) {
  tmp <- chr_perm[2:1001,i]
  if(chr_median[i]>=0) chr_perm[1, i] <- length(tmp[tmp>=chr_median[i]])/1000
  if(chr_median[i]<0) chr_perm[1, i] <- length(tmp[tmp<=chr_median[i]])/1000
    
}
chr_perm[1,] <- p.adjust(chr_perm[1,], method = 'fdr')
chr_plot <- data.frame(chr=as.factor(names(chr_median)),
                       median=chr_median, # medians to plot
                       se_up=chr_median+chr_se,
                       se_down=chr_median-chr_se,
                       sign=t(ifelse(chr_perm[1,]<0.05,'*','')) # FDR-corrected p-values
)
colnames(chr_plot)[5] <- 'sign'
check <- chr_plot[order(chr_plot$median, decreasing = F),1]
chr_plot$chr <- factor(chr_plot$chr, levels = chr_plot[order(chr_plot$median),1])

pltChr <- ggplot(chr_plot,aes(y=median,x=chr))+
  geom_bar(data=chr_plot, aes(y=median,x=chr,fill=(chr_plot$median+1+abs(min(chr_plot$median)))),stat = 'identity') +
  geom_text(aes(label=sign),hjust=ifelse(chr_plot$median<0,2,-2),vjust=0.75,cex=10) +
  scale_fill_gradientn(colours = c("blue","white","red"),limits=c(0,2*abs(min(chr_plot$median))))+
  geom_pointrange(aes(ymin=se_down,ymax=se_up),colour="black")+
  coord_flip()+
  geom_hline(yintercept=0)+
  theme_minimal() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())+
  xlab('Chromosome') +
  ylab('Median Rank')
ggsave('D:/毕业论文temp/resultsPlt/enrichChr.pdf', pltChr)
# Decoding

library(neurobase)
atlas <- readnii("D:/毕业论文temp/decoding/HCP-MMP1_2mm.nii.gz")


niiName <- ls()[grep("Plt", ls())]
for (n in 1:length(niiName)) {
  
  atlascp <- atlas
  atlascp[atlascp] <- 0
  nii <- get(niiName[n])
  for (r in 1:nrow(nii)) {
    
    atlascp[atlas==nii$rank[r]] <- nii$time.tvalue.Corr[r]
    
  }
  
  writenii(atlascp, paste0("D:/毕业论文temp/decoding/", niiName[n], ".nii"))
}

# development change

reDf3VolumeFM <- reDf3Volume[reDf3Volume$Sex==-1,]
reDf3VolumeM <- reDf3Volume[reDf3Volume$Sex==1,]

plotDf <- reDf3Volume[,c(1:7, 91)]
plotDf[,4:8] <- scale(plotDf[,4:8])
plotDf$Sex1 <- 0.345+1*0.07+0.014*reDf3Volume$centertime-0.07*1*reDf3Volume$centertime
plotDf$Sex0 <- 0.345-0.014*reDf3Volume$centertime+0.07*reDf3Volume$centertime
ggplot(plotDf, aes(centertime, lh_L_7AL_ROI_volume))+
  geom_point(alpha = 0.3)+
  geom_line(aes(group = ID), alpha = 0.3, size = 1)+
  geom_line(aes(y=Sex1), color = "red", size = 2)+
  geom_line(aes(y=Sex0), color = "yellow", size = 2)

#
for (t in 50:410) {
  time1 <- reDf3Volume[reDf3Volume$label=="Time1", t]
  time2 <- reDf3Volume[reDf3Volume$label=="Time2", t]
  AA <- t.test(time1, time2)
  print(AA[["p.value"]])
}
subTemp <- calDataArea
calMI <- function(data){
  # data <- calDataArea
  subID <- unique(data$ID)
  subNum <- length(subID)
  brainIndex <- grep("ROI", colnames(subTemp))
  indexBeg <- brainIndex[1]
  indexEnd <- brainIndex[length(brainIndex)]
  
  df <- matrix(0, subNum, 362)
  colnames(df) <- c("ID", "Sex", colnames(data[brainIndex]))
  
  for (s in 1:subNum) {
    
    
    subTemp <- data[data$ID==subID[s],]
    
    timePointNum <- nrow(subTemp)
    M <- timePointNum-1
    MI <- rep(0, M)
    
    for (r in 1:length(brainIndex)) {
      
      print(paste0("sub: ", s, "  region: ", r))
      for (m in 1:(timePointNum-1)) {
        MI[m] <- (subTemp[m+1, brainIndex[r]]-subTemp[m, brainIndex[r]])/subTemp[m+1, "seq"]
      }
      MIMean <- sum(MI)/M
      df[s, r+2] <- MIMean
    }

    df[s, 1] <- subID[s]
    df[s, 2] <- unique(subTemp$Sex)
    
  }
  df <- data.frame(df)
  df <- na.omit(df)
}

calDataAreaMI <- calMI(calDataArea)
calDataThicknessMI <- calMI(calDataThickness)
calDataVolumeMI <- calMI(calDataVolume)
calDataMSMI <- calMI(calDataMS)
calDataMCMI <- calMI(calDataMeanCurve)
calDataGCMI <- calMI(calDataGausCurve)
## select baseline index
idSelected <- unique(calDataAreaMI$ID)
getBaselineInfo <- function(x){
  
  dfTmp <- data.frame(array(0,c(length(idSelected), 362)))
  # x <- calDataArea
  brainIndex <- grep("ROI", colnames(x))
  colnames(dfTmp) <- c('ID', 'Sex', paste0(colnames(data[brainIndex]), '_base'))
  for (s in 1:length(idSelected)) {
    print(paste0(idSelected[s], '   ', s, '/', length(idSelected)))
    dfTmp[s, 1] <- idSelected[s]
    dfTmp[s, 2] <- unique(x[x$ID==idSelected[s], 'Sex'])
    dfTmp[s, 3:362] <- x[x$ID==idSelected[s]&x$label=='Time1', brainIndex]
  }
  return(dfTmp)
}
calDataAreaBase <- getBaselineInfo(calDataArea)
calDataThicknessBase <- getBaselineInfo(calDataThickness)
calDataVolumeBase <- getBaselineInfo(calDataVolume)
calDataMSBase <- getBaselineInfo(calDataMS)
calDataMCBase <- getBaselineInfo(calDataMeanCurve)
calDataGCBase <- getBaselineInfo(calDataGausCurve)



write.csv(calDataAreaBase, "D:/毕业论文temp/outputPath/MI/calDataAreaBase.csv")
write.csv(calDataThicknessBase, "D:/毕业论文temp/outputPath/MI/calDataThicknessBase.csv")
write.csv(calDataVolumeBase, "D:/毕业论文temp/outputPath/MI/calDataVolumeBase.csv")
write.csv(calDataMSBase, "D:/毕业论文temp/outputPath/MI/calDataMSBase.csv")
write.csv(calDataMCBase, "D:/毕业论文temp/outputPath/MI/calDataMCBase.csv")
write.csv(calDataGCBase, "D:/毕业论文temp/outputPath/MI/calDataGCBase.csv")

write.csv(calDataAreaMI, "D:/毕业论文temp/outputPath/MI/calDataAreaMI.csv")
write.csv(calDataThicknessMI, "D:/毕业论文temp/outputPath/MI/calDataThicknessMI.csv")
write.csv(calDataVolumeMI, "D:/毕业论文temp/outputPath/MI/calDataVolumeMI.csv")
write.csv(calDataMSMI, "D:/毕业论文temp/outputPath/MI/calDataMSMI.csv")
write.csv(calDataMCMI, "D:/毕业论文temp/outputPath/MI/calDataMCMI.csv")
write.csv(calDataGCMI, "D:/毕业论文temp/outputPath/MI/calDataGCMI.csv")

##  
list.files('D:/毕业论文temp/outputPath/ML')
msWeight <- read.table('D:/毕业论文temp/outputPath/ML/calDataMSMIWeight.txt')
areaWeight <- read.table('D:/毕业论文temp/outputPath/ML/calDataAreaMIWeight.txt')
volumeWeight <- read.table('D:/毕业论文temp/outputPath/ML/calDataVolumeMIWeight.txt')
thicknessWeight <- read.table('D:/毕业论文temp/outputPath/ML/calDataThicknessMIWeight.txt')
mcWeight <- read.table('D:/毕业论文temp/outputPath/ML/calDataMCMIWeight.txt')
gcWeight <- read.table('D:/毕业论文temp/outputPath/ML/calDataGCMIWeight.txt')
msPlt <- cbind(msWeight, resSex4calDataMS$Sex.tvalue)
areaPlt <- cbind(areaWeight, resSex4calDataArea$Sex.tvalue)
volumePlt <- cbind(volumeWeight, resSex4calDataVolume$Sex.tvalue)
thicknessPlt <- cbind(thicknessWeight, resSex4calDataThickness$Sex.tvalue)
cor.test(msWeight$V1, resSex4calDataMS$Sex.tvalue)
cor.test(areaWeight$V1, resSex4calDataArea$Sex.tvalue)
cor.test(volumeWeight$V1, resSex4calDataVolume$Sex.tvalue)
cor.test(thicknessWeight$V1, resSex4calDataThickness$Sex.tvalue)


ggplot(msPlt, aes(V1, resSex4calDataMS$Sex.tvalue))+
  geom_point()+
  stat_smooth(method = 'lm', size = 2)+
  labs(x = "LR weight")+
  annotate("text", label = "r = 0.406, p < 0.0001", x = 0.6, y = -4, size = 6, fontface = 'italic')
ggsave('D:/毕业论文temp/resultsPlt/msWeight_tmap4.pdf', width = 4, height = 5)

## bootstrap results t-test
library(reshape2)
library(ggpubr)
bootMI <- read.table('D:/毕业论文temp/outputPath/ML/bootInter.txt')
accMI <- read.table('D:/毕业论文temp/outputPath/ML/dataInterAcc.txt')
bootBase <- read.table('D:/毕业论文temp/outputPath/ML/bootBase.txt')
accBase <- read.table('D:/毕业论文temp/outputPath/ML/dataBaseAcc.txt')
bootEx <- read.table('D:/毕业论文temp/outputPath/ML/bootEx.txt')
accEx <- read.table('D:/毕业论文temp/outputPath/ML/dataExAcc.txt')
bootMI$V1 <- as.factor(bootMI$V1)
bootBase$V1 <- as.factor(bootBase$V1)
bootEx$V1 <- as.factor(bootEx$V1)

quantileCheck <- function(x){
  tmpArray <- array(0,c(nrow(x),951))
  for (r in 1:nrow(x)) {
    tmp <- t(x[r,2:1001])
    tmp <- tmp[order(tmp)]
    tmp <- tmp[26:975]
    tmpArray[r,2:951] <- tmp
  }
  tmpArray <- data.frame(tmpArray)
  tmpArray[,1] <- x[,1]
  return(tmpArray)
}
bootMI <- quantileCheck(bootMI)
bootBase <- quantileCheck(bootBase)
bootEx <- quantileCheck(bootEx)

bootMI <- melt(bootMI, 'X1', variable.name = 'bootNum', value.name = 'accValue')
bootBase <- melt(bootBase, 'X1', variable.name = 'bootNum', value.name = 'accValue')
bootEx <- melt(bootEx, 'X1', variable.name = 'bootNum', value.name = 'accValue')


ebBottom <- function(x){
  return(quantile(x, 0.025))
}

ebTop <- function(x){
  return(quantile(x, 0.975))
}
myLevels <- c('calDataMSMI', 'calDataAreaMI', 'calDataThicknessMI', 
              'calDataVolumeMI', 'calDataMCMI', 'calDataGCMI')
bootMI$X1 <- factor(bootMI$X1, levels = myLevels)
myCompar <- list(c('calDataMSMI', 'calDataAreaMI'), c('calDataMSMI', 'calDataThicknessMI'),
                 c('calDataMSMI', 'calDataVolumeMI'), c('calDataMSMI', 'calDataMCMI'), c('calDataMSMI', 'calDataGCMI'))
ggplot(bootMI, aes(X1, accValue))+
  stat_summary(geom = 'bar', fun = 'mean', aes(fill = X1))+
  geom_jitter(alpha = 0.3, size = 0.8)+
  stat_summary(geom = 'errorbar', 
               fun.max = ebTop,
               fun.min = ebBottom,
               width = 0.2, size = 1.2)+
  stat_compare_means(comparisons = myCompar, method = 't.test', bracket.size = 1, size = 5)+
  annotate('point', x=0.8, y=0.635,size=5)+
  annotate('point', x=1.8, y=0.620,size=5)+
  annotate('point', x=2.8, y=0.619,size=5)+
  annotate('point', x=3.8, y=0.616,size=5)+
  annotate('point', x=4.8, y=0.605,size=5)+
  annotate('point', x=5.8, y=0.599,size=5)+
  theme(legend.position = 'none',
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
ggsave('D:/毕业论文temp/resultsPlt/predTtestMI.pdf')
myCompar2 <- list(c('calDataMSBase', 'calDataAreaBase'), c('calDataMSBase', 'calDataThicknessBase'),
                 c('calDataMSBase', 'calDataVolumeBase'), c('calDataMSBase', 'calDataMCBase'), c('calDataMSBase', 'calDataGCBase'))
bootBase$X1 <- factor(bootBase$X1, levels = c('calDataMSBase', 'calDataAreaBase', 'calDataThicknessBase', 
                                              'calDataVolumeBase', 'calDataMCBase', 'calDataGCBase'))
ggplot(bootBase, aes(X1, accValue))+
  stat_summary(geom = 'bar', fun = 'mean', aes(fill = X1))+
  geom_jitter(alpha = 0.3, size = 0.8)+
  stat_summary(geom = 'errorbar', 
               fun.max = ebTop,
               fun.min = ebBottom,
               width = 0.2, size = 1.2)+
  stat_compare_means(comparisons = myCompar2, method = 't.test', bracket.size = 1, size = 5)+
  annotate('point', x=0.8, y=0.670,size=5)+
  annotate('point', x=1.8, y=0.689,size=5)+
  annotate('point', x=2.8, y=0.701,size=5)+
  annotate('point', x=3.8, y=0.715,size=5)+
  annotate('point', x=4.8, y=0.662,size=5)+
  annotate('point', x=5.8, y=0.670,size=5)+
  theme(legend.position = 'none',
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
ggsave('D:/毕业论文temp/resultsPlt/predTtestBase.pdf')
myCompar3 <- list(c('subDataExMS', 'subDataExArea'), c('subDataExMS', 'subDataExThickness'),
                  c('subDataExMS', 'subDataExVolume'), c('subDataExMS', 'subDataExMC'), c('subDataExMS', 'subDataExGC'))
bootEx$X1 <- factor(bootEx$X1, levels = c('subDataExMS', 'subDataExArea', 'subDataExThickness', 
                                          'subDataExVolume', 'subDataExMC', 'subDataExGC'))
ggplot(bootEx, aes(X1, accValue))+
  stat_summary(geom = 'bar', fun = 'mean', aes(fill = X1))+
  geom_jitter(alpha = 0.3, size = 0.8)+
  stat_summary(geom = 'errorbar', 
               fun.max = ebTop,
               fun.min = ebBottom,
               width = 0.2, size = 1.2)+
  stat_compare_means(comparisons = myCompar3, method = 't.test', bracket.size = 1, size = 5)+
  annotate('point', x=0.8, y=0.701,size=5)+
  annotate('point', x=1.8, y=0.705,size=5)+
  annotate('point', x=2.8, y=0.710,size=5)+
  annotate('point', x=3.8, y=0.725,size=5)+
  annotate('point', x=4.8, y=0.687,size=5)+
  annotate('point', x=5.8, y=0.667,size=5)+
  theme(legend.position = 'none',
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
ggsave('D:/毕业论文temp/resultsPlt/predTtestEx.pdf')  

## seq3 and all time cortest
# Sex.tvalue or time.tvalue
allArea <- resSex4calDataArea$Sex.tvalue
allThick <- resSex4calDataThickness$Sex.tvalue
allVolume <- resSex4calDataVolume$Sex.tvalue
allMS <- resSex4calDataMS$Sex.tvalue
allMC <- resSex4calDataMeanCurve$Sex.tvalue
allGC <- resSex4calDataGausCurve$Sex.tvalue
permId <- read.csv('D:/毕业论文temp/PLS/permId1000.csv')
measuresPerm <- c('Area', 'Thickness', 'Volume', 'MS', 'MeanCurve', 'GausCurve')
allPerm <- c('allArea', 'allThick', 'allVolume', 'allMS', 'allMC', 'allGC')
no_overlap <- length(seq)-1
overlap <- length(seqWindow)
slide <- nrow(seqSlide)
corrStore <- data.frame(seq=rep(paste0('seq', rep(1:slide)),each = 6),
                        measures=rep(c('Area', 'Thickness', 'Volume', 'MS', 'meanCurv', 'gausCurv'), slide),
                        corr=rep(0, slide*6), pSpin=rep(0, slide*6))  ## (length(seq)-1) or length(seqWindow)
for (m in 1:6) {
  dfTmp <- measuresPerm[m]
  allTmp <- get(allPerm[m])
  for (s in 1:slide) {
    cmd <- paste0('seq', s, ' <- resSex',s , dfTmp,'$Sex.tvalue') ## resSex or windresSex
    eval(parse(text = cmd))
  }
  # seq1 <- get(paste0('resSex1', dfTmp))$time.tvalue
  # seq2 <- get(paste0('resSex2', dfTmp))$time.tvalue
  # seq3 <- get(paste0('resSex3', dfTmp))$time.tvalue
  corPermTest <- array(0, c(1000, slide))
  for (p in 1:ncol(permId)) {
    idIdx <- permId[,p]
    allP <- allTmp[idIdx]
    for (i in 1:ncol(corPermTest)) {
      corPermTest[p, i] <- cor(get(paste0('seq', i)), allP, method = 'spearman')
    }
    # corPermTest[p, 1] <- cor(seq1, allP, method = 'spearman')
    # corPermTest[p, 2] <- cor(seq2, allP, method = 'spearman')
    # corPermTest[p, 3] <- cor(seq3, allP, method = 'spearman')
    
  }
  for (p in 1:ncol(corPermTest)) {
    p1 <- length(corPermTest[corPermTest[,1]>=cor(seq1, allTmp, method = 'spearman'), 1])/1000
    cmd <- paste0('p', p," <- length(corPermTest[corPermTest[,",p,"]>=cor(seq", p, ", allTmp, method = 'spearman'),", p, "])/1000")
    eval(parse(text = cmd))
    
  }
  # p1 <- length(corPermTest[corPermTest[,1]>=cor(seq1, allTmp, method = 'spearman'), 1])/1000
  # p2 <- length(corPermTest[corPermTest[,2]>=cor(seq2, allTmp, method = 'spearman'), 2])/1000
  # p3 <- length(corPermTest[corPermTest[,3]>=cor(seq3, allTmp, method = 'spearman'), 3])/1000
  # print(paste0(dfTmp, '_corr: ', cor(seq1, allTmp, method = 'spearman'), '  pSpin: ', p1))
  # print(paste0(dfTmp, '_corr: ', cor(seq2, allTmp, method = 'spearman'), '  pSpin: ', p2))
  # print(paste0(dfTmp, '_corr: ', cor(seq3, allTmp, method = 'spearman'), '  pSpin: ', p3))
  
  for (f in 1:slide) {
    corrStore[m+(f-1)*6, 3] <- cor(get(paste0('seq',f)), allTmp, method = 'spearman')
    corrStore[m+(f-1)*6, 4] <- get(paste0('p', f))
  }
  
  # corrStore[m, 3] <- cor(seq1, allTmp, method = 'spearman')
  # corrStore[m+6, 3] <- cor(seq2, allTmp, method = 'spearman')
  # corrStore[m+12, 3] <- cor(seq3, allTmp, method = 'spearman')
  # corrStore[m, 4] <- p1
  # corrStore[m+6, 4] <- p2
  # corrStore[m+12, 4] <- p3
}
cor.test(windresSex2Area$Sex.tvalue, resSex4calDataArea$Sex.tvalue, method = 'spearman')

corrStore$seq <- factor(corrStore$seq, levels = c((paste0('seq', rep(1:slide)))))
plt2 <- ggplot(corrStore, aes(seq, corr))+
  geom_point(aes(color = measures))+
  geom_line(aes(group = measures, color = measures))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust = 0.5))+
  ylab("correlation with all seq")
ggsave('D:/毕业论文temp/resultsPlt/seqCorr_noOverlaptimebysex.pdf',plt2,height = 3, width = 4)
hist(corPermTestArea[,1])
mean(corPermTestArea[,3])
cor.test(resSex4calDataMS$Sex.tvalue, resSex3Area$Sex.tvalue)
cor.test(resSex4calDataMS$Sex.tvalue, resSex3Thickness$Sex.tvalue)
cor.test(resSex4calDataMS$Sex.tvalue, resSex3Volume$Sex.tvalue)

cor.test(resSex3MS$time.tvalue, resSex3Volume$time.tvalue)
cor.test(resSex2MS$time.tvalue, resSex2Volume$time.tvalue)
cor.test(resSex1MS$time.tvalue, resSex1Volume$time.tvalue)

## map MS result to network level
regionAssign <- read.table('D:/毕业论文temp/resultsPlt/cortex_parcel_network_assignments.txt')
tmp <- resSex4calDataMS
Network <- c("Visual1","Visual2","Somatomotor","Cingulo-Opercular","Dorsal-Attention",
             "Language","Frontoparietal","Auditory","Default","Posterior-Multimodal",
             "Ventral-Multimodal","Orbito-Affective")
myCols <- c("#f93bfb","#981497","#FF0000","#01fc00","#fffb00","#069a98",
            "#377b00","#b05a24","#00fdfd","#fa9e00","#0037fd","#5a35fd")
tmp$network <- factor(regionAssign[,1], levels = c(rep(1:12)), labels = Network)
tmp$col <- factor(regionAssign[,1], levels = c(rep(1:12)), labels = myCols)
tmp$network <- factor(tmp$network, levels = c('Visual1', 'Posterior-Multimodal', 'Frontoparietal', 'Visual2',
                                              'Orbito-Affective', 'Dorsal-Attention', 'Default', 'Ventral-Multimodal',
                                              'Somatomotor', 'Auditory', 'Language', 'Cingulo-Opercular'))
voiPlt <- ggplot(tmp, aes(network, Sex.tvalue))+
  geom_violin(aes(fill = network), scale = 'width')+
  geom_boxplot(width = 0.1)+
  theme_classic() +
  scale_fill_manual(values = myCols)+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
ggsave('D:/毕业论文temp/resultsPlt/voiPltMycol.pdf', width = 9.35)
cor.test(resSex4calDataMS$Sex.tvalue, resSex4calDataMS$time.tvalue)
time_timebysexPlt <- ggplot(tmp, aes(Sex.tvalue, time.tvalue))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = 'lm', alpha = 0.5)
ggsave('D:/毕业论文temp/resultsPlt/time_timebysexPlt.pdf', time_timebysexPlt, width = 4, height = 5)


## group pattern MS
check <- glasser
label <- str_remove(colnames(calDataMS[50:409]), '_ROI_.*')
label <- str_replace(label, '\\.', '-')
tmpCalMS <- subset(calDataMS, calDataMS$label=='Time1'&calDataMS$site==1)
pattern <- colMeans(calDataMS[,50:409])
patternPlt <- data.frame(ms=pattern, label=label)
msPattern <- ggseg(patternPlt,atlas = glasser, 
       mapping = aes(fill = ms),colour = "black",
       show.legend = T, 
       position = "stacked") +
  scale_fill_gradient2(midpoint = 0,mid = "green",low = "blue",high = "red")+
  labs(title = 'MS pattern', fill = "ms-value")+
  theme_void()+
  theme(text=element_text(size=16,  family="serif"))
ggsave('D:/毕业论文temp/resultsPlt/msPattern.pdf', msPattern, width = 4, height = 5)

## baseline and intercept pattern

permId <- read.csv('D:/毕业论文temp/PLS/permId1000.csv')
corrPatternArea <- getBaseSlopeCorr(calDataArea, permId)
corrPatternThickness <- getBaseSlopeCorr(calDataThickness, permId)
corrPatternVolume <- getBaseSlopeCorr(calDataVolume, permId)
corrPatternMC <- getBaseSlopeCorr(calDataMeanCurve, permId)
corrPatternGC <- getBaseSlopeCorr(calDataGausCurve, permId)
corrPatternMS <- getBaseSlopeCorr(calDataMS, permId)

corrPatternAll <- c('corrPatternArea', 'corrPatternThickness', 'corrPatternVolume', 
                    'corrPatternMC', 'corrPatternGC', 'corrPatternMS')
if(1){
  for (p in 1:length(corrPatternAll)) {
    tmpName <- corrPatternAll[p]
    tmp <- get(tmpName)
    ggseg (tmp, atlas = glasser, 
           mapping = aes(fill = corr), 
           show.legend = T, 
           position = "stacked") +
      scale_fill_gradient2(midpoint = 0, mid = "white", low = "blue", high = "red")+
      theme_void()+
      theme(text=element_text(size=12,  family="serif"))
    ggsave(paste0('D:/毕业论文temp/resultsPlt/base_long_corr/', tmpName, '.pdf'))
  }
  
  for (p in 1:length(corrPatternAll)) {
    tmpName <- corrPatternAll[p]
    tmp <- get(tmpName)
    tmp$p_adj <- p.adjust(tmp$pvalue, method = 'fdr')
    tmp <- tmp[tmp$p_adj<=0.05,]
    ggseg (tmp, atlas = glasser, 
           mapping = aes(fill = corr), 
           show.legend = T, 
           position = "stacked") +
      scale_fill_gradient2(midpoint = 0, mid = "white", low = "blue", high = "red")+
      theme_void()+
      theme(text=element_text(size=12,  family="serif"))
    ggsave(paste0('D:/毕业论文temp/resultsPlt/base_long_corr/', tmpName, '_adj.pdf'))
  }
  
  for (p in 1:length(corrPatternAll)) {
    tmpName <- corrPatternAll[p]
    tmp <- get(tmpName)
    ggseg (tmp, atlas = glasser, 
           mapping = aes(fill = base), 
           show.legend = T, 
           position = "stacked") +
      scale_fill_gradient2(midpoint = 0, mid = "white", low = "blue", high = "red")+
      theme_void()+
      theme(text=element_text(size=12,  family="serif"))
    ggsave(paste0('D:/毕业论文temp/resultsPlt/base_long_corr/', tmpName, '_base.pdf'))
  }
  
  for (p in 1:length(corrPatternAll)) {
    tmpName <- corrPatternAll[p]
    tmp <- get(tmpName)
    ggseg (tmp, atlas = glasser, 
           mapping = aes(fill = slope), 
           show.legend = T, 
           position = "stacked") +
      scale_fill_gradient2(midpoint = 0, mid = "white", low = "blue", high = "red")+
      theme_void()+
      theme(text=element_text(size=12,  family="serif"))
    ggsave(paste0('D:/毕业论文temp/resultsPlt/base_long_corr/', tmpName, '_slope.pdf'))
  }
}

