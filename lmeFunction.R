lmePipeline <- function(df, outputPath, times){
  # df <- reDf1
  # df <- mcData
  # library packages
  if(require(tidyverse)==0){
    install.packages("tidyverse")
    library(tidyverse)
  }
  if(require(lmerTest)==0){
    install.packages("lmerTest")
    library(lmerTest)
  }

  #df <- group_20Vol
  

  df[,-c(1:3,418:427)] <- scale(df[,-c(1:3,418:427)])
  ModelData <- df
  
  Brain <- colnames(ModelData[,50:409])
  Var <- c("Sex","Beta","E","O","Com","Agen","Ope","Int","time") #????alphaָplasticity??betaָstability
  ResultDf <- matrix(0,length(Brain),3*length(Var))
  colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"),paste(Var,"df",sep = "-"))
  rownames(ResultDf) <- colnames(ModelData[,50:409])
  n_sites <- length(unique(ModelData$site))
  for (vv in 1:length(Var)){
    lmstore <- list()
    # centertime
    Inter <- paste0("centertime*",Var)
    Coef <- paste0("centertime:",Var)
    if (vv == 9) Coef[vv] <- "centertime"
    print(Coef[vv])
    
    if(n_sites>=2){
      for (aa in 1:360) {
        print(paste(Var[vv],aa))
         ll <- switch(vv,
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Beta*centertime+site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Alpha*centertime+site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+O*centertime+C*centertime+A*centertime+N*centertime+site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+C*centertime+E*centertime+A*centertime+N*centertime+site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Agen*centertime+O*centertime+C*centertime+A*centertime+N*centertime++site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Com*centertime+O*centertime+C*centertime+A*centertime+N*centertime+site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Int*centertime+site+C*centertime+E*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Ope*centertime+site+C*centertime+E*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean*centertime+site+(1|ID)"))
        )
        # lmstore[[aa]] <- lmer(ll,data = ModelData, control = lmerControl(optimizer ="Nelder_Mead"))
        lmstore[[aa]] <- lmer(ll,data = ModelData)
      }
    }else {
      for (aa in 1:360) {
        print(paste(Var[vv],aa))
        ll <- switch(vv,
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Beta*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Alpha*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+O*centertime+C*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+C*centertime+E*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Agen*centertime+O*centertime+C*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Com*centertime+O*centertime+C*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Int*centertime+C*centertime+E*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+Ope*centertime+C*centertime+E*centertime+A*centertime+N*centertime+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean*centertime+(1|ID)"))
        )
        lmstore[[aa]] <- lmer(ll,data = ModelData)
      }
    }
    
    if(vv==1) Coef[vv] = "Sex:centertime"
    PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef[vv],]))

    ResultDf[,vv] <- PVector[5,]
    ResultDf[,vv+length(Var)] <- PVector[4,]
    ResultDf[,vv+(2*length(Var))] <- PVector[3,]
  }
  
  ResultDf <- data.frame(ResultDf)
  PlotLabel <- str_replace(rownames(ResultDf),"_ROI.*","")
  ResultDf$label <- PlotLabel
  ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
  
  dataType <- unique(str_remove(rownames(ResultDf), ".*ROI_"))
  if(dir.exists(paste0(outputPath, dataType))==F) dir.create(paste0(outputPath, dataType))
  write.csv(ResultDf,paste0(outputPath, dataType, "/", dataType, times, "ResUncorrected.csv"),row.names = F)
  # FDR-correct
  StatsName <- paste0(dataType,Var,"Stats")
  StatsNameF <- paste0(dataType,Var,"StatsFDR")
  returnDf <- ResultDf
  for (vv in 1:length(Var)) {
    cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+length(Var),"],label=ResultDf$label)")
    cmd2 <- paste0("returnDf[,", vv,"] <- ", StatsName[vv],"[,1]")
    tvalueCorr <- paste0(colnames(returnDf)[vv+9],".Corr")
    cmd3 <- paste0("returnDf$", tvalueCorr," <- 1")
    eval(parse(text = cmd))
    eval(parse(text = cmd2))
    eval(parse(text = cmd3))
    for (i in 1:nrow(ResultDf)) {
      if (returnDf[i,vv+9]<=0) returnDf[i,tvalueCorr] <- -qt(returnDf[i,vv]/2, returnDf[i,vv+18],lower=FALSE)
      if (returnDf[i,vv+9]>0) returnDf[i,tvalueCorr] <- qt(returnDf[i,vv]/2, returnDf[i,vv+18],lower=FALSE)
    }
  }
  # Extract sig regions
  for (ss in 1:length(Var)) {
    if (nrow(get(StatsName[ss]))>=1) {
      cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
      eval(parse(text = cmd))
    }
  }
  # Can not survival after FDR-correct 
  temp <- data.frame(BehaviorDim=StatsNameF,FDR=1)
  for (ss in 1:length(Var)) {
    if(nrow(get(StatsNameF[ss]))==0) {
      # Print regions can not pass FDR-correct
      temp[ss,2] <- 0
      # Threshold the original p-value at 0.001
      cmd1 <- paste0(StatsNameF[ss]," <- subset(ResultDf",",ResultDf[",ss,"]<0.05)")
      # Select p&t values and labels
      cmd2 <- paste0(StatsNameF[ss]," <- ",StatsNameF[ss],"[,c(",ss,",",ss,"+length(Var),",2*length(Var)+1,")]")
      # Rename colnames
      cmd3 <- paste0("colnames(",StatsNameF[ss],") <- ","colnames(",StatsName[ss],")")
      eval(parse(text = cmd1))
      eval(parse(text = cmd2))
      eval(parse(text = cmd3))
    }
  }
  write.table(temp,paste0(outputPath, dataType, "/", dataType, times, "FDR.txt"),quote = F,row.names = F)

  return(returnDf)
}

lmePipelineSex <- function(df, outputPath, times){
  # df <- calDataArea
  # df <- mcData
  # library packages
  # outputPath <- outputPathVal
  if(require(tidyverse)==0){
    install.packages("tidyverse")
    library(tidyverse)
  }
  if(require(lmerTest)==0){
    install.packages("lmerTest")
    library(lmerTest)
  }
  
  #df <- group_20Vol
  
  df[,-c(1:3,417)] <- scale(df[,-c(1:3,417)])
  ModelData <- df
  
  Brain <- colnames(ModelData[,50:409])
  Var <- c("Sex", "time") #????alphaָplasticity??betaָstability
  ResultDf <- matrix(0,length(Brain),3*length(Var))
  AICBIC <- matrix(0, 360, 2*length(Var)+1)
  
  colnames(AICBIC) <- c(paste0(rep(Var, length(Var)), rep(c("AIC", "BIC"), each = length(Var))), "pval")
  colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"),paste(Var,"df",sep = "-"))
  rownames(ResultDf) <- colnames(ModelData[,50:409])
  n_sites <- length(unique(ModelData$site))
  for (vv in 1:length(Var)){
    lmstore <- list()
    # centertime
    Inter <- paste0("centertime*",Var)
    Coef <- paste0("centertime:",Var)
    if (vv == 2) Coef[vv] <- "centertime"
    print(Coef[vv])
    
    if(n_sites>=2){
      for (aa in 1:360) {
        print(paste(Var[vv],aa))
        ll <- switch(vv,
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+site+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+centertime+site+(1|ID)"))
        )
        # lmstore[[aa]] <- lmer(ll,data = ModelData, control = lmerControl(optimizer ="Nelder_Mead"))
        lmstore[[aa]] <- lmer(ll,data = ModelData)

      }
      cmd <- paste0("lmstore", Var[vv], " <- lmstore")
      eval(parse(text = cmd))
      
    }else {
      for (aa in 1:360) {
        print(paste(Var[vv],aa))
        ll <- switch(vv,
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean+",Inter[vv],"+(1|ID)")),
                     as.formula(paste0(Brain[aa], "~ Sex+age_mean*centertime+(1|ID)"))
        )
        lmstore[[aa]] <- lmer(ll,data = ModelData)
      }
    }

    if(vv==1) Coef[vv] = "Sex:centertime"
    PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef[vv],]))
    
    ResultDf[,vv] <- PVector[5,]
    ResultDf[,vv+length(Var)] <- PVector[4,]
    ResultDf[,vv+(2*length(Var))] <- PVector[3,]
  }
  # for (a in 1:360) {
  #   print(paste0("anova: ", a))
  #   # a <- 15
  #   # a <- 1
  #   ANO <- anova(lmstoreSex[[a]], lmstoretime[[a]])
  #   AICBIC[a, 1] <- ANO["lmstoreSex[[a]]", "AIC"]
  #   AICBIC[a, 2] <- ANO["lmstoretime[[a]]", "AIC"]
  #   AICBIC[a, 3] <- ANO["lmstoreSex[[a]]", "BIC"]
  #   AICBIC[a, 4] <- ANO["lmstoretime[[a]]", "BIC"]
  #   AICBIC[a, 5] <- ANO["lmstoreSex[[a]]", "Pr(>Chisq)"]
  # }
  ResultDf <- data.frame(ResultDf)
  PlotLabel <- str_replace(rownames(ResultDf),"_ROI.*","")
  ResultDf$label <- PlotLabel
  ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
  
  dataType <- unique(str_remove(rownames(ResultDf), ".*ROI_"))
  file.path(outputPath)
  
  # dir create
  outdir <- unlist(strsplit(outputPath, "/"))
  create <- outdir[1]
  for (d in 2:length(outdir)) {
    create <- paste0(create, "/", outdir[d])
    if(!dir.exists(create)) dir.create(create)
  }
  if(!dir.exists(paste0(create, "/Sex"))) dir.create(paste0(create, "/Sex"))
  if(!dir.exists(paste0(outputPath, "Sex/", dataType))) dir.create(paste0(outputPath, "Sex/", dataType))
  write.csv(ResultDf,paste0(outputPath, "Sex/", dataType, "/", dataType, times, "ResUncorrected.csv"),row.names = F)
  # FDR-correct
  StatsName <- paste0(dataType,Var,"Stats")
  StatsNameF <- paste0(dataType,Var,"StatsFDR")
  returnDf <- ResultDf
  for (vv in 1:length(Var)) {
    cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+length(Var),"],label=ResultDf$label)")
    cmd2 <- paste0("returnDf[,", vv,"] <- ", StatsName[vv],"[,1]")
    tvalueCorr <- paste0(colnames(returnDf)[vv+2],".Corr")
    cmd3 <- paste0("returnDf$", tvalueCorr," <- 1")
    eval(parse(text = cmd))
    eval(parse(text = cmd2))
    eval(parse(text = cmd3))
    for (i in 1:nrow(ResultDf)) {
      if (returnDf[i,vv+2]<=0) returnDf[i,tvalueCorr] <- -qt(returnDf[i,vv]/2, returnDf[i,vv+4],lower=FALSE)
      if (returnDf[i,vv+2]>0) returnDf[i,tvalueCorr] <- qt(returnDf[i,vv]/2, returnDf[i,vv+4],lower=FALSE)
    }
  }
  # Extract sig regions
  for (ss in 1:length(Var)) {
    if (nrow(get(StatsName[ss]))>=1) {
      cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
      eval(parse(text = cmd))
    }
  }
  # Can not survival after FDR-correct 
  temp <- data.frame(BehaviorDim=StatsNameF,FDR=1)
  for (ss in 1:length(Var)) {
    if(nrow(get(StatsNameF[ss]))==0) {
      # Print regions can not pass FDR-correct
      temp[ss,2] <- 0
      # Threshold the original p-value at 0.001
      cmd1 <- paste0(StatsNameF[ss]," <- subset(ResultDf",",ResultDf[",ss,"]<0.05)")
      # Select p&t values and labels
      cmd2 <- paste0(StatsNameF[ss]," <- ",StatsNameF[ss],"[,c(",ss,",",ss,"+length(Var),",2*length(Var)+1,")]")
      # Rename colnames
      cmd3 <- paste0("colnames(",StatsNameF[ss],") <- ","colnames(",StatsName[ss],")")
      eval(parse(text = cmd1))
      eval(parse(text = cmd2))
      eval(parse(text = cmd3))
    }
  }
  write.table(temp,paste0(outputPath, "Sex/",dataType, "/", dataType, times, "FDR.txt"),quote = F,row.names = F)
  
  returnDf <- cbind(returnDf, AICBIC)
  return(returnDf)
}


# demPlot
demPlot <- function(df){
  if(require(tidyverse)==0){
    install.packages("tidyverse")
    library(tidyverse)
  }
  Voldata <- df
  #Voldata <- ModelData
  subID <- Voldata[!duplicated(Voldata$ID),c(1,4)]
  idReOder <- subID[order(subID$Age1),]
  idReOder$rank <- 1:nrow(idReOder)
  colnames(idReOder)[2] <- "Age1R"
  plotData <- inner_join(Voldata,idReOder, by="ID")
  mycols<-c("blue","purple","red","orange")
  demography <- ggplot(plotData,aes(Age1,rank,group=ID))+
    geom_point(aes(colour=label),size = 2.2,alpha = 0.8)+
    geom_line(aes(colour=label), alpha = 0.5)+
    labs(x="Age of Scanning(years)",y = "Subject",colour = "Scanning")+
    scale_colour_manual(values = mycols)+
    theme(axis.text = element_text(size = 12),text = element_text(size = 12),
          axis.line = element_line(color = "black",size = 1),axis.ticks = element_line(size = 1),
          panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())#+
  geom_vline(xintercept = 12,color="red",size=2)
  return(demography)
}

# split data
dataSplit <- function(subGbb, subSlim, percentage){
  # calculate sub number
  nSubGbb <- length(subGbb)
  nSubSlim <- length(subSlim)
  # random
  randIndexGbb <- sample(1:nSubGbb, replace = F)
  randIndexSlim <- sample(1:nSubSlim, replace = F)
  # extract index
  n_extSubGbb <- round(nSubGbb * percentage)
  n_extSubSlim <- round(nSubSlim * percentage)
  
  extIndexGbb <- randIndexGbb[1:n_extSubGbb]
  extIndexSlim <- randIndexSlim[1:n_extSubSlim]
  
  testSetIndexGbb <- subGbb[extIndexGbb]
  testSetIndexSlim <- subSlim[extIndexSlim]
  
  testSetIndex <- list(testSetIndexGbb = testSetIndexGbb, testSetIndexSlim = testSetIndexSlim,
                       nSubGbb = nSubGbb, nSubSlim = nSubSlim,
                       n_extSubGbb = n_extSubGbb, n_extSubSlim = n_extSubSlim)
  return(testSetIndex)
}

# apply split
applySplitSingle <- function(spliteIndexSlim, df){
  indexTest <- spliteIndexSlim
  #?ϲ?????ģ̬֮????df
  trainSet <- df
  for (i in 1:length(indexTest)) {
    trainSet <- subset(trainSet, trainSet$ID != indexTest[i])
  }
  
  indexAll <- unique(df$ID)
  indexTrain <- indexAll
  for (i in 1:length(indexTest)) {
    indexTrain <- subset(indexTrain,indexTrain!=indexTest[i])
  }
  
  testSet <- df
  for (i in 1:length(indexTrain)) {
    testSet <- subset(testSet, testSet$ID != indexTrain[i])
  }
  train_test_split <- list(trainSet = trainSet, testSet = testSet)
  return(train_test_split)
}

applySplit <- function(spliteIndexSlim, spliteIndexGbb, df){
  indexTest <- c(spliteIndexSlim, spliteIndexGbb)
  #?ϲ?????ģ̬֮????df
  trainSet <- df
  for (i in 1:length(indexTest)) {
    trainSet <- subset(trainSet, trainSet$ID != indexTest[i])
  }
  
  indexAll <- unique(df$ID)
  indexTrain <- indexAll
  for (i in 1:length(indexTest)) {
    indexTrain <- subset(indexTrain,indexTrain!=indexTest[i])
  }
  
  testSet <- df
  for (i in 1:length(indexTrain)) {
    testSet <- subset(testSet, testSet$ID != indexTrain[i])
  }
  train_test_split <- list(trainSet = trainSet, testSet = testSet)
  return(train_test_split)
}
# calculate MS
calMs <- function(ctData, gcData, gmData, mcData, saData){
  regIndex <- grep(colnames(gcData),pattern = "ROI")
  startIndex <- regIndex[1]
  endIndex <- regIndex[length(regIndex)]
  
  n <- ncol(ctData[,startIndex : endIndex])
  if (n != 360) print("region numbers != 360")
  
  ctZscore <- data.frame(scale(t(ctData[,startIndex: endIndex])))
  gcZscore <- data.frame(scale(t(gcData[,startIndex: endIndex])))
  gmZscore <- data.frame(scale(t(gmData[,startIndex: endIndex])))
  mcZscore <- data.frame(scale(t(mcData[,startIndex: endIndex])))
  saZscore <- data.frame(scale(t(saData[,startIndex: endIndex])))
  
  sub <- array(0, dim = c(5, nrow(ctZscore), ncol(ctZscore)))
  subMsMtx <- list()
  regionalMs <- data.frame(array(0, dim = c(ncol(ctZscore), nrow(ctZscore))))
  colnames(regionalMs) <- gsub("_ROI.*", "", colnames(ctData)[startIndex: endIndex])
  regionalMs$ID <- ctData$ID
  regionalMs$label <- ctData$label
  for (i in 1:dim(sub)[3]) {
    print(paste("sub", i, ctData$ID[i], sep = "-"))
    sub[1,,i] <- ctZscore[,i]
    sub[2,,i] <- gcZscore[,i]
    sub[3,,i] <- gmZscore[,i]
    sub[4,,i] <- mcZscore[,i]
    sub[5,,i] <- saZscore[,i]
    
    c1 <- cor(sub[,,i])
    subMsMtx[[i=paste0(ctData$ID[i],ctData$label[i])]] <- c1
    
    diag(c1) <- 0
    ms <- rep(0,360)
    for (j in 1:nrow(c1)) {
      ms[j] <- sum(c1[,j])/(nrow(c1)-1)
    }
    regionalMs[i,1:360] <- ms
  }
  reList <- list(regionalMs = regionalMs, subMsMtx = subMsMtx)
  return(reList)
}

calLongAndCrossCov <- function(x){
  # x <- reDf1Area
  # x <- reDfWindow.reDf7Subc
  # x <- reDfWindow.reDf7MS
  
  ## 我不该把这部分写在这里，但是我懒得改了所以我决定结尾重新算一次 ##
  ageMeanSub <- aggregate(x$Age1, by = list(x$ID), mean)
  colnames(ageMeanSub) <- c("ID", "ageMeanSub")
  indexCol <- rep(1:ncol(x))
  rmIndex <- c(indexCol[colnames(x)=="AgeMeanSub"], indexCol[colnames(x)=="age_mean"])
  x <- x[, -rmIndex]
  temp <- inner_join(x, ageMeanSub)
  
  ageMeanSubAll <- mean(x$Age1)
  temp$age_mean <- temp$ageMeanSub-ageMeanSubAll
  temp$centertime <- temp$Age1 - temp$ageMeanSub
  temp$centertime2 <- (temp$Age1 - temp$ageMeanSub)^2
  #######################################################################
  
  subNum <- unique(temp$ID)
  rmmmmm <- rep(0, nrow(temp))
  rmmmmm2 <- rep(0, nrow(temp))
  for (s in 1:length(subNum)) {
    # s <- 101
    # s <- 138
    tempSub <- temp[temp$ID == subNum[s],]
    rowNum <- nrow(tempSub)
    
    
    if (rowNum > 2) {
      if (rowNum == 3) {
        tempSeq <- rep(0, 3)
        tempSeq[2] <- tempSub[2, "Age1"] - tempSub[1, "Age1"]
        tempSeq[3] <- tempSub[3, "Age1"] - tempSub[2, "Age1"]
        
        tempIndex <- c(1:3)
        rmIndex <- tempIndex[tempSeq == max(tempSeq)]
        
        rmSub <- tempSub[rmIndex, c("ID", "label")]
        tempRowNum <- 1:nrow(temp)
        rmmmmm[s] <- tempRowNum[(temp[, "ID"]==rmSub$ID&temp[, "label"]==rmSub$label)]
      } else if(rowNum == 4) {
        tempSeq <- rep(0, 4)
        tempSeq[2] <- tempSub[2, "Age1"] - tempSub[1, "Age1"]
        tempSeq[3] <- tempSub[3, "Age1"] - tempSub[2, "Age1"]
        tempSeq[4] <- tempSub[4, "Age1"] - tempSub[3, "Age1"]
        
        tempIndex <- c(1:4)
        saveIndex1 <- tempIndex[tempSeq == sort(tempSeq)[2]]
        saveIndex2 <- tempIndex[tempSeq == sort(tempSeq)[2]]-1
        rmIndex <- c(5-saveIndex1, 5-saveIndex2)
        
        rmSub <- tempSub[rmIndex, c("ID", "label")]
        tempRowNum <- 1:nrow(temp)
        rmmmmm2[s] <- tempRowNum[(temp[, "ID"]==rmSub[1,"ID"]&temp[, "label"]==rmSub[1,"label"])]
        rmmmmm2[s+1] <- tempRowNum[(temp[, "ID"]==rmSub[2,"ID"]&temp[, "label"]==rmSub[2,"label"])]
      }
    }
  }
  rmmmmm <- rmmmmm[rmmmmm!=0]
  rmmmmm2 <- rmmmmm2[rmmmmm2!=0]
  rmFinal <- c(rmmmmm, rmmmmm2)
  
  if (length(rmFinal)!=0) {
    temp <- temp[-rmFinal, ]
  }
  reOrderLabel <- unique(temp$ID)
  for (r in 1:length(reOrderLabel)) {
    rowNum <- length(temp[temp$ID == subNum[r],"label"])
    temp[temp$ID == subNum[r],"label"] <- paste0(rep("Time", rowNum), rep(1 : rowNum))
  }
  
  ageMeanSub <- aggregate(temp$Age1, by = list(temp$ID), mean)
  colnames(ageMeanSub) <- c("ID", "ageMeanSub")
  indexCol <- rep(1:ncol(temp))
  rmIndex <- c(indexCol[colnames(temp)=="AgeMeanSub"], indexCol[colnames(temp)=="age_mean"])
  temp2 <- inner_join(temp, ageMeanSub)
  
  ageMeanSubAll <- mean(temp$Age1)
  temp$age_mean <- temp$ageMeanSub-ageMeanSubAll
  temp$centertime <- temp$Age1 - temp$ageMeanSub
  temp$centertime2 <- (temp$Age1 - temp$ageMeanSub)^2
  return(temp2)
}

calLongAndCrossCov2 <- function(x){
  # x <- reDf3MeanCurve
  # x <- reDfWindow.reDf7Subc
  # x <- reDfWindow.reDf7MS
  
  ageMeanSub <- aggregate(x$Age1, by = list(x$ID), mean)
  colnames(ageMeanSub) <- c("ID", "ageMeanSub")
  indexCol <- rep(1:ncol(x))
  rmIndex <- c(indexCol[colnames(x)=="AgeMeanSub"], indexCol[colnames(x)=="age_mean"])
  x <- x[, -rmIndex]
  temp <- inner_join(x, ageMeanSub)
  
  ageMeanSubAll <- mean(x$Age1)
  temp$age_mean <- temp$ageMeanSub-ageMeanSubAll
  temp$centertime <- temp$Age1 - temp$ageMeanSub
  temp$centertime2 <- (temp$Age1 - temp$ageMeanSub)^2
  

  return(temp)
}

getBaseSlopeCorr <- function(x, permId){
  tmpT1 <- x[x$time>=2&x$label=='Time1',]
  tmpT2 <- x[x$time>=2&x$label=='Time2',]
  tmpT3 <- x[x$time>=2&x$label=='Time3',]
  tmpT4 <- x[x$time>=2&x$label=='Time4',]
  tmp <- rbind(tmpT1, tmpT2, tmpT3, tmpT4)
  tmp <- tmp[order(tmp$ID, tmp$label),]
  
  dup <- grep('TRUE', duplicated(tmp$ID))
  tmp <- tmp[-(dup-1),]
  selectInd <- rownames(tmp)
  selectID <- tmp[selectInd,'ID']
  idAll <- unique(x$ID)
  selectID_anti <- setdiff(idAll, selectID)
  
  x_long <- tmp
  x_base <- x
  for (r in 1:length(selectID_anti)) {
    x_base <- subset(x_base, x_base$ID!=selectID_anti[r])
  }
  x_base <- x_base[x_base$label=='Time1',]
  
  region_long <- x_long[,50:409]
  region_base <- x_base[,50:409]
  Slope <- (region_long-region_base)/x_long$time
  base_Slope_corr <- diag(cor(region_base, Slope, method = 'spearman'))
  corrPerm <- array(0, c(1000, 360))
  sigArr <- array(0, c(1, 360))
  for (p in 1:ncol(permId)) {
    print(paste0('now performing spin test... ', p, '/', ncol(permId)))
    region_base_perm <- region_base[, permId[, p]]
    base_Slope_corr_perm <- diag(cor(region_base_perm, Slope, method = 'spearman'))
    corrPerm[p,] <- base_Slope_corr_perm
  }
  for (s in 1:ncol(sigArr)) {
    permTmp <- corrPerm[, s]
    if(base_Slope_corr[s]<=0) sigArr[, s] <- length(permTmp[permTmp<=base_Slope_corr[s]])/1000
  }
  returnDf <- data.frame(label=colnames(x)[50:409],
                         corr=base_Slope_corr,
                         pvalue=t(sigArr),
                         slope=colMeans(Slope),
                         base=colMeans(region_base))
  returnDf$label <- str_remove(returnDf$label, '_ROI.*')
  returnDf$label <- str_replace(returnDf$label, '\\.', '-')
  return(returnDf)
}

# nolinear test

lmePipelineTest <- function(df, aim){
  # df <- calDataArea
  # df <- mcData
  # library packages
  # outputPath <- outputPathVal
  if(require(tidyverse)==0){
    install.packages("tidyverse")
    library(tidyverse)
  }
  if(require(lmerTest)==0){
    install.packages("lmerTest")
    library(lmerTest)
  }
  
  #df <- group_20Vol
  
  df[,-c(1:3,417)] <- scale(df[,-c(1:3,417)])
  ModelData <- df
  
  Brain <- colnames(ModelData[,50:409])
  if (aim == 'time') Var <- c('centertime','centertime2')
  if (aim == 'Sex') Var <- c('Sex:centertime', 'Sex:centertime2')
   #????alphaָplasticity??betaָstability
  ResultDf <- matrix(0,length(Brain),3*length(Var))
  AICBIC <- matrix(0, 360, 2*length(Var)+1)
  
  colnames(AICBIC) <- c(paste0(rep(Var, length(Var)), rep(c("AIC", "BIC"), each = length(Var))), "pval")
  colnames(ResultDf) <- c(paste(Var,"pvalue",sep = "-"),paste(Var,"tvalue",sep = "-"),paste(Var,"df",sep = "-"))
  rownames(ResultDf) <- colnames(ModelData[,50:409])
  n_sites <- length(unique(ModelData$site))
  for (vv in 1:length(Var)){
    lmstore <- list()
    # centertime
    if (aim == 'time') {
      Coef <- Var[vv]
      if(n_sites>=2){
        for (aa in 1:360) {
          print(paste(Var[vv],aa))
          ll <- switch(vv,
                       as.formula(paste0(Brain[aa], "~ Sex+age_mean+site+centertime+(1|ID)")),
                       as.formula(paste0(Brain[aa], "~ Sex+age_mean+site+centertime2+centertime+(1|ID)"))
          )
          # lmstore[[aa]] <- lmer(ll,data = ModelData, control = lmerControl(optimizer ="Nelder_Mead"))
          lmstore[[aa]] <- lmer(ll,data = ModelData)
          
        }
        cmd <- paste0("lmstore", Var[vv], " <- lmstore")
        eval(parse(text = cmd))
        
        PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef,]))
        
        ResultDf[,vv] <- PVector[5,]
        ResultDf[,vv+length(Var)] <- PVector[4,]
        ResultDf[,vv+(2*length(Var))] <- PVector[3,]
      }
    }
      if (aim == 'Sex') {
        Coef <- Var[vv]
        if(n_sites>=2){
          for (aa in 1:360) {
            print(paste(Var[vv],aa))
            ll <- switch(vv,
                         as.formula(paste0(Brain[aa], "~ Sex+age_mean+site+centertime+centertime:Sex+Sex+(1|ID)")),
                         as.formula(paste0(Brain[aa], "~ Sex+age_mean+site+centertime2+centertime+centertime:Sex+
                                           centertime2:Sex+Sex+(1|ID)"))
            )
            # lmstore[[aa]] <- lmer(ll,data = ModelData, control = lmerControl(optimizer ="Nelder_Mead"))
            lmstore[[aa]] <- lmer(ll,data = ModelData)
            
          }
          tmp <- str_remove(Var[vv], ':')
          cmd <- paste0("lmstore", tmp, " <- lmstore")
          eval(parse(text = cmd))
          
          PVector <- do.call(cbind,lapply(lmstore, function(x) summary(x)$coefficients[Coef,]))
          
          ResultDf[,vv] <- PVector[5,]
          ResultDf[,vv+length(Var)] <- PVector[4,]
          ResultDf[,vv+(2*length(Var))] <- PVector[3,]
      }

    }
  }
  for (a in 1:360) {
    print(paste0("anova: ", a))
    # a <- 15
    # a <- 1
    ANO <- anova(lmstorecentertime[[a]], lmstorecentertime2[[a]])
    AICBIC[a, 1] <- ANO["lmstorecentertime[[a]]", "AIC"]
    AICBIC[a, 2] <- ANO["lmstorecentertime2[[a]]", "AIC"]
    AICBIC[a, 3] <- ANO["lmstorecentertime[[a]]", "BIC"]
    AICBIC[a, 4] <- ANO["lmstorecentertime2[[a]]", "BIC"]
    AICBIC[a, 5] <- ANO["lmstorecentertime2[[a]]", "Pr(>Chisq)"]
  }
  AICBIC[,5] <- p.adjust(AICBIC[,5],method = 'fdr')
  ResultDf <- data.frame(ResultDf)
  PlotLabel <- str_replace(rownames(ResultDf),"_ROI.*","")
  ResultDf$label <- PlotLabel
  ResultDf$label <- str_replace(ResultDf$label,"\\.","-")
  
  dataType <- unique(str_remove(rownames(ResultDf), ".*ROI_"))
  # file.path(outputPath)
  
  # dir create
  # outdir <- unlist(strsplit(outputPath, "/"))
  # create <- outdir[1]
  # for (d in 2:length(outdir)) {
  #   create <- paste0(create, "/", outdir[d])
  #   if(!dir.exists(create)) dir.create(create)
  # }
  # if(!dir.exists(paste0(create, "/Sex"))) dir.create(paste0(create, "/Sex"))
  # if(!dir.exists(paste0(outputPath, "Sex/", dataType))) dir.create(paste0(outputPath, "Sex/", dataType))
  # write.csv(ResultDf,paste0(outputPath, "Sex/", dataType, "/", dataType, times, "ResUncorrected.csv"),row.names = F)
  # FDR-correct
  VarTmp <- str_remove(Var, ':')
  StatsName <- paste0(dataType,VarTmp,"Stats")
  StatsNameF <- paste0(dataType,VarTmp,"StatsFDR")
  returnDf <- ResultDf
  for (vv in 1:length(Var)) {
    cmd <- paste0(StatsName[vv]," <- data.frame(pvalue=p.adjust(ResultDf[,",vv,"],method = 'BH'),tvalue=ResultDf[,",vv+length(Var),"],label=ResultDf$label)")
    cmd2 <- paste0("returnDf[,", vv,"] <- ", StatsName[vv],"[,1]")
    tvalueCorr <- paste0(colnames(returnDf)[vv+2],".Corr")
    cmd3 <- paste0("returnDf$", tvalueCorr," <- 1")
    eval(parse(text = cmd))
    eval(parse(text = cmd2))
    eval(parse(text = cmd3))
    for (i in 1:nrow(ResultDf)) {
      if (returnDf[i,vv+2]<=0) returnDf[i,tvalueCorr] <- -qt(returnDf[i,vv]/2, returnDf[i,vv+4],lower=FALSE)
      if (returnDf[i,vv+2]>0) returnDf[i,tvalueCorr] <- qt(returnDf[i,vv]/2, returnDf[i,vv+4],lower=FALSE)
    }
  }
  # Extract sig regions
  for (ss in 1:length(Var)) {
    if (nrow(get(StatsName[ss]))>=1) {
      cmd <- paste0(StatsNameF[ss]," <- subset(",StatsName[ss],",",StatsName[ss],"[,1]<0.05)")
      eval(parse(text = cmd))
    }
  }
  # Can not survival after FDR-correct 
  temp <- data.frame(BehaviorDim=StatsNameF,FDR=1)
  for (ss in 1:length(Var)) {
    if(nrow(get(StatsNameF[ss]))==0) {
      # Print regions can not pass FDR-correct
      temp[ss,2] <- 0
      # Threshold the original p-value at 0.001
      cmd1 <- paste0(StatsNameF[ss]," <- subset(ResultDf",",ResultDf[",ss,"]<0.05)")
      # Select p&t values and labels
      cmd2 <- paste0(StatsNameF[ss]," <- ",StatsNameF[ss],"[,c(",ss,",",ss,"+length(Var),",2*length(Var)+1,")]")
      # Rename colnames
      cmd3 <- paste0("colnames(",StatsNameF[ss],") <- ","colnames(",StatsName[ss],")")
      eval(parse(text = cmd1))
      eval(parse(text = cmd2))
      eval(parse(text = cmd3))
    }
  }
  # write.table(temp,paste0(outputPath, "Sex/",dataType, "/", dataType, times, "FDR.txt"),quote = F,row.names = F)
  
  returnDf <- cbind(returnDf, AICBIC)
  return(returnDf)
  
  
}