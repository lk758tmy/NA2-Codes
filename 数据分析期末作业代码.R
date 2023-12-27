library(jsonlite)
setwd("/home/758_lk/Desktop/数据分析/")
lfile <- list.dirs(path="data_accounts")
#初始化

gameList <- vector()
for(i in seq_along(lfile)){
  if(i==1) next
  data <- readLines(con<-file(paste0(lfile[i],"/playing.txt"),encoding="GBK"))
  data <- fromJSON(data)
  close(con)
  for(j in 1:nrow(data[[1]])){
    gameName <- data[[1]][j,1]
    if(!(gameName %in% gameList))
      gameList <- c(gameList,gameName)
  }
}
gameCnt <- length(gameList)
#得到所有棋的名字

playCntTable <- data.frame()
winRateTable <- data.frame()
for(i in seq_along(lfile)){
  if(i==1) next
  data <- readLines(con<-file(paste0(lfile[i],"/playing.txt"),encoding="GBK"))
  data <- fromJSON(data)
  close(con)
  playCnt <- rep(0,gameCnt)
  winRate <- rep(0,gameCnt)
  maxPlayCnt=0;
  for(j in 1:nrow(data[[1]])){
    k <- match(data[[1]][j,1],gameList)
    maxPlayCnt=max(maxPlayCnt,data[[1]][j,2])
    playCnt[k] <- data[[1]][j,2]
    winRate[k] <- data[[1]][j,3]/data[[1]][j,2]
  }
  for(j in 1:gameCnt){
    playCnt[j] <- playCnt[j]/maxPlayCnt
  }
  playCntTable <- rbind(playCntTable,playCnt)
  winRateTable <- rbind(winRateTable,winRate)
}
colnames(playCntTable) <- gameList
colnames(winRateTable) <- gameList
jointTable <- rbind(playCntTable,winRateTable)
playCntTable <- t(playCntTable)
winRateTable <- t(winRateTable)
jointTable <- t(jointTable)
#得到所有棋手的对局信息

#以上就完成了所有数据读入的部分，得到两个二维列表

library(stats)
#wss <- (nrow(data)-1)*sum(apply(winRateTable,2,var))
#for (i in 2:15) 
#  wss[i] <- sum(kmeans(winRateTable,centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Cluster",ylab="WSS")
result <- kmeans(playCntTable,centers=7)
print(result)
result <- kmeans(winRateTable,centers=7)
print(result)
result <- kmeans(jointTable,centers=7)
print(result)
#K均值算法（无监督）

classMGStyle=list(
  list(1,1),list(2,1),list(3,1),list(4,1),list(5,1),list(6,1),list(7,2),
  list(10,2),list(11,2),list(16,1),list(19,1),list(20,1),list(21,2),list(22,1),
  list(28,2),list(30,2),list(35,1)
)
#监督学习算法的初始化，手动切分打标记

playCntTrain <- data.frame()
winRateTrain <- data.frame()
classTrain <- c()
rowName <- c()
for(i in 1:length(classMGStyle)){
  id <- unlist(classMGStyle[[i]][1])
  rowName <- c(rowName,gameList[id])
  classTrain <- c(classTrain,classMGStyle[[i]][2])
  playCntTrain <- rbind(playCntTrain,playCntTable[id,])
  winRateTrain <- rbind(winRateTrain,winRateTable[id,])
}
rownames(playCntTrain) <- rowName
rownames(winRateTrain) <- rowName
classTrain <- unlist(classTrain)
jointTrain <- cbind(playCntTrain,winRateTrain)

playCntTask <- data.frame()
winRateTask <- data.frame()
trainedGames <- rowName
rowName <- c()
for(i in 1:gameCnt){
  if(i %in% trainedGames) next
  rowName <- c(rowName,gameList[i])
  playCntTask <- rbind(playCntTask,playCntTable[i,])
  winRateTask <- rbind(winRateTask,winRateTable[i,])
}
rownames(playCntTask) <- rowName
rownames(winRateTask) <- rowName
jointTask <- cbind(playCntTask,winRateTask)

colnames(playCntTrain) <- (1:ncol(playCntTrain))
colnames(playCntTask) <- (1:ncol(playCntTask))
colnames(winRateTrain) <- (1:ncol(winRateTrain))
colnames(winRateTask) <- (1:ncol(winRateTask))
colnames(jointTrain) <- (1:ncol(jointTrain))
colnames(jointTask) <- (1:ncol(jointTask))
#监督学习算法的初始化，生成训练集和待测集

library(class)
result <- knn(train=playCntTrain,test=playCntTask,cl=classTrain,k=5)
print(result)
result <- knn(train=winRateTrain,test=winRateTask,cl=classTrain,k=5)
print(result)
result <- knn(train=jointTrain,test=jointTask,cl=classTrain,k=5)
print(result)
#先来个k近邻

library(e1071)
model <- svm(playCntTrain,classTrain,kernel="linear",scale=FALSE)#,type="C")
#print(predict(model,playCntTrain))
result <- predict(model,playCntTask)
print(result)
model <- svm(winRateTrain,classTrain,kernel="linear",scale=FALSE)#,type="C")
result <- predict(model,winRateTask)
print(result)
model <- svm(jointTrain,classTrain,kernel="linear",scale=FALSE)#,type="C")
result <- predict(model,jointTask)
print(result)
#再来个SVM

#最后来个决策树
library(mlpack)
result <- decision_tree(training=playCntTrain,labels=as.matrix(classTrain),
                       minimum_leaf_size=4,test=playCntTask)
print(t(result$probabilities))
result <- decision_tree(training=winRateTrain,labels=as.matrix(classTrain),
                       minimum_leaf_size=4,test=winRateTask)
print(t(result$probabilities))
result <- decision_tree(training=jointTrain,labels=as.matrix(classTrain),
                        minimum_leaf_size=4,test=jointTask)
print(t(result$probabilities))
