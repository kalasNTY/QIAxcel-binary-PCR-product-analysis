#Test if specified wells are present in the raw data. If not throw error
testly <- function(table,chara) {
   if(is.numeric(table[,chara])==TRUE){
     
   }
}  

#Test if tableratio contains at least one row
testly2 <- function(table) {
  y <- apply(table,2,function(x) if(x[1]==x[2]) c(x[1],NA) else c(x[1],x[2]))
  invy <- t(y)
  invy <- data.frame(invy[,1],as.numeric(as.character(invy[,2])))
}  

#Function that generates a bar graph using mean and SD from replicates
plotty <- function(table,output) {
  require(ggplot2)
  png(paste("barplot-sd-",output, ".png", sep=""))
  plot <- ggplot(table, aes(x = Sample, y = AVR)) +  
    geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
    geom_errorbar(aes(ymin=AVR-SD, ymax=AVR+SD)) +
    ggtitle("QIAxcel: exon inclusion ratio") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank())
  dev.off()
  print(plot)
}
#function that inverses and formats the output table
formatty <- function(table,rownumber){
  table <- table[,-c(1)]
  table <- t(table)
  colnames(table) <- c("AVR","SD")
  rownames(table) <- c(1:rownumber)
  table <- as.data.frame(table)
  table$Sample <- c(1:rownumber)
  table <- na.omit(table[,c(3,1,2)])
}


#Takes qiaxcel data string from csv file (e.g. "data.csv")  and looks for amplification products
#of defined length (prdshort,prdlong) in a range of wells defined and outputs the ration of the longer
#versus the sum of both products (long + short)
qcell <- function(flname, wells = "all", prdshort, prdlong, margin=6, repl = 1, round = 6){

  
  #Import csv file and store in data frame
  dfqia <- read.csv(flname,header=FALSE)  
  nrows <- nrow(dfqia)
  
  #Find the column ^that contains "Property" in the top and then find row where size is located.
  #Conc (ng/ul) is assumed to be the row below size (as there are two rows called "conc" grepl is suboptimal)
  prop <- which(apply(dfqia, 2, function(x) any(grepl("Property", x))))
  prop <- as.numeric(prop)
  wellwell <- which(apply(dfqia, 2, function(x) any(grepl("Well", x))))
  wellwell <- as.numeric(wellwell)
  prop2 <- which(apply(dfqia, 1, function(x) any(grepl("size", x))))
  prop2 <- prop2[1]-1
  prop3 <- prop2+1

   
  #create new data frame containing amplification product lenght and concentration
  dfqiawell <- split(dfqia,dfqia[,wellwell])
  outputbp <<- sapply(dfqiawell, function(x) ifelse(x[prop2,]!=0,x[prop2,],NA))
  outputconc <<- sapply(dfqiawell, function(x) ifelse(x[prop3,]!=0,x[prop3,],NA))
  outputbp <- outputbp[-c(1:5),]
  outputconc <- outputconc[-c(1:5),]
  outputbp <- as.data.frame(outputbp)
  outputconc <- as.data.frame(outputconc)

  #find the amplification products in the data set +/- a margin
  lowmarg <- prdshort - margin
  himarg <- prdshort + margin
  products <- c(1:nrow(outputbp))
  nwell <- c(1:ncol(outputbp))

  h <- 1
  
  for (y in nwell) {
 
    n <- 1
    for(i in products){
      
      value <- outputbp[n,h]
      if(is.na(value)){
        
        next()
        }
      if(value >= lowmarg && value <= himarg){
        outputbp[n,h] <- prdshort
      } 
      n <- n+1
    }
  h <- h+1
  }
  
  lowmarg <- prdlong - margin
  himarg <- prdlong + margin
  h <- 1
  
  for (y in nwell) {
    
    n <- 1
    for(i in products){
      
      value <- outputbp[n,h]
      if(is.na(value)){
        
        next()
      }
      if(value >= lowmarg && value <= himarg){
        outputbp[n,h] <- prdlong
      } 
      n <- n+1
    }
    h <- h+1
  }
  
  #Create a vector of wells to be analyzed based on the input
    if(wells[1]=="all"){
      wells <- c(sprintf("A%02d", 1:12),sprintf("B%02d", 1:12),sprintf("C%02d", 1:12),sprintf("D%02d", 1:12),
      sprintf("E%02d", 1:12),sprintf("F%02d", 1:12),sprintf("G%02d", 1:12),sprintf("H%02d", 1:12))
      } else if(grepl('row',wells[1])==TRUE) {
        rowno <- as.numeric(substr(wells,4,4))
        alfabet <- c("A%02d","B%02d","C%02d","D%02d","E%02d","F%02d","G%02d","H%02d")
        key <- alfabet[rowno]
        wells <- sprintf(key, 1:12)
    } else if(grepl('column',wells[1])==TRUE) {
      colno <- as.numeric(substr(wells,7,8))
      wells <- sprintf("%s%02d", c(LETTERS[1:8]),colno)
    }
  
  #Go through all wells specified and calculate the product ratios
  tableratio <- matrix(nrow=2)
  n <- 1
  for (i in wells) {
    well <- wells[n]
  
    #Test If the well is  present in the raw data (testly): if not skip well and go to next one
    possibleError <- tryCatch(
      testly(outputbp,well),
      error=function(e) e
    )
    #If the specified well is present in the raw data, then calculate ratios
    if(!inherits(possibleError, "error")){
      
    
      bp <- as.numeric(outputbp[,well])
      conc <- as.numeric(outputconc[,well])
      wellpeaks <- cbind(bp,conc)
      wellpeaks <- as.data.frame(wellpeaks)
      colnames(wellpeaks) <- c("bp","conc")
      wellpeaks <- na.omit(wellpeaks)
      shortconc <- wellpeaks[wellpeaks$bp==prdshort,"conc"]
      longconc <- wellpeaks[wellpeaks$bp==prdlong,"conc"]
      ratio <- longconc/(longconc+shortconc)
      ratio <- round(ratio,round)
    
      colratio <- c(well,ratio)
      tableratio <- data.frame(tableratio,colratio)
      n <- n+1
    }
  }
  tableratio <- tableratio[,-c(1)]
  
  #if no wells are present in tableratio then throw error and omit formating
   possError <- tryCatch(
    testly2(tableratio),
    error=function(e) e
  )
 
   #If the specified well is present in the raw data, then calculate ratios
  if(!inherits(possError, "error")){
    tableratio <- apply(tableratio,2,function(x) if(x[1]==x[2]) c(x[1],NA) else c(x[1],x[2])) 
    invtable <- t(tableratio)
    tableratio <- data.frame(invtable[,1],as.numeric(as.character(invtable[,2])))
    colnames(tableratio) <- c(1:ncol(tableratio))
    rownames(tableratio) <- c(1:nrow(tableratio))
    
  } else if(inherits(possError, "error")){
    
    return(tableratio)
  }
   
   
  #Calculate average and standard deviation for replicates if repl > 1
  if(repl==2){
    samptab <- matrix(nrow=2)
    halftab <- nrow(tableratio)/2
    samples <- c(1:halftab)
    x1 <- 1
    x2 <- 2
    for (g in samples) {
      averdupl <- tableratio[c(x1,x2),2]
      aver <- mean(averdupl)
      aver <- round(aver,round)
      std <- sd(averdupl)
      std <- round(std,round)
      colsamp <- c(aver,std)
      samptab <- data.frame(samptab,colsamp)
      x1 <- x1+2
      x2 <- x2+2
    }
     
      samptab <- formatty(samptab,halftab)
      plotty(samptab,flname)
      write.csv(samptab, file = paste("qcell-output-",repl,"-icates-",flname, sep=""))
      return(samptab)

  } else if(repl==3){
    samptab <- matrix(nrow=2)
    halftab <- nrow(tableratio)/3
    samples <- c(1:halftab)
    x1 <- 1
    x2 <- 2
    x3 <- 3
    for (g in samples) {
      averdupl <- tableratio[c(x1,x2,x3),2]
      aver <- mean(averdupl)
      aver <- round(aver,round)
      std <- sd(averdupl)
      std <- round(std,round)
      colsamp <- c(aver,std)
      samptab <- data.frame(samptab,colsamp)
      x1 <- x1+3
      x2 <- x2+3
      x3 <- x3+3
    }
    
    samptab <- formatty(samptab,halftab)
    
    plotty(samptab,flname)
    write.csv(samptab, file = paste("qcell-output-",repl,"-icates-",flname, sep=""))
    return(samptab)
    
   } else
     
     return(tableratio)
}