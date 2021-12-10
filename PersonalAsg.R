#FUNCTION sameShop: true if shop same for product i and j
SameShop <- function(product_i,product_j){
  if(product_i[["shop"]] == product_j[["shop"]]  ){
    return(TRUE)
  }else {
    return(FALSE)
  }
}

#FUNCTION diffBrand: Compare if different brands
DiffBrand <- function(product_i,product_j){
  i_brand = BrandName(product_i,BrandVec)
  j_brand = BrandName(product_j,BrandVec)
  if(is.null(i_brand)| is.null(j_brand)){
    return(FALSE)
  }else{
    if(i_brand != j_brand){
      return(TRUE)
    }else {
      return(FALSE)
    }
  }
}

#FUNCTION BrandName: get brandname for product 
BrandName <- function(product, brand_vector){
  
  brand_name = product$featuresMap$Brand
  
  if(is.null(brand_name)){
    title = product$title 
    words_list =  str_split(title, " ")
    for(w in words_list[[1]]){
      find_brand = match(w, brand_vector)
      if(!is.na(find_brand)){
        brand_name = brand_vector[find_brand]
        break
      }
    }
  }
  if(!is.null(brand_name)){
    brand_name = AdjustName(brand_name)
  }
  return(brand_name)
}

#FUNCTION AdjustName: ensure names are comparable
AdjustName <- function(brand){
  
  brand = tolower(brand)
  if(brand == "lg electronics"){brand = "lg"}
  else if(brand == "sceptre inc"){brand = "sceptre"}
  else if(brand == "jvc tv"){brand = "jvc"}
  return(brand)
  
}

#FUNCTION SplitString: split the string in a loop
SplitString <-function(string){
  substrings = vector()
  for (i in 1:(nchar(string)-1)) {
    substrings = c( substrings, substr(string, i, i+1))
  }
  return(substrings)
}

#FUNCTION CalcSim: calculate the Sørensen–Dice similarity
CalcSim <- function(string1,string2,SUBSTRINGSIZE){
  string1 = tolower(str_replace_all(string1, " ", ""))
  string2 = tolower(str_replace_all(string2, " ", ""))
  
  substrings1 = SplitString(string1)
  substrings2 = SplitString(string2)
  
  nt = 0 # number of character bigrams found in both strings
  for (i in 1:length(substrings1)) {
    substring_i = substrings1[i]
    for (j in 1:length(substrings2)) {
      substring_j = substrings2[j]
      if(substring_i == substring_j ){
        nt = nt + 1
      }
    }
  }
  Similarity = 2*nt / (length(substrings1) +  length(substrings2))
  return(Similarity)
}

#FUNCTION ExMW: all models words from the values of attributes from product p
ExMW <- function(nmk){
  ModelWords = vector()
  for (i in 1:length(nmk)) {
    regResult = gregexpr("(^\\d+(\\.\\d+)?[a-zA-Z\"]+$|ˆ\\d+(\\.\\d+)?$)",nmk[i])
    if( attr(regResult[[1]], 'match.length') == -1 ){
      next
    }
    wordlist = regmatches(nmk[i], regResult)
    for (j in 1:length(wordlist[[1]])) {
      Word = wordlist[[1]][j]
      if(!any(ModelWords==Word)){
        ModelWords = c(ModelWords, Word)
      }
    }
  }
  return(ModelWords)
}

#FUNCTION MW: percentge of matching model words from two sets of model words
MW <-function(set1, set2){
  if(length(set1)==0|length(set2)==0){
    return(0)
  }
  
  counter = 0
  for (i in 1:length(set1)) {
    for (j in length(set2)) {
      if(set1[i] == set2[j]){
        counter = counter+1
      }
    }
  }
  percentage = counter*2 /(length(set1) + length(set2)) * 100
  
  return(percentage)
}

#FUNCTION MSM: Multi-component Similarity Method
MSM <- function(ProductData, CandidateMatrix, MU, GAMMA){
  n = length(ProductData)
  disSim = matrix(data = 0, nrow = n, ncol = n)
  disSim[lower.tri(disSim)] <- NA 

  for (i in 1:(n-1)) {
    for (j in i:n) {
      #writeLines(paste("(i,j)", i,j))
      product_i = TVData[[i]]
      product_j = TVData[[j]]
      #Check if same shop, different brand or not candidate
      if (SameShop(product_i,product_j) | CandidateMatrix[i,j] == 0 |  DiffBrand(product_i,product_j)){
        disSim[i,j] = Inf
      }else{
        sim = 0
        avgSim = 0
        m = 0 #number of matches
        w = 0 # weight of matches
        nmkf_i = TVData[[i]][["featuresMap"]] 
        nmkf_j = TVData[[j]][["featuresMap"]]
        nmk_i = nmkf_i
        nmk_j = nmkf_j
        #apply KVP
        for (q in 1:length(nmkf_i)) {
          for (r in 1:length(nmkf_j)) {
            keySim = CalcSim(names(nmkf_i)[q],names(nmkf_j)[r])
            if(keySim > GAMMA){
              valueSim = CalcSim(nmkf_i[[q]],names(nmkf_j)[[r]])
              weight = keySim
              sim = sim + weight*valueSim
              m = m + 1
              w = w+weight
              nmk_i = nmk_i[-q]
              nmk_j = nmk_j[-r]
            }
          }
        }
        if(w > 0){
          avgSim = sim / w
        }
        #apply MW
        if(length(nmk_i==0)|length(nmk_j==0)){
          mwPerc = 0
        }else{
          mwPerc = MW(ExMW(nmk_i),ExMW(nmk_j))
        }
        #calculate the new weight
        theta1 = (1-MU)*m/ min (length(nmkf_i),length(nmkf_j))
        theta2 = 1-MU-theta1
        hSim = theta1*avgSim + theta2*mwPerc
        disSim[i,j] = 1 - hSim
      }
    }
  }
  return(disSim)
}

#FUNCTION DuplicateSelection: examine dissemilarity matrix of MSM
DuplicateSelection <- function(MSM, EPSILON){
  n = dim(MSM)[2]
  DuplicateMatrix = matrix(data = 0, nrow = n, ncol = n)
  DuplicateMatrix[lower.tri(DuplicateMatrix)] <- NA 
  for (i in 1:n) {
    for (j in i:n) {
      if(MSM[i,j]<EPSILON)
        DuplicateMatrix[i,j] = 1
    }
  }
  return(DuplicateMatrix)
}

#FUNCTION GetTvTitles: calculate all integral 
BandSizes <-function(HASHNUMBER ){
  temp = c(1:HASHNUMBER)
  bandsizes = 1
  for (i in 2:length(temp)) {
    if(HASHNUMBER %%temp[i] ==0)
      bandsizes = c(bandsizes,temp[i])
  }
  return(bandsizes)
}

#FUNCTION GetTvTitles: Extract all TVtitles in one list
GetTVTitles <- function(TVData){
  #Get a list of product titles
  allTVtitles <- list()
  for (i in 1:length(TVData)){
    allTVtitles[i] = TVData[[i]][["title"]]
  }
  return(allTVtitles)
}

#FUNCTION BrandVector: get all the brands
BrandVector <- function(ProductData){
  
  brand_vector = c()
  for(i in 1:length(ProductData)){
    brand = ProductData[[i]]$featuresMap$Brand
    if(!is.null(brand)){
      brand_vector = append(brand_vector, brand)
    }
  }
  return(unique(brand_vector))
}

#FUNCTION GetBinaryVec: Get a uniform vector representation for each product
GetBinaryVec <- function(TVtitles,BrandVector){
  # -Input -> allTVtitles: list of product titles 
  # -Output -> BinaryVecComplete: matrix of binary vector for each product where row are the model words and columns the product 
  
  #get set of all model words from titles of all product descriptions
 
  Dictionary = vector()
  for (i in 1:length(TVtitles)) {
    Titlereg = gregexpr("[a-zA-Z0-9](([0-9]+[^0-9, ]+)|([^0-9, ]+[0-9]+))[a-zA-Z0-9] ]*",TVtitles[i])
    if( attr(Titlereg[[1]], 'match.length')[1] == -1 ){
      next
    }
    Titlewords = regmatches(TVtitles[i],Titlereg)
    for (j in 1:length(Titlewords[[1]])) {
      Word = Titlewords[[1]][j]
      if(!any(Dictionary==Word)){
        Dictionary = c(Dictionary, Word)
      }
    }
  }
  Dictionary = c(Dictionary,BrandVector)
  #print(Dictionary)
  
  #get binary vec for each product
  BinaryVecComplete = matrix(data = 0, nrow = length(Dictionary), ncol = length(TVtitles))
  for (i in 1:length(TVtitles)) {
    productTitle = TVtitles[[i]]
    for (j in 1:length(Dictionary)) {
      if(grepl(Dictionary[j],productTitle)){
        BinaryVecComplete[j,i] = 1
      } 
    }
  }
  

  return(BinaryVecComplete)
}

#FUNCTION Hashfunction: compute hash values
Hashfunction <- function(a,b,r){
  return(  (a*r+b) %%  PRIME)
}

#FUNCTION MinHashingEfficient:  apply minhashing 
MinHashingEfficient <- function(BinaryMatrix,HASHNUMBER){
  # Input -> BinaryMatrix uniform matrix representation for each product
  # Output -> signature_matrix

  #create hashfunction coefficients
  dim_BM = dim(BinaryMatrix)  #check if prime number sufficient large
  if(dim_BM[1]>PRIME){
    print("Warning: prime number is smaller then adviced.")
  }
  
  hashcoefficients = vector() 
  for (i in 1:HASHNUMBER) {
    a = as.integer(runif(1, min = 1, max = PRIME))
    b = as.integer(runif(1, min = 0, max = PRIME))
    #b = 0
    pair = c(a,b)
    hashcoefficients = cbind(hashcoefficients,pair)
  }
  
  #print(hashcoefficients)
  #fill in signature matrix
  signature_matrix = matrix(data = Inf, nrow = HASHNUMBER, ncol = dim_BM[2]) 
  hashfunctionresult = vector(mode = "integer", length = HASHNUMBER)
  for (r in 1:dim_BM[1]) {
    #compute hashfunction
    for (h in 1:HASHNUMBER) {
      hashfunctionresult[h] = Hashfunction(hashcoefficients[1,h],hashcoefficients[2,h],r)
    }
    #print(hashfunctionresult)
    for (c in 1:dim_BM[2]) {
      if(BinaryMatrix[r,c] == 1){
        for (h in 1:HASHNUMBER) {
          if( hashfunctionresult[h] < signature_matrix[h,c] )
            signature_matrix[h,c] = hashfunctionresult[h]
        }
      }
    }
  }
  
  return(signature_matrix)
}

#FUNCTION LSH: Local Sensitivity Hashing
LSH <- function(signature_matrix,b){
  # Input -> signature_matrix 
  #       -> b: is the amount of bands chosen to divide
  # Output -> CandidateMatrix binary matrix with potential product duplicates
  dim_SM = dim(signature_matrix)
  r = dim_SM[1]/b
  #print(r)
  CandidateMatrix = matrix(data = 0, nrow = dim_SM[2], ncol = dim_SM[2])
  CandidateMatrix[lower.tri(CandidateMatrix)] <- NA 
  band_matrix = signature_matrix[1:r,]
  for (i in 1:b){
    #assign first band or next
    if(i == 1 ){
      band_matrix = signature_matrix[1:r,]
    }else{
      band_matrix = signature_matrix[((i-1)*r+1):(r*i),]
    }
    
    #for every column check all to the right if equal
    for (j in 1:(dim_SM[2]-1)) {
      col_fixed = band_matrix[,j]

      for (k in (j+1):dim_SM[2]) {
        #writeLines(paste("compare", col_fixed , "with", band_matrix[,k]))
        if(identical(col_fixed, band_matrix[,k]) ){
          #writeLines(paste("found duplciate", j,k))
          CandidateMatrix[j,k] = 1 #potential duplicate product j and k
        }
      }
    }
  }
  return(CandidateMatrix)
}

#FUNCTION LSH: Local Sensitivity Hashing
LSHonerow <- function(signature_matrix,b){
  # Input -> signature_matrix 
  #       -> b: is the amount of bands chosen to divide
  # Output -> CandidateMatrix binary matrix with potential product duplicates
  dim_SM = dim(signature_matrix)
  r = 1
  #print(r)
  CandidateMatrix = matrix(data = 0, nrow = dim_SM[2], ncol = dim_SM[2])
  CandidateMatrix[lower.tri(CandidateMatrix)] <- NA 
  band_matrix = vector()
  for (i in 1:b){
    band_matrix = signature_matrix[i,]
    #for every column check all to the right if equal
    for (j in 1:(dim_SM[2]-1)) {
      col_fixed = band_matrix[j]
      for (k in (j+1):dim_SM[2]) {
        #writeLines(paste("compare", col_fixed , "with", band_matrix[,k]))
        if(identical(col_fixed, band_matrix[k]) ){
          #writeLines(paste("found duplciate", j,k))
          CandidateMatrix[j,k] = 1 #potential duplicate product j and k
        }
      }
    }
  }
  return(CandidateMatrix)
}

#FUNCTION BootstrapMethod: generate bootstrapdata
BootStrapData <-function(TVData){
  
  #list_of_products_train = unique(sample(TVData, 1624, replace = T))
  list_of_products_train = unique(sample(TVData, length(TVData), replace = T))
  list_of_products_test = setdiff(TVData, list_of_products_train)
  return(list('TVData_BS_Train' = list_of_products_train, 'TVData_BS_Test' = list_of_products_test))
}

#FUNCTION EvaluationMatrix: construct the tp,np,fp,fn values
EvaluationMatrix = function(neighbor_matrix, list_of_products){
  #print(dim(neighbor_matrix))
  #print(length(list_of_products))
  TP = 0; FP = 0; TN = 0; FN = 0
  n_products = length(list_of_products)
  
  # make a adjusted dissimilarity matrix 
  for (i in 1:n_products){
    # compare product i with every product j 
    for(j in 1:n_products){
      # if this holds than product i and j are similar 
      if(j>i && neighbor_matrix[i,j]==1){
        # find out if they are truly similar
        if(list_of_products[[i]]$modelID == list_of_products[[j]]$modelID){
          TP = TP + 1
        } else{FP = FP + 1}
      }
      # else the products are not classified as similar
      else if(j>i){
        if(list_of_products[[i]]$modelID == list_of_products[[j]]$modelID){
          FN = FN + 1
        } else{TN = TN + 1}
      }
    }
  }
  
  writeLines(paste("confustion TP FP TN FN", TP,FP,TN,FN))
  return(list('TP'= TP, 'FP'=FP, 'TN'=TN, 'FN'=FN))
  
}

#FUNCTION Evaluation: compute evaluation values
Evaluation <-function(EvaluationMatrix,B_LSH,n,band){
  if(B_LSH == TRUE){
    PQ = ((EvaluationMatrix$TP) / (n*band))
  }else{
    PQ = (EvaluationMatrix$TP) / (EvaluationMatrix$TP + EvaluationMatrix$FP)
  }
  
  PC = EvaluationMatrix$TP / (EvaluationMatrix$TP + EvaluationMatrix$FN)
  F1 = 2/((1/PQ)+(1/PC))
  FC = ((sum(1:(n-2))+n-2)*band) /((sum(1:(n-2))+n-2)*HASHNUMBER)
  #print(EvaluationMatrix$TP)
  #print((EvaluationMatrix$TP + EvaluationMatrix$FP))
  #print((EvaluationMatrix$TP + EvaluationMatrix$FN))
  
  return(list('PQ'= PQ, 'PC'=PC, 'F1'=F1, 'FractionComp' = FC))
}

#FUNCTION EstimateAlgorithm: execute the algorithm
EstimateAlgorithm <- function(TVData_BS){
  bandsizes = BandSizes(HASHNUMBER)
  print(bandsizes)
  allTVtitles = GetTVTitles(TVData_BS)
  BinaryMatrix = GetBinaryVec(allTVtitles,BrandVec)
  SignatureValues = MinHashingEfficient(BinaryMatrix,HASHNUMBER)
  
  TotalEvaluation_BS <- matrix(data = NA, nrow = length(bandsizes), ncol = 5)
  for (b in 1:length(bandsizes)) {
    band = bandsizes[b]
    writeLines(paste("amount of bands", band))
    #print(b)
    #print(length(bandsizes))
    if(b == (length(bandsizes))){
      PotentialDup = LSHonerow(SignatureValues,band)
    }else{
      PotentialDup = LSH(SignatureValues,band)
    }
    ConfusionMatLSH = ConfusionMatrix(PotentialDup, TVData_BS)
    EvaluationLSH  = Evaluation(ConfusionMatLSH,TRUE,length(TVData_BS),band)
    #print(EvaluationLSH)
    TotalEvaluation_BS[b,1:4] = unlist(EvaluationLSH,)
    
    MSMVector = MSM(TVData_BS, PotentialDup, MU, GAMMA)
    DuplicateMatrix = DuplicateSelection(MSMVector,EPSILON)
    ConfusionMatMSM = ConfusionMatrix(DuplicateMatrix, TVData_BS)
    EvaluationMSM  = Evaluation(ConfusionMatMSM,FALSE,length(TVData_BS),band)
    TotalEvaluation_BS[b,5] =  EvaluationMSM$F1
    
    colnames(TotalEvaluation_BS) <- c("PQ","PC","F1Star",'FractionComp', "F1")
  }
  print(TotalEvaluation_BS)
  return(list('Evalation' = TotalEvaluation_BS, 'MU' = MU, 'GAMMA' = GAMMA, 'EPSILON' = EPSILON))
}

### RUNNING CODE
#Define required packages and constants
library(rjson) # import data
library(stringr) # split strings
library("ggplot2")

PRIME <- 1637 #317 #53 #1637 #317  #109  #17 #22447 #2269
SUBSTRINGSIZE <- 3 #substring when using Sørensen–Dice similarity
BOOTSTRAPCOUNT <- 5 #do 5 bootstraps

#Import data and seperate shops
TVMergedData <- fromJSON(file="TVs-all-merged.json") 
TVData <- list()
counter = 1
for (i in 1:length(TVMergedData)){
  temp = length(TVMergedData[[i]])
  if(temp != 1){
    for (j in 1:temp) {
      TVData[counter+j-1] = TVMergedData[[i]][j]
    }
    counter = counter + temp
  }else{
    TVData[counter] = TVMergedData[[i]]
    counter = counter+1
  }
}
remove(TVMergedData,i,j,temp,counter)
#TVData = TVData[1:300] ##ADJUST when testing

#Generate GridSearch plain
gammavec =  seq(0.5, 0.7, by=0.1)
muvec = seq(0.5, 0.7, by =0.1)
epsilonvec = seq(0.95, 1.05, by=0.1)
gridsearchvaluescoef = matrix(data = NA, nrow = (length(gammavec)*length(muvec)*length(epsilonvec)), ncol = 3)
counter = 1
for (g in 1:length(gammavec)) {
  for (m in 1:length(muvec)) {
    for (e in 1:length(epsilonvec)) {
      gridsearchvaluescoef[counter,] = c(gammavec[g],muvec[m],epsilonvec[e])
      counter = counter + 1
    }
  }
}
remove(e,epsilonvec,g,gammavec,m,muvec,counter)
  
#Apply algorithm to train set to find optimal parameters
TVData_BSTest_all = list()
for (g in 1:dim(gridsearchvaluescoef)[1]) {
  GAMMA = gridsearchvaluescoef[g,1]
  MU = gridsearchvaluescoef[g,2]
  EPSILON = gridsearchvaluescoef[g,3]
  
  for (i in 1:BOOTSTRAPCOUNT){
    
    TVData_BS_complete = BootStrapData(TVData)
    TVData_BSTrain = TVData_BS_complete$TVData_BS_Train
    TVData_BSTest_all[[length(TVData_BSTest_all)+1]] =  TVData_BS_complete$TVData_BS_Test #safe test set for next step
    writeLines(paste("BOOTSTRAP", i,"on", length(TVData_BSTrain),"products"))
    BrandVec = BrandVector(TVData_BSTrain)
    HASHNUMBER <- as.integer(length(TVData_BSTrain)*0.5)
    TotalEvaluation_BS = EstimateAlgorithm (TVData_BSTrain)
    #TotalEvaluation_list[[length(TotalEvaluation_list)+1]] = TotalEvaluation_BS
    if( i ==1 ){
      TotalEvaluation_Averaged_Train = TotalEvaluation_BS$Evalation
    }else{
      TotalEvaluation_Averaged_Train = TotalEvaluation_Averaged_Train+TotalEvaluation_BS$Evalation
    }
    #print("total evaluation is now ")
    #print(TotalEvaluation_Averaged)
  }
  TotalEvaluation_Averaged_Train[,1] = TotalEvaluation_Averaged_Train[,1] / BOOTSTRAPCOUNT
  TotalEvaluation_Averaged_Train[,2] = TotalEvaluation_Averaged_Train[,2] / BOOTSTRAPCOUNT
  TotalEvaluation_Averaged_Train[,3] = TotalEvaluation_Averaged_Train[,3] / BOOTSTRAPCOUNT
  TotalEvaluation_Averaged_Train[,5] = TotalEvaluation_Averaged_Train[,5] / BOOTSTRAPCOUNT
}
print("Last evaluation ")
print(TotalEvaluation_Averaged_Train)
  
#Calculate optima values for GAMMA, MU, EPSILON and apply on test set
GAMMA <- 0.6
MU <- 0.5
EPSILON <- 1.05

#Apply algorithm to test set for evaluation  
for (i in 1:BOOTSTRAPCOUNT){
    TVData_BSTest =  TVData_BSTest_all[i]
    writeLines(paste("BOOTSTRAP", i,"on", length(TVData_BSTest),"products"))
    BrandVec = BrandVector(TVData_BSTest)
    HASHNUMBER <- as.integer(length(TVData_BSTest)*0.5)
    TotalEvaluation_BS = EstimateAlgorithm (TVData_BSTest)
    if( i ==1 ){
      TotalEvaluation_Averaged = TotalEvaluation_BS$Evalation
    }else{
      TotalEvaluation_Averaged = TotalEvaluation_Averaged+TotalEvaluation_BS$Evalation
    }
    #print("total evaluation is now ")
    #print(TotalEvaluation_Averaged)
  }
TotalEvaluation_Averaged = TotalEvaluation_Averaged / BOOTSTRAPCOUNT
print("Last evaluation ")
print(TotalEvaluation_Averaged)
    
#Draw ggplot2 scatterplot with smooth curve
df = as.data.frame(TotalEvaluation_Averaged)

names(df) = c("Pair Quality", 'Pair Completeness', 'F1*-Measure', 
                      'Fraction of Comparisons', 'F1-Measure')
                                
pq <- ggplot(df, aes(`Fraction of Comparisons`,`Pair Quality`)) + 
  geom_line()
pq + theme_classic()

pc <- ggplot(df, aes(`Fraction of Comparisons`,`Pair Completeness`)) + 
    geom_line()
pc + theme_classic()

f1s <- ggplot(df, aes(`Fraction of Comparisons`,`F1*-Measure`)) + 
      geom_line()
f1s + theme_classic()

f1 <- ggplot(df, aes(`Fraction of Comparisons`,`F1-Measure`)) + 
      geom_line()
f1 + theme_classic()







