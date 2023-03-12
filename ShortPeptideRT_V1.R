
#------------------------------------------------------------------
# Title: ShortPeptideRT 
# sub: Predicts the retention time of short peptides, including homologous structures.  
#
# Author: Boudewijn Hollebrands
# Date: 12 March 2023
# Version: 1.0
#
#
# Adress: https://github.com/BoudewijnHollebrands/ShortPeptideRT
# 
#
# Description: For each peptide sequence a set of descriptors is calculated
# that are used to predict the retention time using a pre-trained model.
#  
# 1. A set of non-sequence specific descriptors are calculated.
# 2. The ASP-descriptors are calculated  
# 3. The descriptors are used as input for a pre-trained support vector regression model.
#---------------------------------------------------------------------

# Change the sequences between "" to predict other Retention times. 
# for example: input_peptide <- c("WNPV","WNVP") 

input_peptide <- c("WNPV","WNVP")  #For these sequences the RT is calculated.


#Amino acid indices obtained from AAindexDB
#https://www.genome.jp/entry/aaindex:BULH740102
#https://www.genome.jp/entry/aaindex:QIAN880111
aa_index <- data.frame(aa = c("Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr"),
                       code = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
                       BULH740102 = c("0.691","0.624","0.558","0.632","0.756","0.592","0.646","0.809","0.767","0.842","0.709","0.596","0.730","0.649","0.728","0.594","0.655","0.777","0.743","0.743"),
                       QIAN880111 = c("0.21","-0.12","-0.58","-0.23","-0.06","-0.15","0.37","0.31","0.28","0.70","0.61","-0.04","-1.03","0.13","0.07","-0.28","-0.25","0.00","0.21","0.16")
)

#convert columns to numeric
aa_index$BULH740102 <- as.numeric(aa_index$BULH740102)
aa_index$QIAN880111 <- as.numeric(aa_index$QIAN880111)


# Load pre-trained prediction model
fit_svm = readRDS("060123_fit_svm.rda")


# calculate descriptors for input_peptide 

for (x in 1:length(input_peptide)){     #for each peptide sequence in input_peptide

######## 1. calculate non-sequence descriptors #######
length <- Peptides::lengthpep(input_peptide[x]) 
mw <- Peptides::mw(input_peptide[x])
LogP <- Peptides::hydrophobicity(input_peptide[x], scale = "AbrahamLeo")
pI <- Peptides::pI(input_peptide[x])
vhse <- Peptides::vhseScales(input_peptide[x])
input <- cbind(Length= length,Mw= mw,LogP=LogP,pI= pI, t(unlist(vhse)))

######### 2. calculate ASP descriptors #########

# peptide N-terminal values
N_term <- substr(input_peptide[x],1,1)
N1_bul <- aa_index[aa_index$code == N_term,3] 
N1_Qian <- aa_index[aa_index$code == N_term,4] 

# peptide C-terminal values
C_term <- substr(input_peptide[x],length,length)
C1_bul<-aa_index[aa_index$code == C_term,3] 
C1_Qian<-aa_index[aa_index$code == C_term,4] 

# average value of amino acids between the N-terminal and C-terminal amino acids
mid <- substr(input_peptide[x],2,(length-1))
mid_bul <- 0
mid_Qian <- 0

for (i in 1:nchar(mid)){
mid_bul <- mid_bul+(aa_index[aa_index$code == substr(mid,i,i),3])
mid_Qian <- mid_Qian +(aa_index[aa_index$code == substr(mid,i,i),4])
}
mid_bul <- mid_bul/(length-2)
mid_Qian <- mid_Qian/(length-2)

# average value of all amino acids
av_bul <- 0
av_Qian <- 0

for (i in 1:nchar(input_peptide[x])){
  av_bul <- av_bul+ (aa_index[aa_index$code == substr(input_peptide[x],i,i),3])
  av_Qian <- av_Qian +(aa_index[aa_index$code == substr(input_peptide[x],i,i),4])
}

av_bul <- av_bul/(length)
av_Qian <- av_Qian/(length)

# here all the calculated descriptors are merged to single data frame
input <- c(input,av_bul,av_Qian,N1_bul,N1_Qian,mid_bul,mid_Qian,C1_bul,C1_Qian)
names <- c("length","mw","LogP","pI","V1","V2","V3","V4","V5","V6","V7","V8","1","2","3","4","5","6","7","8")
input <-t(as.data.frame(input))
colnames(input) <- names



# Print input peptide sequence        /  predict by model "fit_svm" using "input" as input-values
print(paste0(input_peptide[x],"       ",(round(predict(fit_svm,input),2))))
}
