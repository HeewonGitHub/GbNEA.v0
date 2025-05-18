library(HDeconometrics)

# Input 1 
# Gene network of Phenotype C: 
NW_C<-read.table("NW_C.csv",sep=",")
# Input 2 : Gene network of Phenotype N
NW_N<-read.table("NW_N.csv",sep=",")
# Input 3 : Genes in a specific pathway
PW_genes<-read.table("PW_genes.csv",sep=",")[,1]


##### 1. Compute Regulatory effect 
RE_C<-matrix(numeric(length(unique(NW_C[,1]))*2),ncol=2)
rownames(RE_C)<-unique(NW_C[,1])
RE_N<-matrix(numeric(length(unique(NW_N[,1]))*2),ncol=2)
rownames(RE_N)<-unique(NW_N[,1])

for (r in 1:length(unique(NW_C[,1]))){
if (nrow(RE_C)>=r){
RE_C[unique(NW_C[,1])[r],1]<-sum(as.numeric(NW_C[NW_C[,1]==unique(NW_C[,1])[r],"RE"]))
}
if (nrow(RE_N)>=r){
RE_N[unique(NW_N[,1])[r],1]<-sum(as.numeric(NW_N[NW_N[,1]==unique(NW_N[,1])[r],"RE"]))
}
}
C_RG<-unique(c(rownames(RE_C),rownames(RE_N)))
RGmtx<-matrix(numeric(length(C_RG)*3),ncol=3)
rownames(RGmtx)<-C_RG
RGmtx[rownames(RE_C),1]<-RE_C[,1]
RGmtx[rownames(RE_N),2]<-RE_N[,1]
colnames(RGmtx)<-c("C","N","Dif")
RGmtx[,3]<-RGmtx[,1]-RGmtx[,2]

####### 2. Compute Enrichment score based on regulatory effect
SC_RG<-matrix(numeric(length(RGmtx[,3])*6),ncol=6)
rownames(SC_RG)<-rownames(RGmtx)
colnames(SC_RG)<-c("R","IN","NiN","Hit","Miss","SC")
SC_RG[,1]<-RGmtx[,3]
SC_RG<-SC_RG[order(SC_RG[,1],decreasing=TRUE),]
SC_RG[is.na(match(rownames(SC_RG),PW_genes))==FALSE,2]<-1
SC_RG[is.na(match(rownames(SC_RG),PW_genes))==TRUE,3]<-1
GST<-seq(0,0,length=length(SC_RG))
GST[SC_RG[,2]==1]<-"2"
N<-nrow(SC_RG)
Nh<-sum(SC_RG[,2])
NR<-sum(SC_RG[,2]) ## p==0 case
if (NR!=0){
SC_RG[,4]<-cumsum(SC_RG[,2])/NR
}
SC_RG[,5]<-cumsum(SC_RG[,3])/(N-Nh)
SC_RG[,6]<-SC_RG[,4]-SC_RG[,5]
### Enrichment score based on Regulatory effect
ES_RG<-SC_RG[,6][order(abs(SC_RG[,6]),decreasing=TRUE)[1]]

####### 3. Compute Jaccard Distance 
RG_T<-unique(c(NW_C[,1],NW_N[,1]))
JI<-matrix(numeric(length(RG_T)*3),ncol=3)
rownames(JI)<-RG_T
for (r in 1:nrow(JI)){
JI[r,1]<-1-
length(intersect(NW_C[NW_C[,1]==RG_T[r],2],NW_N[NW_N[,1]==RG_T[r],2]))/
length(unique(c(NW_C[NW_C[,1]==RG_T[r],2],NW_N[NW_N[,1]==RG_T[r],2])))
}

####### 4. Compute Enrichment Score based on Regulatory effect and Jaccard distance
SC_RGJI<-matrix(numeric(length(RGmtx[,3])*6),ncol=6)
rownames(SC_RGJI)<-rownames(RGmtx)
colnames(SC_RGJI)<-c("R","IN","NiN","Hit","Miss","SC")
SC_RGJI[rownames(RGmtx),1]<-JI[rownames(RGmtx),1]*RGmtx[,3]
SC_RGJI<-SC_RGJI[order(SC_RGJI[,1],decreasing=TRUE),]
SC_RGJI[is.na(match(rownames(SC_RGJI),PW_genes))==FALSE,2]<-1
SC_RGJI[is.na(match(rownames(SC_RGJI),PW_genes))==TRUE,3]<-1
GST<-seq(0,0,length=length(SC_RGJI))
GST[SC_RGJI[,2]==1]<-"2"

N<-nrow(SC_RGJI)
Nh<-sum(SC_RGJI[,2])
NR<-sum(SC_RGJI[,2]) ## p==0 case
if (NR!=0){
SC_RGJI[,4]<-cumsum(SC_RGJI[,2])/NR
}
SC_RGJI[,5]<-cumsum(SC_RGJI[,3])/(N-Nh)
SC_RGJI[,6]<-SC_RGJI[,4]-SC_RGJI[,5]
### Enrichment score based on Regulatory effect and Jaccard Distance
ES_RGJD<-SC_RGJI[,6][order(abs(SC_RGJI[,6]),decreasing=TRUE)[1]]


# Output 1: Enrichment score based on Regulatory effect
ES_RG
# Output 2: Enrichment score based on Regulatory effect and Jaccard Distance
ES_RGJD
