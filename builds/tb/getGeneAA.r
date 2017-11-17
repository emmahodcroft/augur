#In order to create an reference/outgroup .gb file with translations,
#we need to find the genes and then add the translations to the .gb file!

setwd("C:/Users/Emma/bash/githubE/augur/builds/tb/metadata")

#Now that code has run, this file is 'tb_outgroup_orig.gb'
tbR <- read.delim("tb_outgroup.gb", header=F, skip=31, quote="")
tbR <- tbR$V1
tbR <- as.character(tbR)

##
#This part goes through and finds all of the CDS sections which have both
#a gene name AND a protein.id - because you can get protein.id WITHOUT it 
#being a gene - at the moment these are excluded - too complicated

#This is a bit bulky but the challenge is identifying the CDS 'section' when
#the next thing might be 'gene' or 'tRNA' or whatever. And then only looking
#within this section for the 'gene' and 'protein_id' - or else it could pull
#the gene from one CDS and the protein_id from the next CDS!

CDS <- grep("    CDS",tbR)

genes <- c()
proteins <- c()
unmatchedGenes <- c()	#to check
unmatchedProteins <- c()	#to check	
emptyCDS <- c()		#to check

for(i in 1:length(CDS)){
	curCDS <- CDS[i]
	
	if(i == length(CDS)){	#if at last CDS in the file
		nextCDS <- length(tbR)
		endCDS <- length(tbR)
	} else {
		nextCDS <- CDS[i+1]
		#this finds the 'start' of the next gene/CDS/tRNA/whatever
		#don't search for genes/proteins beyond here
		endCDS <- grep("^     [[:alpha:]]", tbR[curCDS+1:nextCDS])[1] + curCDS
	}
	
	gen <- tbR[curCDS:endCDS][grep("/gene=",tbR[curCDS:endCDS])]
	prot <- tbR[curCDS:endCDS][grep("/protein_id=",tbR[curCDS:endCDS])]
	
	if(length(prot)==0 | length(gen)==0){
		if(length(prot)==0 & length(gen)!=0 ){
			unmatchedGenes <- c(unmatchedGenes, gen)
		}
		if(length(gen)==0 & length(prot)!=0 ){
			unmatchedProteins <- c(unmatchedProteins, prot)
		}
		if( length(gen)==0 & length(prot)==0){
			emptyCDS <- c(emptyCDS, curCDS)
		}
	} else  {
		genes <- c(genes, gen)
		proteins <- c(proteins, prot)
	}
	
}

#for checking
length(genes)
length(proteins)
length(unmatchedGenes)
length(unmatchedProteins)
length(emptyCDS)

#for the current TB reference these should be the numbers
# > length(genes)
# [1] 1673
# > length(proteins)
# [1] 1673
# > length(unmatchedGenes)
# [1] 0
# > length(unmatchedProteins)
# [1] 2667
# > length(emptyCDS)
# [1] 0


###
#Now write out the protein IDS to a text file so we can get them from 
#genbank batch entrez - I'm sure you could automate this.

protNums <- gsub("                     /protein_id=", "", proteins)
protNums <- gsub("\"","",protNums)

write.table(protNums, file="protAccess.txt", quote=F, row.names=F, col.names=F)

#######
#######
#In here, go to genbank batch entrez and put in the .txt file
#download the results to a FASTA file with name 'geneProt.fasta' inside metadata folder
#######
#######

##
#Now read in the protein sequences and add them to the original .gb file

library(seqinr)
geneP <- read.fasta("geneProt.fasta", seqtype="AA", as.string=T)

tbR <- read.delim("tb_outgroup.gb", header=F, quote="")
tbR <- tbR$V1
tbR <- as.character(tbR)
tbRBackup <- tbR

#Find where all of the proteins/genes that we found should go in the file, roughly
#(find the right CDS section)
locs <- c()
for(i in 1:length(proteins)){
	prot <- proteins[i]
	curLoc <- grep(prot, tbR)
	locs <- c(locs, curLoc)
}

#Write out each file of the .gb file back out, but if we get to the location where
#a protein_id we've found it, write out the translation instead, then continue
#this takes a while
for(i in 1:length(tbR)){
	if(i%in%locs){
		write.table(tbR[i],file="tbr.gb",quote=F,row.names=F,col.names=F,append=T)
		j <- which(locs==i)
		
		if( attr(geneP[[j]],"name")!=protNums[j] ) {
			print("Not matching what we thought!!")
			break
		}
		trans <- paste("                     /translation=\"",geneP[[j]][1],"\"",sep="",collapse="")
		write.table(trans, file="tbr.gb",quote=F,row.names=F,col.names=F,append=T)
		
	} else {
		write.table(tbR[i],file="tbr.gb",quote=F,row.names=F,col.names=F,append=T)
	}
}

