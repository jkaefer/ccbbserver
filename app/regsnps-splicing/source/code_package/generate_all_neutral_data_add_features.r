##parameters
desp	= commandArgs(TRUE)[1]
job_loc = commandArgs(TRUE)[2]
print('test') #4-2-2019
class_lab = 0;
##add RNA 2nd structure and structure distance features
#######-------------------------------------------- main ----------------------------------------------------
f_proximity 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.proximity", sep="")
write.csv(f_proximity, file = "my_data.csv");
f_junc 		= paste(job_loc, "/", desp, ".snp.ref.mut.seq.junction-score", sep="")
f_junc_m 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.junction-score.mut", sep="")
f_evo 		= paste(job_loc, "/", desp, ".snp.ref.mut.seq.mutBed.phylop", sep="")
f_sfrs1 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.SRSF1.score", sep="")
f_sfrs2 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.SRSF2.score", sep="")
f_sfrs5 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.SRSF5.score", sep="")
f_sfrs6 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.SRSF6.score", sep="")
f_len 		= paste(job_loc, "/", desp, ".snp.ref.mut.seq.exon-intron.len", sep="")
f_postp 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.changeBinding.pvalue", sep="")  ## added
f_mag 		= paste(job_loc, "/", desp, ".snp.ref.mut.seq.logor", sep="")
f_clus_score	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.motif_cluster.score", sep="")
f_ptm 		= paste(job_loc, "/", desp, ".snp.ref.mut.seq.ptm-matrix.dbPTM", sep="")
f_disorder 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.disorder-matrix", sep="")
f_pfam 	= paste(job_loc, "/", desp, ".snp.ref.mut.seq.pfam-matrix", sep="")
f_ef		= paste(job_loc, "/", desp, ".snp.ref.mut.seq.efscore.value", sep="")
f_2nd		= paste(job_loc, "/", desp, ".combined.Struct.compare.score", sep="")

print('test1') #4-2-2019
n.proximity 		= read.table(f_proximity, header=F, sep='\t', quote="")
n.junction.score 	= read.table(f_junc,  header=F, sep='\t', quote="")
n.junction.score.mut	= read.table(f_junc_m,header=F, sep='\t', quote="")
n.evo 			= read.table(f_evo,   header=F, sep='\t', quote="")
n.sfrs1 		= read.table(f_sfrs1, header=F, sep='\t', quote="")
n.sfrs2 		= read.table(f_sfrs2, header=F, sep='\t', quote="")
n.sfrs5 		= read.table(f_sfrs5, header=F, sep='\t', quote="")
n.sfrs6 		= read.table(f_sfrs6, header=F, sep='\t', quote="")
n.exon.intron.len 	= read.table(f_len,   header=F, sep='\t', quote="")
n.efscore		= read.table(f_ef,    header=F, sep='\t', quote="")
n.struct		= read.table(f_2nd,	  header=F, sep='\t', quote="")
n.postp 		= read.table(f_postp, header=T, sep='\t', quote="")
n.magnitude 		= read.table(f_mag,   header=T, sep='\t', quote="")
data 			= read.table(f_clus_score, header=T,  sep='\t', quote="")
n.cluster.score 	= rowSums(data[, 445:885] - data[, 4:444])
n.cluster.score.ref 	= rowSums(data[, 445:885])
n.ptm 			= read.table(f_ptm, header=T, sep='\t', quote="")
n.disorder.score 	= read.table(f_disorder, header=F, sep='\t', quote="")
n.pfam 			= read.table(f_pfam, sep='\t', header=T, quote="")

print('test2') #4-2-2019
n.pfam[,395] 	= sapply(n.pfam[,395], f<-function(x)gsub( "[a-zA-Z\\_0-9\\.\\-]+:", "",x))
n.pfam		= n.pfam[!duplicated(n.pfam[,1]),]
n.ptm  		= n.ptm[!duplicated(n.ptm[,1]),]
n.disorder.score= n.disorder.score[!duplicated(n.disorder.score[,1]),]

n.snp.id = n.disorder.score[,1]
n.ptm 	 = n.ptm[n.ptm[,1] %in% n.snp.id, ]

n.pfam 	= n.pfam[n.pfam[,1] %in% n.snp.id,]
n.snp.id= sapply(n.snp.id, f<-function(x)gsub(':[+-]$', '',x,perl=T)  )

neu.asa = read.table(paste(job_loc, "/",  desp, ".asa.ss.max_prob.data",sep=""), sep='\t', quote="", header=F)
neu.asa = neu.asa[neu.asa[,1] %in% n.snp.id, ]

n.idx 	     = match(n.snp.id, n.proximity[,3])
n.asa.idx    = match(n.snp.id, neu.asa[,1])

n.proximity 		= n.proximity[n.idx, ]
n.junction.score 	= n.junction.score[n.idx, ]
n.junction.score.mut	= n.junction.score.mut[n.idx, ]
n.evo 			= as.vector(unlist(n.evo))
n.evo 			= n.evo[n.idx]
n.sfrs1 		= n.sfrs1[n.idx,]
n.sfrs2 		= n.sfrs2[n.idx,]
n.sfrs5 		= n.sfrs5[n.idx,]
n.sfrs6 		= n.sfrs6[n.idx,]
n.exon.intron.len 	= n.exon.intron.len[n.idx,]
n.postp 		= n.postp[n.idx,]
n.magnitude 		= n.magnitude[n.idx,]
n.cluster.score 	= n.cluster.score[n.idx]
n.cluster.score.ref 	= n.cluster.score.ref[n.idx]
neu.asa			= neu.asa[n.asa.idx,]
n.efscore		= n.efscore[n.idx,]
n.struct		= n.struct[n.idx ,]

n.magnitude  = sapply(n.magnitude, round, digits =4)
n.postp      = sapply(n.postp,     round, digits =4)

data.neu  = cbind(n.exon.intron.len[,4:6], n.junction.score[,4:5], n.junction.score.mut[,4:5], n.sfrs1[,4:5],n.sfrs2[,4:5],n.sfrs5[,4:5],n.sfrs6[,4:5], n.proximity[,3:6], n.evo, n.magnitude, n.postp, n.cluster.score, n.cluster.score.ref, n.pfam[,396], rowSums(n.ptm[,3:dim(n.ptm)[2]]), n.disorder.score[,3:14])

nonempty.idx = which(rowSums(neu.asa[,2:dim(neu.asa)[2]])!=0)
data.neu  = cbind(data.neu, neu.asa, n.efscore[,3], n.struct)
data.neu  = data.neu[nonempty.idx,]
data.neu  = data.neu[, -(dim(data.neu)[2]-dim(neu.asa)[2] -1)]

write.table(file=paste(job_loc, "/", desp,".data",sep=""), data.neu, sep=',', quote=F, col.names=F, row.names=F)

#-------------------------------------------------------------------------

##remove SNPs on splicing site:
n.ss_idx = which(data.neu[,17] == 0)
n.ss_idx = c(n.ss_idx, which(data.neu[,18] <= 2))
print('got_here')
data.neu.ss = data.neu[n.ss_idx, ]
if(length(n.ss_idx) >0){
	print('in writing loop')
	#write.table(file=paste(job_loc, "/", desp,".weka.data.ss.info",sep=""), data.neu.ss, sep=',', quote=F, col.names=F, row.names=F)
	#data.neu.ss = data.neu.ss[,-c(16,17,18,19)]## still keep proximity for SNPs on splicing site, but not used for predict/train
	data.neu.ss[,6] = data.neu.ss[,6] - data.neu.ss[,4]
	data.neu.ss[,7] = data.neu.ss[,7] - data.neu.ss[,5]
	for (i in 1:dim(data.neu.ss)[1]){
		if(data.neu.ss[i,6] == 0 & data.neu.ss[i,7] !=0 ){
			data.neu.ss[i,6] = data.neu.ss[i,7]
		}
	}
	data.neu.ss = data.neu.ss[,-7]
	write.table(file=paste(job_loc, "/", desp,".weka.data.ss.info",sep=""), data.neu.ss, sep=',', quote=F, col.names=F, row.names=F)
	data.neu.ss = data.neu.ss[,-c(15,16,17,18)]
}

data.class   = rep(class_lab, dim(data.neu.ss)[1])
out.data.neu = cbind(data.neu.ss, data.class)
write.table(file=paste(job_loc, "/", desp,".weka.data.ss.info",sep=""), out.data.neu, sep=',', quote=F, col.names=F, row.names=F)#remove proximity 
wekafl = paste(job_loc, "/", desp,".weka.data.ss.pred.arff", sep="")
arff.head = read.table('./regsnps-splicing/source/code_package/model/weka.head.ss', header=F, quote="", sep='\n')
write.table(file= wekafl, arff.head, quote=F, row.names=F, col.names=F)
write.table(file= wekafl, out.data.neu, sep=',', quote=F, row.names=F, col.names=F, append=T)


## filter selected features, for SNPs on exon body
if(length(n.ss_idx) >0){
    data.neu.inner = data.neu[-n.ss_idx, ]
}else{
    data.neu.inner = data.neu
}

if(length(n.ss_idx) != dim(data.neu)[1]){
	write.table(file=paste(job_loc, "/", desp,".weka.data.inner.info",sep=""), data.neu.inner, sep=',', quote=F, col.names=F, row.names=F)
	data.neu.inner = data.neu.inner[,-c(6,7,16,19)]
}

data.class   = rep(class_lab, dim(data.neu.inner)[1])
out.data.neu = cbind(data.neu.inner, data.class)
#write.table(file=paste(job_loc, "/", desp,".weka.data.inner.info",sep=""), data.neu.inner, sep=',', quote=F, col.names=F, row.names=F)

wekafl = paste(job_loc, "/", desp, ".weka.data.inner.pred.arff", sep="")
arff.head = read.table('./regsnps-splicing/source/code_package/model/weka.head.inner', header=F, quote="", sep='\n')
write.table(file= wekafl, arff.head, quote=F, row.names=F, col.names=F)
write.table(file= wekafl, out.data.neu, sep=',', quote=F, row.names=F, col.names=F, append=T)

rm(list = setdiff(ls(), lsf.str()))
























