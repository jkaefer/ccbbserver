string_len<-function(stri, char)
{
	len = vector()
	idx = 0;
	record = 0
	for (i in 1:length(stri)){
		if(stri[i] == char & record == 0 ){
			record=1
			idx = idx+1
			len[idx] = 1
		}else if(stri[i] == char & record ==1){
			len[idx] = len[idx] +1
		}else if (stri[i] != char){
			record = 0
		}
	}
	if(length(len) ==0){
		len= 0
	}
	return(len)
}


asa_ss_feature<-function(anno, frag, flname)
{
	info = matrix(0, ncol= 16, nrow =dim(frag)[1] )
	
	for (i in 1:dim(frag)[1]){
	fl.n = paste("/data/pro/as/data/HGMD/exon_SNP/yunlong_exonic_splicing/spineX/combined_spXout/", frag[i,1], ".spXout", sep="")
	info[i,1]  = as.vector(unlist(anno[i,1]))
	if(file.exists(fl.n)){
		spx.data = read.table(fl.n, sep="", header=F, skip=1,quote="")
		colnames(spx.data)=c("index" , "AA" , "SS"  , "phi1"   , "psi1"  ,   "P_E"    ,   "P_C"   ,    "P_H"  ,   "phi0" , "psi0" , "ASA"  ,  "S_pk" , "S_SS"  , "pk_phi" , "pk_psi" ,"pkc_phi"  , "pkc_ps")
		ed = as.numeric(frag[i,3])
		if(ed > dim(spx.data)[1]){
			ed = dim(spx.data)[1]
		}
		ss = spx.data[seq(as.numeric(frag[i,2]), ed,1),]
		info[i,2] = length(which(ss[,3] == 'C'))/dim(ss)[1]
		info[i,3] = length(which(ss[,3] == 'H'))/dim(ss)[1]
		info[i,4] = length(which(ss[,3] == 'E'))/dim(ss)[1]
		
		info[i,5] = max(string_len(ss[,3], 'C'))
		info[i,6] = min(string_len(ss[,3], 'C'))
		info[i,7] = mean(string_len(ss[,3],'C'))
		
		info[i,8] = max( string_len(ss[,3], 'H'))
		info[i,9] = min( string_len(ss[,3], 'H'))
		info[i,10] = mean(string_len(ss[,3],'H'))	
			
		info[i,11] = max(string_len(ss[,3], 'E'))
		info[i,12] = min(string_len(ss[,3], 'E'))
		info[i,13] = mean(string_len(ss[,3],'E'))
	}
	norm_asa = paste("/data/pro/as/data/HGMD/exon_SNP/yunlong_exonic_splicing/spineX/combined_spXout/", frag[i,1], ".norm", sep="")
	if(!file.exists(norm_asa)){
		system();
	}
		n_asa = read.table(norm_asa, sep="", header=F, quote="")
		info[i,14]  = mean(n_asa[seq(as.numeric(frag[i,2]), as.numeric(frag[i,3]),1), 3]) # ASA
		info[i,15]  = min( n_asa[seq(as.numeric(frag[i,2]), as.numeric(frag[i,3]),1), 3])  # ASA
		info[i,16]  = max( n_asa[seq(as.numeric(frag[i,2]), as.numeric(frag[i,3]),1), 3])  # ASA

	
	
	
	}
	write.table(file=flname, info, sep='\t', quote=F, col.names=F, row.names=F)
	
}

asa_ss_feature_ss_prob<-function(anno, frag, flname)## the max prob among P_E, P_H, P_C
{
	info = matrix(0, ncol= 16, nrow =dim(frag)[1] )
	
	for (i in 1:dim(frag)[1]){
		fl.n = paste("/data/pro/as/data/HGMD/exon_SNP/yunlong_exonic_splicing/spineX/combined_spXout/", frag[i,1], ".spXout", sep="")
		info[i,1]  = as.vector(unlist(anno[i,1]))
		if(file.exists(fl.n)){
			spx.data = read.table(fl.n, sep="", header=F, skip=1,quote="")
			colnames(spx.data)=c("index" , "AA" , "SS"  , "phi1"   , "psi1"  ,   "P_E"    ,   "P_C"   ,    "P_H"  ,   "phi0" , "psi0" , "ASA"  ,  "S_pk" , "S_SS"  , "pk_phi" , "pk_psi" ,"pkc_phi"  , "pkc_ps")
			ed = as.numeric(frag[i,3])
			if(ed > dim(spx.data)[1]){
				ed = dim(spx.data)[1]
			}
			ss = spx.data[seq(as.numeric(frag[i,2]), ed,1), 6:8]
			info[i,2] = max(apply(ss,1, max))
			info[i,3] = min(apply(ss,1, max))
			info[i,4] = round(mean(apply(ss,1, max)),4)
			
			info[i,5]  = round(mean(ss[,1]),4) # beta sheet
			info[i,6]  = min( ss[,1])  # beta sheet
			info[i,7]  = max( ss[,1])  # beta sheet
			info[i,8]  = round(mean(ss[,2]),4) # coil
			info[i,9]  = min( ss[,2])  # coil
			info[i,10] = max( ss[,2])  # coil
			info[i,11] = round(mean(ss[,3]),4) # helix
			info[i,12] = min( ss[,3])  # helix
			info[i,13] = max( ss[,3])  # helix
			
			norm_asa = paste(fl.n, ".norm", sep="")
			if(!file.exists(norm_asa)){
				cmd = paste('python convert_rASA.py ', fl.n,' > ', norm_asa , sep='')
				try(system(cmd))
			}
			n_asa = read.table(norm_asa, sep="", header=F, quote="")
			info[i,14]  = round(mean(n_asa[seq(as.numeric(frag[i,2]), ed,1), 3]),4) # ASA
			info[i,15]  =        min(n_asa[seq(as.numeric(frag[i,2]), ed,1), 3])  # ASA
			info[i,16]  =        max(n_asa[seq(as.numeric(frag[i,2]), ed,1), 3])  # ASA

		}
	
	}
	#info[,2:16] = round(info[,2:16], 4)
	write.table(file=flname, info, sep='\t', quote=F, col.names=F, row.names=F)
	
}

desp = commandArgs(TRUE)[1] 

#neu.anno = read.table("/data/pro/as/data/HGMD/exon_SNP/yunlong_exonic_splicing/all_data/new_RBP/neutral.snp.ref.mut.seq.disorder-matrix", header=F, sep='\t', quote="")
neu.anno = read.table(paste(desp,"snp.ref.mut.seq.disorder-matrix",sep="."), header=F, sep='\t', quote="")
neu.pro.frag = matrix(as.vector(unlist(strsplit(as.vector(unlist(neu.anno[,2])), ":"))), ncol=3, byrow=T)
neu.info = matrix(0, ncol= 13, nrow =dim(neu.pro.frag)[1] )

####
#asa_ss_feature(hgmd.anno, hgmd.pro.frag, "hgmd.asa.ss.data")
#asa_ss_feature(neu.anno,  neu.pro.frag,  "neu.asa.ss.data")

#asa_ss_feature_ss_prob(hgmd.anno, hgmd.pro.frag, "hgmd.asa.ss.max_prob.data")
#asa_ss_feature_ss_prob(neu.anno,  neu.pro.frag,  "neu.asa.ss.max_prob.data")

asa_ss_feature_ss_prob(neu.anno,  neu.pro.frag,  file = paste(desp, ".asa.ss.max_prob.data", sep=""))


