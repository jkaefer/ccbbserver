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



asa_ss_feature_ss_prob<-function(anno, frag, flname)## the max prob among P_E, P_H, P_C
{
    info = matrix(0, ncol= 16, nrow =dim(frag)[1] )
    
    for (i in 1:dim(frag)[1]){
		fl.n   = paste("/N/slate/jkaefer/splicingdb/combined_spXout/", frag[i,1], ".spXout", sep="")
		asa.fl = paste("/N/slate/jkaefer/splicingdb/combined_spXout/", frag[i,1], ".spXout.norm", sep="")
		dis.fl = paste("/N/slate/jkaefer/splicingdb/spineD_disorderScore/", frag[i,1], ".spd", sep="")

		info[i,1]  = as.vector(unlist(anno[i,1]))
		if(file.exists(fl.n)){
			if(!file.exists(asa.fl)){
				cmd = paste('python ./regsnps-splicing/source/code_package/ASASS/convert_rASA.py ', fl.n,' > ', asa.fl , sep='')
				try(system(cmd))
			}

			spx.data = read.table(fl.n, sep="", header=F, skip=1,quote="")
			spd.data = read.table(dis.fl, sep="", header=F,quote="")
			asa.data = read.table(asa.fl, sep="", header=F,quote="")
			
			colnames(spx.data)=c("index" , "AA" , "SS"  , "phi1"   , "psi1"  ,   "P_E"    ,   "P_C"   ,    "P_H"  ,   "phi0" , "psi0" , "ASA"  ,  "S_pk" , "S_SS"  , "pk_phi" , "pk_psi" ,"pkc_phi"  , "pkc_ps")
			ed = as.numeric(frag[i,3])
			if(ed > dim(spx.data)[1]){
				ed = dim(spx.data)[1]
			}
			ss  = spx.data[seq(as.numeric(frag[i,2]), ed,1), 6:8]
			dis = spd.data[seq(as.numeric(frag[i,2]), ed,1), ]
			#ss = ss*(1-dis[,3])
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
			
			info[i,14]  = mean(asa.data[seq(as.numeric(frag[i,2]), ed,1), 3])  # ASA
			info[i,15]  = min( asa.data[seq(as.numeric(frag[i,2]), ed,1), 3])  # ASA
			info[i,16]  = max( asa.data[seq(as.numeric(frag[i,2]), ed,1), 3])  # ASA
					
		}
	
    }
    #info[,2:16] = round(info[,2:16], 4)
    write.table(file=flname, info, sep='\t', quote=F, col.names=F, row.names=F)
    
}
message('test') #4-2-2019
desp    = commandArgs(TRUE)[1] 
job_loc = commandArgs(TRUE)[2]

neu.anno = read.table(paste(job_loc, "/", desp,".snp.ref.mut.seq.disorder-matrix", sep=""), header=F, sep='\t', quote="")
neu.pro.frag = matrix(as.vector(unlist(strsplit(as.vector(unlist(neu.anno[,2])), ":"))), ncol=3, byrow=T)
neu.info = matrix(0, ncol= 13, nrow =dim(neu.pro.frag)[1] )

asa_ss_feature_ss_prob(neu.anno,  neu.pro.frag,  paste(job_loc, "/", desp, ".asa.ss.max_prob.data", sep=""))

message('test2') #4-2-2019


rm(list = setdiff(ls(), lsf.str()))

