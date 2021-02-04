#load('binding.var.RData') # binding score mean and variance
load('/N/slate/jkaefer/splicingdb/binding.var.all.RData') # binding score mean and variance

## Information Content, for prior of binding probability
calcu_ic<-function(pwm){
    #zero_cols = colSums(pwm==0);#how many zero for each column
    #idx = which(zero_cols>0);
    #if(length(idx)>0){
    #    colsum = colSums(pwm);
    #    add_col = sqrt(colsum[idx]);
    #    for(i in 1:length(add_col)){
    #        pwm[,idx[i]] = pwm[,idx[i]]+add_col[i]/4*matrix(data=1,ncol=1,nrow=4);
    #    }
    #}
    #pwm = apply(pwm,2,f<-function(vec){vec/sum(vec)})# convert count to frequency
    #bg = c(0.3, 0.2, 0.2, 0.3)
	bg = c(.25, .25, .25, .25);
	#t = log2(pwm/0.25);
    #t = pwm*t;#p*log2(p/0.25)
	fq = (2^pwm)*bg
	t = fq*pwm;
    t =colSums(t);
    ic = sum(t)
    return(ic);
}

## Determine beta distribution parameters based on mode and a quantile value qtl_val
## mode: (a-1)/(a+b-2)  max: max value to search for a
prec_beta_prior_params<-function(mode,quantile,qtl_val,max){# 1/2^ic, 0.05, mode/10, max=10
    f<-function(x,mode,quantile,qtl_val){
        pbeta(qtl_val,x,((x-1)/(mode)-(x-1)+1))-quantile; 
    }
    r = uniroot(f,c(1.1, max),tol=1e-4,mode=mode,quantile=quantile,qtl_val=qtl_val); # find root which makes f = 0
    a = r$root; b=(a-1)/mode-(a-1)+1;#relationship between alpha, beta and mode
    return(c(a,b));
}

## posterior probability of binding change when the prior follows a beta distribution specified by a and b
betaprior_cbpostprob_bayes<-function(score,bmean,bsd,nbmean,nbsd,a,b)
{
	if(max(score) < (bmean-3*bsd)){
		return(0);
	}else{
		ref_score <- score[1];
		alt_score <- score[2];
		b_ref <- pnorm(ref_score,bmean,bsd);
		b_alt <- pnorm(alt_score,bmean,bsd);
		nb_ref <- 1-pnorm(ref_score,nbmean,nbsd);
		nb_alt <- 1-pnorm(alt_score,nbmean,nbsd);

		f<-function(x){ ## Beta distribution, mode=IC based pb, 2.5% is pb/10
			## P(A=B|SA)*P(R=NB|SR)+P(A=NB|SA)*P(R=B|SR)
			## P(B)*(1-P(B))(P(SA|A=B)P(SR|RNB)+P(SA|A=NB)P(SR|R=B))/(P(SA)*P(SR))
			#options(digits=20);
			x*(1-x)*(b_ref*nb_alt+nb_ref*b_alt)/((x*b_ref+(1-x)*nb_ref)*(x*b_alt+(1-x)*nb_alt))*dbeta(x,shape1=a,shape2=b);
		}
		#cat("Beta, center at 0.01, alpha=11,beta=991\n")
		post=integrate(f,0,1,rel.tol=1.e-10,stop.on.error=FALSE);
		return(post$value);
	}
}

## posterior probability of binding given observed binding score when the prior follows a beta distribution specified by a and b
betaprior_postprob_bayes<-function(score,bmean,bsd,nbmean,nbsd,a,b){
    conb_s <- pnorm(score,bmean,bsd);
    connb_s <- 1-pnorm(score,nbmean,nbsd);
    
    f<-function(x){ ## Beta distribution, mode=IC based pb, 2.5% is pb/10
    ## P(B|S)=P(S|B)*P(B)/(P(S|B)*P(B)+P(S|NB)*P(NB));
    ## P(B|S)=P(S|B)*P(B)/((P(S|B)-P(S|NB))*P(B)+P(S|NB));
    #options(digits=20);
    conb_s*x/((conb_s-connb_s)*x+connb_s)*dbeta(x,shape1=a,shape2=b);
    }
    #cat("Beta, center at 0.01, alpha=11,beta=991\n")
    post=integrate(f,0,1,rel.tol=1.e-10,stop.on.error=FALSE);
    return(post$value);
}
betaprior_postprob_bayes_vec<-function(score_vec,bmean,bsd,nbmean,nbsd,a,b){
    s=vector();
    s[1]=betaprior_postprob_bayes(score_vec[1],bmean,bsd,nbmean,nbsd,a,b);
    s[2]=betaprior_postprob_bayes(score_vec[2],bmean,bsd,nbmean,nbsd,a,b);
    return(s);
}

#ratio.filename	= "/home/zhangxin/Project/indels/motif/pwm_f.53.list";
#ratio.filename = "pwm_f.53.list";
ratio.filename = "./code_package/db/pwm.all.list";
pfm.ratio.fl    = read.table(ratio.filename,strip.white = TRUE)
pfm.ratio.fl    = as.vector(unlist(pfm.ratio.fl))

#-------------
ic = rep(0, length(pfm.ratio.fl))

for (i in 1:length(pfm.ratio.fl)) # for each pssm
{
    pfm.ratio = as.matrix(read.table(pfm.ratio.fl[i]))  # ratio file 
#	pfm.ratio= (2^pfm.ratio)*0.25
    ic[i] = calcu_ic(pfm.ratio)                              # 1/2^ic is the frequency that any sequence of a RBP will be seen in genome     
}



#tg.bind.score    = read.table("/data/pro/indels/1000Genome/motif/tg.indel.intronexon.anno.mut.cor.bed.seq.motif_change.out", quote="", sep="\t", header=T)
#tg.bind.score = read.table("test.tg.indel.vcf.exonic.intronexon.anno.mut.cor.bed.seq.motif_change.out",quote="", sep="\t", header=T)
fprefix = commandArgs(TRUE)[1]
#fprefix = "hgmd.vcf";

sfl = paste(fprefix, ".motif_change.out", sep='' )
tg.bind.score = read.table(sfl,quote="", sep="\t", header=F, skip=1)
#pwm.list = as.vector(unlist(read.table("/data/pro/indels/motif/pwm")))
pwm.list = as.vector(unlist(read.table("./code_package/db/pwm.all")))

#tg.bind.score = tg.bind.score[1:(56*10),]

#logOR and pvalue for 1000Genome, change 53 -> 57
n_pwm = 201
indel.list.tg = "0"
indel.logor.tg.yy   = matrix(0, dim(tg.bind.score)[1]/n_pwm, n_pwm)
indel.invert.binding.pvalue.tg = matrix(0, dim(tg.bind.score)[1]/n_pwm, n_pwm)

for (i in 1:(dim(tg.bind.score)[1]/n_pwm))#each indel
{

    block.data.tg = as.matrix(tg.bind.score[(n_pwm*(i-1)+1) : (n_pwm*i),])
    indel.list.tg = c(indel.list.tg, block.data.tg[1,4])
    k=0
    
    for (j in 1:(dim(block.data.tg)[1]))#each motif for one indel
    {
     #  motif.id = which(pwm.list == block.data.tg[j,1])
	pn = block.data.tg[j,1]
	#pn = sub('/data/pro/as/data/pwms_all_motifs/', '', pn) # 11/24/2014
	#motif.id = which(grep(block.data.tg[j,1], pwm.list)) # modified 10/12/2013
	motif.id = grep(pn, pwm.list)

	if(length(motif.id) != 0 )
	#if(j != 6 & j!= 11 & j != 55)
	{
	    sv.binding.mean     = binding.mean[motif.id]
        sv.binding.var 	= binding.var[motif.id]
        sv.nonbinding.mean  = nonbinding.mean[motif.id]
        sv.nonbinding.var   = nonbinding.var[motif.id]
          
        k = k+1
	    sv.score.mut = as.numeric(block.data.tg[j, 8])
	    sv.score.ref = as.numeric(block.data.tg[j, 10])

	   #larger.score = pmax(sv.score.ref, sv.score.mut)
	   #pro.larger.score.binding = pnorm(larger.score, mean = sv.binding.mean, sd = sv.binding.var, lower.tail = TRUE, log.p = FALSE)
	        
            
            #HYY's method
	    pro.ref.binding         = pnorm( sv.score.ref, mean = sv.binding.mean, sd = sv.binding.var, lower.tail = TRUE, log.p = FALSE)
	    pro.ref.nonbinding.yy   = pnorm( sv.score.ref, mean = sv.nonbinding.mean, sd = sv.nonbinding.var, lower.tail = FALSE, log.p = FALSE)
	    pro.mut.binding         = pnorm( sv.score.mut, mean = sv.binding.mean, sd = sv.binding.var, lower.tail = TRUE, log.p = FALSE)
	    pro.mut.nonbinding.yy   = pnorm( sv.score.mut, mean = sv.nonbinding.mean, sd = sv.nonbinding.var, lower.tail = FALSE, log.p = FALSE)
	    odd.ref.yy = (pro.ref.binding) / pro.ref.nonbinding.yy
	    odd.mut.yy = (pro.mut.binding) / pro.mut.nonbinding.yy
	    log.odd.ratio.indel = log2(odd.mut.yy / odd.ref.yy)
	        
        #indel.logor.tg.yy[i,j] =  log.odd.ratio.indel*pro.larger.score.binding
        indel.logor.tg.yy[i,j] =  log.odd.ratio.indel
            
        #Pvalue
        mode = 1/2^ic[motif.id]
        quantile = 0.005
        qtl_val = mode/10
        max=10
        beta_par = prec_beta_prior_params(mode, quantile, qtl_val, max)
        indel.invert.binding.pvalue.tg[i,j] = betaprior_cbpostprob_bayes(c(sv.score.ref ,sv.score.mut), sv.binding.mean , sv.binding.var , sv.nonbinding.mean , sv.nonbinding.var , beta_par[1] , beta_par[2])
		#s=vector();
		#s[1] = betaprior_postprob_bayes(sv.score.mut, sv.binding.mean, sv.binding.var, sv.nonbinding.mean, sv.nonbinding.var, beta_par[1],beta_par[2])
		#s[2] = betaprior_postprob_bayes(sv.score.ref, sv.binding.mean, sv.binding.var, sv.nonbinding.mean, sv.nonbinding.var, beta_par[1],beta_par[2])
		#indel.invert.binding.pvalue.tg[i,j] = s[1]-s[2];
    }
  }
}


indel.list.tg = indel.list.tg[-1]

#pro = as.matrix(read.table("RBPDB_v1.2.2_proteins_2011-07-19.tdt", quote="", sep="\t", header=F))
#exp = read.table("RBPDB_v1.2.2_experiments_2011-07-19.tdt", quote="", sep="\t", header=F)
#protexp = read.table("RBPDB_v1.2.2_protExp_2011-07-19.tdt", quote="", sep="\t", header=F)
#err.pwm = c(6,11,55)

anno = matrix(0, n_pwm, 3 )
aa = block.data.tg[,1]
#aa = aa[-err.pwm]
aa = matrix(unlist(strsplit(aa,'_')), nc = 2, byrow=T)


for (i in 1:dim(anno)[1])
{	
	expid = aa[i, 1]
#	expid = sub('PWMDir/', '', expid)
#	if('SRSF' !=expid){
#		anno[i,1] = expid
#		pro_id = protexp[protexp[,2] == expid,1]
#		anno[i,2] = pro[pro[,1] == pro_id, 5]
#		anno[i,3] = pro[pro[,1] == pro_id, 6]
#	}else{
		anno[i,2] = expid;
		anno[i,1] = '';
		anno[i,3] = '';
#	}
}

#indel.logor.tg.yy = indel.logor.tg.yy[,-err.pwm]
#indel.invert.binding.pvalue.tg = indel.invert.binding.pvalue.tg[, -err.pwm]

rownames(indel.logor.tg.yy) = indel.list.tg
#colnames(indel.logor.tg.yy) = block.data.tg[,1]
colnames(indel.logor.tg.yy) = anno[,2]
rownames(indel.invert.binding.pvalue.tg) = indel.list.tg
#colnames(indel.invert.binding.pvalue.tg) = block.data.tg[,1]
colnames(indel.invert.binding.pvalue.tg) = anno[,2]

sfl = paste(fprefix, ".logor", sep='' )
write.table(file = sfl,indel.logor.tg.yy , quote=F, row.names =F , sep = "\t")

sfl = paste(fprefix, ".changeBinding.pvalue", sep='' )
write.table(file = sfl,indel.invert.binding.pvalue.tg , quote=F, row.names =F , sep = "\t")

#save.image(file='maxOR.r.image')

