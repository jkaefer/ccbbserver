fprefix = commandArgs(TRUE)[1]

#sfl = paste(fprefix, "weka.tg.data", sep='.' )
#tg.data = read.table(sfl, quote="", sep='\t', header=F)

sfl = paste(fprefix, "changeBinding.pvalue", sep='.' )
tg.p   = read.table(sfl, header=T)

sfl = paste(fprefix, "logor", sep='.' )
tg.or   = read.table(sfl, header=T)

y = which(tg.p   >= 0.5, arr.in =T)## use 0.5, Jul, 28,2013
count.change.binding.tg   = rowSums(tg.p[,   1:53] > 0.5)## use 0.5, Jul, 28,2013
affected.rbp = ""
for (i in 1:dim(tg.p)[1])
{
	affected.rbp = c(affected.rbp, paste(colnames(tg.p)[which(tg.p[i, 1:53]>= 0.5)], collapse = ";"))
}
affected.rbp = affected.rbp[-1]
max.or.tg   = apply(abs(tg.or), 1, max)
max.or.tg = round(max.or.tg,2)
tg.data =   cbind( max.or.tg,   count.change.binding.tg)

sfl = paste(fprefix, "magnitude.change.binding", sep='.' )
write.table(tg.data, file = sfl, quote=F, sep='\t', col.names=F, row.names=F )

#sfl = paste(fprefix, "affected.RBP", sep='.' )
#write.table(affected.rbp, file = sfl, quote=F, sep=',', col.names=F, row.names=F )
