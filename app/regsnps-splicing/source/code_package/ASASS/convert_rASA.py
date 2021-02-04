#!/usr/bin/env python
import sys

rnam1_std0 = "ACDEFGHIKLMNPQRSTVWYX"
rnam1_std = rnam1_std0[:20]
rASA_std = (115, 135, 150, 190, 210, 75, 195, 175, 200,
		170, 185, 160, 145, 180, 225, 115, 140, 155, 255, 230)
def aaDefine(rn, itype=0, DEBUG=1):
	if itype == 0: itype = len(rn)
	if itype == 3:
		ir = rnam3_std_dict.get(rn, -1)
	elif itype == 1:
		if rn == 'Z': rn = 'E'
		elif rn == 'B': rn = 'D'
		ir = rnam1_std.find(rn)
	else: die("not known type: %d" % itype)
	if DEBUG and ir < 0:
		print >>sys.stderr, 'unknown res:', rn
	return ir
def rdphipsi(fn):
	seq1, ssec1, phipsi1, Fphipsi = [], [], [], []
	for line in file(fn):
		if line[0] in ">#": continue
		ss = line.split()
		if ss[1] == 'X': continue
		rid = aaDefine(ss[1])
		print ss[0], ss[1], '%.1f' % (float(ss[10])/rASA_std[rid]*100.)
		continue

		seq1.append(aa1_std[rid]); ssec1.append(ss[2])
		asa = min(100, int(round(float(ss[10])/rASA_std[rid]*100.)))
		phi = int(round(float(ss[3])))
		psi = int(round(float(ss[4])))
		iss = int(float(ss[12])*10. - 2)
		iphi = int(abs(float(ss[15])*10. - 5) * 2)
		ipsi = int(abs(float(ss[16])*10. - 5) * 2)
		iss = max(0, min(7, iss))
		iphi = max(0, min(9, iphi))
		ipsi = max(0, min(9, ipsi))
		phipsi1.append( (phi, psi, asa) )
		Fphipsi.append( (iss, iphi, ipsi) )
	return seq1, ssec1, phipsi1, Fphipsi
if len(sys.argv) < 2:
	print >>stderr, 'usage: RUN *.phipsi'
rdphipsi(sys.argv[1])
