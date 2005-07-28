include ${SEPINC}/SEP.top
RESDIR = ./Fig
DATDIR = ./Data
BINDIR = .

RESULTSNR = dirres
RESULTSCR = dirjbm
RESULTSER = dirmcg dirmcd diritr dircvv dirrst

UF90LIBS = ${BEILIB90}

TA2VPLOT = Ta2vplot crowd=.82 >/dev/null
INVRAS = y

default: view

LBL = label2=iter label1=sample label3=amp 
${RESDIR}/dirm%.v : ${DATDIR}/mparab.HH ${BINDIR}/Miss%.x
	<   ${DATDIR}/mparab.HH Dd esize=4 >   mparab.H
	Spike n1=101 k1=50 d1=1 o1=0 n2=1 o2=0 mag=10 ${LBL} > in.H 
	grep n1 in.H >junk.n1
	<mparab.H  Pad par=junk.n1     >junk.pad
	<in.H ${BINDIR}/Miss$*.x filtin=mparab.H  >junk.miss
	<junk.miss Wiggle wherelabel1=t labelsz=10 pclip=100 \
		title=" " > junk.V out=ja.v
	<junk.miss Thplot title=" " axissz=10 > junk.V out=junk.v
	vp_Unrotate junk.v > jb.v
	<junk.miss Reverse axis=2 >junk.H
	< junk.H Transp >junk1.H
	< junk1.H Contour title=" " nc=5 c0=0 labelsz=10 >c.V out=jc.v
	vp_SideBySideIso ja.v jb.v jc.v > ${RESDIR}/dirm$*.v

cgerror_s.H:	cmp_tpm.H CGscan10_w1.P ${BINDIR}/CGscan.x
		${BINDIR}/CGscan.x < cmp_tpm.H par=CGscan10_w1.P \
		miter=2 psun1=0 psun2=0 f_err=cgerror_s.H > cgscan_s.H

cgerror_w.H:	cmp_tpm.H CGscan10_w1.P ${BINDIR}/CGscan.x
		${BINDIR}/CGscan.x < cmp_tpm.H par=CGscan10_w1.P \
		miter=2 psun1=1 psun2=1 f_err=cgerror_w.H > cgscan_w.H

cdscan_w.H cderror_w.H:	cmp_tpm.H CGscan10_w1.P ${BINDIR}/CGscan.x
		${BINDIR}/CGscan.x < cmp_tpm.H par=CGscan10_w1.P \
		miter=10 psun1=2 psun2=0 f_err=cderror_w.H > cdscan_w.H

cdcmp_w.H:	cdscan_w.H CGscan10_w1.P ${BINDIR}/CGscan.x
		${BINDIR}/CGscan.x < cdscan_w.H par=CGscan10_w1.P \
		niter=0 adj=1 psun2=0			> cdcmp_w.H

cgiter.v:	cgerror_s.H cgerror_w.H cgiter.p
		Graph < cgerror_s.H min2=0 max2=0.4 out=junk0.v \
			labelsz=10 titlesz=12 \
			title="Unweighted versus Pseudounitary" > /dev/null
		Graph < cgerror_w.H min2=0 max2=0.4 fastplot=20 \
			out=junk1.v wantaxis=n wanttitle=n > /dev/null 
		vp_Overlay junk0.v junk1.v > cgjunk.v
		vp_annotate <cgjunk.v text=cgiter.p batch=y >cgiter.v

cditer.v:	cgerror_w.H cderror_w.H cditer.p
		Graph < cgerror_w.H min2=0 max2=0.2 out=junk0.v \
			labelsz=10 titlesz=12 \
			title="Gradients versus Directions" > /dev/null
		Graph < cderror_w.H min2=0 max2=0.2 fastplot=20 \
			out=junk1.v wantaxis=n wanttitle=n > /dev/null
		vp_Overlay junk0.v junk1.v > cdjunk.v
		vp_annotate <cdjunk.v text=cditer.p batch=y >cditer.v

${RESDIR}/diritr.v: cgiter.v cditer.v
		vp_SideBySideAniso cgiter.v cditer.v >${RESDIR}/diritr.v


cmp.v:		cmp_tpm.H
		Byte < cmp_tpm.H >cmp.A
		<cmp.A Ta2vplot title="Input CMP gather" \
		out=cmp.v 			> /dev/null

vscan.v:	cdscan_w.H
		Byte < cdscan_w.H >vscan.A
		<vscan.A Ta2vplot title="L2 inverse" \
		out=vscan.v > /dev/null

${RESDIR}/dircvv.v:	cmp.v  vscan.v
		vp_SideBySideAniso cmp.v  vscan.v > ${RESDIR}/dircvv.v

res.v:		cmp_tpm.H cdcmp_w.H
		Add cmp_tpm.H cdcmp_w.H scale=1,-1 > junk.H
		Byte < junk.H clip=184.167 >res.A
		<res.A Ta2vplot \
		title="L2 residual data" out=res.v 	> /dev/null

mod.v:		cdcmp_w.H
		Byte < cdcmp_w.H >mod.A
		<mod.A Ta2vplot \
		title="L2 inverse + modeling" out=mod.v	> /dev/null

${RESDIR}/dirrst.v:	mod.v  res.v
		vp_SideBySideAniso mod.v res.v > ${RESDIR}/dirrst.v


# Tpow the CMP gather 807 at well4 location...
cmp_tpm.H:	${DATDIR}/cmp807_raw.HH
		Tpow < ${DATDIR}/cmp807_raw.HH tpow=2 >tpow.H 
		<tpow.H	Mute vmute=1.3 tramp=0.2 > cmp_tpm.H

int.v: ${DATDIR}/wind18.HH ${BINDIR}/Levint.x  interp.p 
	< ${DATDIR}/wind18.HH Transp | ${BINDIR}/Levint.x \
	par=interp.p warmup=20 cdstep=1 > int.H 
	< int.H Byte pclip=100. > int.A
	< int.A ${TA2VPLOT} title="Conjugate directions" out=int.v

lint.v: ${DATDIR}/wind18.HH ${BINDIR}/Levint.x  interp.p 
	<${DATDIR}/wind18.HH Transp | ${BINDIR}/Levint.x \
	par=interp.p warmup=20 cdstep=0 >int.H 
	< int.H Byte pclip=100. > int.A
	< int.A ${TA2VPLOT} title="Conjugate gradients" out=lint.v

${RESDIR}/dirjbm.v: int.v lint.v
	vp_SideBySideAniso lint.v int.v >${RESDIR}/dirjbm.v

clean: jclean

include ${SEPINC}/SEP.bottom
