include ${SEPINC}/SEP.top

VPLOTLIBF90 = -lvplotf90
UF90LIBS = ${BEILIB90} ${VPLOTLIBF90}

RESDIR = ./Fig

RESULTSNR = amoapp amotta amotra amomod amomig amosmo amosmi \
	    amospe amormo amormi amormm amoflt amoimr amocon \
	    amox12 amolow amosft amochp \
            amoimp 
RESULTSER = imp-noa imp-aal imp-fla imp-dip imp-all

LATOPTS = class=geophysics

PSTEXOPTS = invras=n color=n
BINDIR = .

CONT = Contour nc=50 yreverse=n transp=n >/dev/null
DOTS = Dots gaineach=0 dots=0 seemean=1 strings=1 constsep=1 overlap=5.5 transp=1 >/dev/null

default: view

aamod.H:	${BINDIR}/Image.x aamod.p
	${BINDIR}/Image.x par=aamod.p                  > aamod.H

aamkn.H:	aamod.H ${BINDIR}/Kirchnew.x ${BINDIR}/Halfdif.x
	<aamod.H ${BINDIR}/Halfdif.x adj=0 | \
	${BINDIR}/Kirchnew.x adj=0 v0=3 >aamkn.H

aakne.H:	aamkn.H ${BINDIR}/Kirchnew.x ${BINDIR}/Halfdif.x 
	<aamkn.H ${BINDIR}/Kirchnew.x adj=1 v0=3| \
	${BINDIR}/Halfdif.x adj=1  >aakne.H

aamki.H:	aamod.H ${BINDIR}/Kirchhoff.x ${BINDIR}/Halfdif.x 
	<aamod.H ${BINDIR}/Halfdif.x adj=0 | \
	${BINDIR}/Kirchhoff.x adj=0 v0=3 > aamki.H

aakir.H:	aamki.H ${BINDIR}/Kirchhoff.x ${BINDIR}/Halfdif.x 
	<aamki.H ${BINDIR}/Kirchhoff.x adj=1 v0=3 | \
	${BINDIR}/Halfdif.x adj=1   > aakir.H

aamka.H:	aamod.H ${BINDIR}/Kaafast.x ${BINDIR}/Halfdif.x 
	<aamod.H ${BINDIR}/Halfdif.x adj=0 | \
	${BINDIR}/Kaafast.x aal=2 adj=0 add=0 v0=3 > aamka.H

aakaa.H:	aamka.H ${BINDIR}/Kaafast.x ${BINDIR}/Halfdif.x 
	<aamka.H ${BINDIR}/Kaafast.x aal=2 adj=1 add=0 v0=3 | \
	${BINDIR}/Halfdif.x adj=1   > aakaa.H

${RESDIR}/amomod.v: aamod.H aamki.H aamkn.H aamka.H 
	<aamod.H ${DOTS} out=mod.v
	<aamki.H ${DOTS} out=mki.v
	<aamkn.H ${DOTS} out=mkn.v
	<aamka.H ${DOTS} out=mka.v
	vp_SideBySideAniso mod.v mki.v > top.v
	vp_SideBySideAniso mkn.v mka.v > bot.v
	vp_OverUnderAniso top.v bot.v > ${RESDIR}/amomod.v

${RESDIR}/amomig.v: aamod.H aakir.H aakne.H aakaa.H 
	<aamod.H ${DOTS} out=mod.v
	<aakir.H ${DOTS} out=gki.v
	<aakne.H ${DOTS} out=gkn.v
	<aakaa.H ${DOTS} out=gka.v
	vp_SideBySideAniso mod.v gki.v > top.v
	vp_SideBySideAniso gkn.v gka.v > bot.v
	vp_OverUnderAniso top.v bot.v > ${RESDIR}/amomig.v

${RESDIR}/amoapp.v: app.p ${BINDIR}/App.x 
	${BINDIR}/App.x par=app.p h2=1000  fb=10 >a11.H
	${BINDIR}/App.x par=app.p h2=500 fb=10 >a12.H
	${BINDIR}/App.x par=app.p h2=250  fb=10 >a13.H
	${BINDIR}/App.x par=app.p h2=1000  fb=30 >a21.H
	${BINDIR}/App.x par=app.p h2=500 fb=30 >a22.H
	${BINDIR}/App.x par=app.p h2=250  fb=30 >a23.H
	${BINDIR}/App.x par=app.p h2=1000  fb=60 >a31.H
	${BINDIR}/App.x par=app.p h2=500 fb=60 >a32.H
	${BINDIR}/App.x par=app.p h2=250  fb=60 >a33.H
	<a11.H ${CONT} title=" " out=a11.v labelsz=18 
	<a12.H ${CONT} title=" " out=a12.v labelsz=18 
	<a13.H ${CONT} title=" " out=a13.v labelsz=18 
	<a21.H ${CONT} title=" " out=a21.v labelsz=18 
	<a22.H ${CONT} title=" " out=a22.v labelsz=18 
	<a23.H ${CONT} title=" " out=a23.v labelsz=18 
	<a31.H ${CONT} title=" " out=a31.v labelsz=18 
	<a32.H ${CONT} title=" " out=a32.v labelsz=18 
	<a33.H ${CONT} title=" " out=a33.v labelsz=18 
	vp_SideBySideAniso a11.v a12.v a13.v >top.v
	vp_SideBySideAniso a21.v a22.v a23.v >mid.v
	vp_SideBySideAniso a31.v a32.v a33.v >bot.v
	vp_OverUnderAniso top.v mid.v bot.v > ${RESDIR}/amoapp.v

${RESDIR}/amotta.v:	app.p ${BINDIR}/Time.x ${BINDIR}/Amp.x
	${BINDIR}/Time.x par=app.p h2=500  fb=30 >time.H
	<time.H   ${CONT} title="AMO traveltime" out=time.v
	${BINDIR}/Amp.x par=app.p h2=500  fb=30 >amp.H
	<amp.H   ${CONT} title="AMO amplitude" out=amp.v
	vp_SideBySideAniso time.v amp.v >${RESDIR}/amotta.v

LABELS = labels='in:C in:C`C in:CC`C in:C`CC`C in' labelsz=10
DOTARGS = yreverse=1 overlap=.8 constsep=1 Xdots=0 connect=0
${RESDIR}/amoflt.v : ${BINDIR}/Causint.x 
	${BINDIR}/Causint.x hello > junk.H
	<junk.H Dots ${LABELS} ${DOTARGS} >/dev/null out=${RESDIR}/amoflt.v

mod.H : ${BINDIR}/Sigmoid.x     kirsig.p     3.p
	 ${BINDIR}/Sigmoid.x par=kirsig.p par=3.p > junk.H
	<junk.H Smooth par=kirsig.p diff1=1 | Smooth par=kirsig.p >mod.H

mmod.H: mod.H ${BINDIR}/Halfdif.x
	<mod.H ${BINDIR}/Halfdif.x adj=0  >mmod.H

dki.H:  mmod.H ${BINDIR}/Kirchhoff.x
	<mmod.H ${BINDIR}/Kirchhoff.x adj=0   > dki.H

kir.H:	 dki.H ${BINDIR}/Kirchhoff.x ${BINDIR}/Halfdif.x
	<dki.H  ${BINDIR}/Kirchhoff.x adj=1 | \
	${BINDIR}/Halfdif.x adj=1   > kir.H

dka.H:  mmod.H ${BINDIR}/Kaafast.x
	<mmod.H ${BINDIR}/Kaafast.x aal=2 adj=0 add=0 velocity=0.75 > dka.H

kaa.H:	 dka.H ${BINDIR}/Kaafast.x ${BINDIR}/Halfdif.x
	<dka.H  ${BINDIR}/Kaafast.x aal=2 adj=1 add=0 velocity=0.75| \
	${BINDIR}/Halfdif.x adj=1  > kaa.H

dkn.H:  mmod.H ${BINDIR}/Kirchnew.x
	<mmod.H ${BINDIR}/Kirchnew.x adj=0 > dkn.H

kne.H:	 dkn.H ${BINDIR}/Kirchnew.x ${BINDIR}/Halfdif.x
	<dkn.H  ${BINDIR}/Kirchnew.x adj=1 | \
	${BINDIR}/Halfdif.x adj=1   > kne.H

BT = Byte| Ta2vplot labelsz=12 >/dev/null
WINT = Byte| Window n1=100 n2=100 f1=50 f2=50 |Ta2vplot labelsz=16 titlesz=20 >/dev/null
${RESDIR}/amosmo.v:  mod.H kne.H
	<mod.H ${BT} title="model" out=${RESDIR}/amosmo.v

${RESDIR}/amosmi.v:  mod.H kir.H kne.H kaa.H
	<mod.H ${WINT} title="zoomed model" out=mod.v
	<kir.H ${WINT} title="nearest neighbor" out=kir.v
	<kne.H ${WINT} title="reciprocity" out=kne.v
	<kaa.H ${WINT} title="triangles" out=kaa.v
	vp_SideBySideAniso  mod.v kir.v >top.v
	vp_SideBySideAniso  kne.v kaa.v >bot.v
	vp_OverUnderAniso top.v bot.v > ${RESDIR}/amosmi.v

${RESDIR}/amormo.v:	mmod.H ${BINDIR}/Kirchnew.x
	<mmod.H ${BINDIR}/Kirchnew.x adj=0 sw=1 > dk1.H
	<mmod.H ${BINDIR}/Kirchnew.x adj=0 sw=2 > dk2.H
	<dk1.H ${BT} out=dk1.v title="flat data"
	<dk2.H ${BT} out=dk2.v title="steep data"
	vp_SideBySideAniso  dk1.v dk2.v > ${RESDIR}/amormo.v

${RESDIR}/amormi.v:	dkn.H ${BINDIR}/Kirchnew.x ${BINDIR}/Halfdif.x
	<dkn.H ${BINDIR}/Kirchnew.x adj=1 sw=1 | \
	${BINDIR}/Halfdif.x adj=1   > kn1.H
	<dkn.H ${BINDIR}/Kirchnew.x adj=1 sw=2 | \
	${BINDIR}/Halfdif.x adj=1   > kn2.H
	<kn1.H ${BT} out=kn1.v title="flat image"
	<kn2.H ${BT} out=kn2.v title="steep image"
	vp_SideBySideAniso  kn1.v kn2.v > ${RESDIR}/amormi.v

${RESDIR}/amormm.v:	dkn.H kne.H
	<dkn.H ${BT} out=dkn.v title="data"
	<kne.H ${BT} out=kne.v title="image"
	vp_SideBySideAniso  dkn.v kne.v > ${RESDIR}/amormm.v

DOTSPEC3 = Dots yreverse=1 labels="m=model:kirchfast:reciprocity:triangles" label1=Hz dots
${RESDIR}/amospe.v:	mod.H kir.H kne.H kaa.H
	<mod.H Spectrum | Scale axis=123  > junk.mod
	<kir.H Spectrum | Scale axis=123  > junk.kir
	<kne.H Spectrum | Scale axis=123  > junk.kne
	<kaa.H Spectrum | Scale axis=123  > junk.kaa
	Merge axis=2 space=n junk.mod junk.kir junk.kne junk.kaa> junk.merge
	<junk.merge  ${DOTSPEC3} > kirknes.V out=${RESDIR}/amospe.v

${RESDIR}/amotra.v:	${BINDIR}/Nmotraj.x 
	${BINDIR}/Nmotraj.x - > amotra.V 
	Cp amotra.V junk.V out=${RESDIR}/amotra.v

${RESDIR}/amoimp.v: ${BINDIR}/Amo.x nmo0.HH amo.p nmo3.HH
	<nmo0.HH ${BINDIR}/Amo.x par=amo.p nop=35 azin=0 azout=30 >amo3.H
	<nmo0.HH Window n1=1 f1=99 |Byte pclip=100 >slice0.A
	<nmo3.HH Window n1=1 f1=99 |Byte pclip=100 >slice3.A
	<amo3.H Window n1=1 f1=99 |Byte pclip=100 >slice.A
	Cat axis=1 slice0.A slice.A slice3.A >tslice.A
	<tslice.A Ta2vplot transp=n yreverse=n  title="Time Slice" >/dev/null out=tslice.v
	<nmo0.HH Window n2=1 f2=32 |Byte pclip=100 >cross0.A
	<nmo3.HH Window n2=1 f2=32 |Byte pclip=100 >cross3.A
	<amo3.H Window n2=1 f2=32 |Byte pclip=100 >cross.A
	Cat axis=2 cross0.A cross.A cross3.A >cline.A
	<cline.A Ta2vplot title="Cross-line" >/dev/null out=cline.v
	<nmo0.HH Window n3=1 f3=32 |Byte pclip=100 >line0.A
	<nmo3.HH Window n3=1 f3=32 |Byte pclip=100 >line3.A
	<amo3.H Window n3=1 f3=32 |Byte pclip=100 >line.A
	Cat axis=2 line0.A line.A line3.A >iline.A
	<iline.A Ta2vplot title="In-line" >/dev/null out=iline.v
	vp_OverUnderAniso iline.v tslice.v cline.v > ${RESDIR}/amoimp.v

imr.H:			${BINDIR}/Amo.x imp.p
	Spike par=imp.p >spike.H
	<spike.H ${BINDIR}/Amo.x mint=100 hin=500 hout=500 azin=0 azout=30 vel=500 >imr.H nop=35 type=0

imr1.H:			${BINDIR}/Amo.x imp.p
	Spike par=imp.p >spike.H
	<spike.H ${BINDIR}/Amo.x mint=100 hin=500 hout=500 azin=0 azout=30 vel=500 >imr1.H nop=35 type=1

imr2.H:			${BINDIR}/Amo.x imp.p
	Spike par=imp.p >spike.H
	<spike.H ${BINDIR}/Amo.x mint=100 hin=500 hout=500 azin=0 azout=30 vel=500 >imr2.H nop=35 type=2

${RESDIR}/amoimr.v:	imr.H
	<imr.H Window n1=100 f1=78 >impr.H
	<impr.H Window n2=1 f2=7 | Wiggle pclip=100 transp=y yreverse=y poly=y title="Crossline 1" labelsz=18 titlesz=20 >/dev/null out=cross1.v
	<impr.H Window n2=1 f2=12 | Wiggle pclip=100 transp=y yreverse=y poly=y title="Crossline 2" labelsz=18 titlesz=20 >/dev/null out=cross2.v
	<impr.H Window n2=1 f2=17 | Wiggle pclip=100 transp=y yreverse=y poly=y title="Crossline 3" labelsz=18 titlesz=20 >/dev/null out=cross3.v
	<impr.H Window n3=1 f3=7 | Wiggle pclip=100 transp=y yreverse=y poly=y title="Inline 1" labelsz=18 titlesz=20 >/dev/null out=in1.v
	<impr.H Window n3=1 f3=12 | Wiggle pclip=100 transp=y yreverse=y poly=y title="Inline 2" labelsz=18 titlesz=20 >/dev/null out=in2.v
	<impr.H Window n3=1 f3=17 | Wiggle pclip=100 transp=y yreverse=y poly=y title="Inline 3" labelsz=18 titlesz=20 >/dev/null out=in3.v
	vp_SideBySideAniso in1.v in2.v in3.v  >in.v
	vp_SideBySideAniso cross1.v cross2.v cross3.v  >cross.v
	vp_OverUnderAniso in.v cross.v >${RESDIR}/amoimr.v

${RESDIR}/amocon.v:	imr.H imr1.H imr2.H
	<imr.H Transp plane=13 |Stack norm=no >foot.H
	<imr1.H Transp plane=13 |Stack norm=no >foot1.H
	<imr2.H Transp plane=13 |Stack norm=no >foot2.H
	<foot.H Scale axis=12| Byte allpos=y|Ta2vplot >/dev/null labelsz=15 titlesz=18 title="total" out=tot.v
	<foot1.H Scale axis=12| Byte allpos=y|Ta2vplot >/dev/null labelsz=15 titlesz=18 title="flat part" out=flat.v 
	<foot2.H Scale axis=12| Byte allpos=y|Ta2vplot >/dev/null labelsz=15 titlesz=18 title="steep part" out=steep.v 
	vp_SideBySideAniso flat.v steep.v tot.v >${RESDIR}/amocon.v
	vp_SideBySideAniso steep.v tot.v >${RESDIR}/amocon.v

CLIPIMP=95
LABT = label1="Time (s)" label2="CMP X (m)" labelsz=8 labelfat=2
LABSL = label1="CMP Y (s)" label2="CMP X (m)" labelsz=8 labelfat=2
VSY=0
VST=6.
VPPST = vpstyle=no big=no 

imp-aal.H: Mig.x trace.HH hdr.H mig.p
	< trace.HH Window min1=0 | Mig.x par=mig.p hdr=hdr.H > $@
	echo "o2=0. o3=0." >>$@

imp-noa.H: Mig.x trace.HH hdr.H mig.p
	< trace.HH Window min1=0 | Mig.x par=mig.p hdr=hdr.H antialias=no > $@
	echo "o2=0. o3=0." >>$@

imp-fla.H: Migflat.x trace.HH hdr.H mig.p
	< trace.HH Window min1=0 | Migflat.x par=mig.p hdr=hdr.H > $@
	echo "o2=0. o3=0." >>$@

imp-dip.H: Migsteep.x trace.HH hdr.H mig.p
	< trace.HH Window min1=0 | Migsteep.x par=mig.p hdr=hdr.H > $@
	echo "o2=0. o3=0." >>$@

imp-all.H: imp-fla.H imp-dip.H
	Add imp-fla.H imp-dip.H > $@

${RESDIR}/imp-%.v: imp-%.H
	<imp-$*.H Window n1=1 min1=1 | Transp plane=12 \
	|Grey pclip=$(CLIPIMP) >t.V \
	title=" " $(LABSL) titlefat=0 grid1=n d1num=1000 
	<imp-$*.H Window n3=1 min3=1000 \
	| Grey pclip=$(CLIPIMP) >y.V \
	title=" " $(LABT) titlefat=0 grid1=n d1num=1 wantaxis2=n
	Vppen <t.V >/dev/null $(VPPST) vshift=$(VST) xscale=.7 yscale=.35 \
	out=t.v
	Vppen <y.V >/dev/null $(VPPST) vshift=$(VSY) scale=.7  \
	out=y.v
	vppen < t.v y.v $(VPPST) erase=once >$@ scale=1.05

hdr.H: hdr.HH
	< hdr.HH Dd esize=4 | Window > $@

clean: jclean

include ${SEPINC}/SEP.bottom
