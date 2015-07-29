include ${SEPINC}/SEP.top

UTUBEFLAGS = see_progress=y
UF90LIBS = ${GEELIB}

RESULTSER = data im1 fm1 fm6 pm1 early1 schwab1 \
	mall sall seabdat seabold seabnew \
	waves autowaves thin42 laplac \
	bin500 bin50 norm500 norm50 \
	elev raindata lapinter lapstat precinter
RESULTSCR = helix helix1


RESDIR = ./Fig
DATDIR = ./Dat

SRCDIR = ../Src
BINDIR = ../Bin/${MTYPE}
OBJDIR = ../Obj/${MTYPE}

RAINDIR = ${DATALIB}/fun/chernobyl

MM = min1=-185.556 max1=193.18275 min2=-127.262 max2=127.25044 

LATOPTS = option=minimal

GRAPH = Graph clip=1 max2=1 min2=-1 crowd=.99 > junk.H
GRAPH2 = Graph >/dev/null dash=1,0
DOTS = Dots >/dev/null dots=2 connect=0 strings=1 labelsz=14
DOTY = Dots dots=0 connect=0 gaineach=0 constsep=1 strings=0 clip=1.2 > /dev/null
LBL = label1="Iterations"
TAVPLOT = Ta2vplot crowd=.85  wantaxis=0 >/dev/null
MPAD = prec=0 padin=2 padout=2
DOTARGS2 = dots=2 connect=0 gaineach=0 overlap=1.30 labelsz=16
DOTARGS = dots=0 strings=0 connect=1 constsep=1  \
	labels="model-space:data-space:de-spike:input" \
	labelsz=6 Xclip=2.5 gaineach=0 overlap=.95
PLOTSIG = Dots ${DOTARGS} > /dev/null
CUBE = Cubeplot flat=n  \
	point1=0.9 point2=0.9  frame1=100 frame2=200 frame3=8 >/dev/null
CUBE1 = Cubeplot flat=n  \
	point1=0.9 point2=0.9  frame1=100 frame2=280 frame3=11 >/dev/null

COLOR=y

GREY = Grey gainpanel=every pclip=100 wantaxis=n >/dev/null \
	crowd=.88 gpow=.7 label1=' ' label2=' '

.PRECIOUS: alias.P data.H head.H

default : view

${RESDIR}/data.v alias.P dat.v: ${BINDIR}/Alias.x ${BINDIR}/Ideal.x
	 ${BINDIR}/Alias.x tmax=80. pairs=1   n1=200 > alias.P
	 ${BINDIR}/Alias.x tmax=80. pairs=0   n1=200 > alias.H
	<alias.H  ${DOTY}  strings=1 out=dat.v
	${BINDIR}/Ideal.x tmax=80. pairs=0 n1=200 > ideal.H
	< ideal.H ${DOTY} connect=1 out=ideal.v
	vp_Overlay dat.v ideal.v > both.v
	vp_SideBySideAniso both.v dat.v > ${RESDIR}/data.v

aliasbin.H: alias.P ${BINDIR}/Bin.x
	< alias.P ${BINDIR}/Bin.x n1=200 > $@

aliasbin.P: aliasbin.H ${BINDIR}/Aliasreg.x
	<  aliasbin.H ${BINDIR}/Aliasreg.x > $@

aliasbin.v:  aliasbin.H  mlines.H ${BINDIR}/Miss.x
	<aliasbin.H ${BINDIR}/Miss.x filt=mlines.H niter=200 prec=1 >int.H
	< int.H ${DOTY}  connect=1 out=$@

${RESDIR}/bm$*.v ${RESDIR}/bm%.v3 bm%.H: aliasbin.H mlines.H \
	${BINDIR}/Miss.x  dat.v
	<aliasbin.H ${BINDIR}/Miss.x mmov=bm$*.H \
	filt=mlines.H niter=301 prec=1 >int.H
	<bm$*.H Transp plane=23     > transp.H
	<transp.H ${GRAPH} out=${RESDIR}/bm$*.v3 fastplot=20
	<int.H ${DOTY}  connect=1 out=int.v
	vp_Overlay dat.v int.v > ${RESDIR}/bm$*.v

${RESDIR}/im%.v ${RESDIR}/im%.v3 im%.H: \
	dat.v alias.P  ${BINDIR}/Invint.x  
	<alias.P ${BINDIR}/Invint.x n1=200 mm=im$*.H > int.H method=1 filter=$*
	<im$*.H Transp plane=23     > transp.H
	<transp.H ${GRAPH} out=${RESDIR}/im$*.v3 fastplot=20
	<int.H ${DOTY}  connect=1 out=int.v
	vp_Overlay dat.v int.v > ${RESDIR}/im$*.v

${RESDIR}/fm%.v ${RESDIR}/fm%.v3 fm%.H: \
	dat.v  alias.P ${BINDIR}/Invint.x  
	<alias.P ${BINDIR}/Invint.x n1=200 mm=fm$*.H > int.H \
	method=2 filter=$* eps=.01
	<fm$*.H Transp plane=23     > transp.H
	<transp.H ${GRAPH} out=${RESDIR}/fm$*.v3 fastplot=20
	<int.H ${DOTY}  connect=1 out=int.v
	vp_Overlay dat.v int.v > ${RESDIR}/fm$*.v

${RESDIR}/pm%.v ${RESDIR}/pm%.v3 pm%.H: \
	dat.v  alias.P ${BINDIR}/Invint.x  
	<alias.P ${BINDIR}/Invint.x n1=200 mm=pm$*.H > int.H \
	method=2 filter=$* pef=1 eps=.01
	<pm$*.H Transp plane=23     > transp.H
	<transp.H ${GRAPH} out=${RESDIR}/pm$*.v3 fastplot=20
	<int.H ${DOTY}  connect=1 out=int.v
	vp_Overlay dat.v int.v > ${RESDIR}/pm$*.v

NITER = 250
ITERATION = 5
${RESDIR}/early%.v : im%.H fm%.H    
	< im$*.H Window n2=1 f2=${ITERATION}   > iwind.H
	< fm$*.H Window n2=1 f2=${ITERATION}   > fwind.H
	< fm$*.H Window n2=1 f2=${NITER}       > icorr.H
	Merge <iwind.H fwind.H         icorr.H space=no axis=2 > junk.H
	< junk.H ${DOTY} connect=1 overlap=2.0 title=" " \
	labels=regul.${ITERATION}:precon.${ITERATION}:niter=${NITER} \
	labelsz=10 out=${RESDIR}/early$*.v

${RESDIR}/schwab%.v :  ${BINDIR}/Modnorm.x  im%.H fm%.H 
	<im$*.H ${BINDIR}/Modnorm.x > junk1.H
	<fm$*.H ${BINDIR}/Modnorm.x > junk2.H
	Merge <junk1.H junk2.H         space=n axis=1 > junk.H
	< junk.H Transp plane=12  > transp.H
	<transp.H Window f1=0 n1=100 > wind.H
	< wind.H Graph ${LBL} symbol='md' symbolsz=10 symbolsz=7 \
	title="Model Residual Vector Length" \
	min2=0 min1=-1 out=${RESDIR}/schwab$*.v >/dev/null

${RESDIR}/beam%.v :  ${BINDIR}/Modnorm.x  bm%.H fm%.H 
	<bm$*.H ${BINDIR}/Modnorm.x > junk1.H
	<fm$*.H ${BINDIR}/Modnorm.x > junk2.H
	Merge <junk1.H junk2.H         space=n axis=1 > junk.H
	< junk.H Transp plane=12  > transp.H
	<transp.H Window f1=0 n1=40 > wind.H
	< wind.H Graph ${LBL} symbol='bd' symbolsz=10 symbolsz=7 \
	title="Model Residual Vector Length" \
	min2=0 min1=-1 out=${RESDIR}/beam$*.v >/dev/null

#-------------- 3D ------------------------------------

head.H:  ${DATDIR}/head.H
	Window < ${DATDIR}/head.H n2=2 > head.H

dat.H:  ${DATDIR}/dat.H
	Window < ${DATDIR}/dat.H n2=150 f2=50 > dat.H

offset.H: ${DATDIR}/head.H
	Window < ${DATDIR}/head.H n2=2 f2=2 > $@

aoff.H: ${DATDIR}/head.H
	Window < ${DATDIR}/head.H n2=1 f2=2 > $@

head1.H: head.H aoff.H ${BINDIR}/Select.x
	< head.H ${BINDIR}/Select.x header=aoff.H hmin=500 hmax=600 > $@

dat1.H: dat.H aoff.H ${BINDIR}/Select.x
	< dat.H ${BINDIR}/Select.x header=aoff.H hmin=500 hmax=600 > $@

${RESDIR}/offset.v:  offset.H
	Window <   offset.H n2=1 f2=0 > offx.H
	Window <   offset.H n2=1 f2=1 > offy.H
	Cmplx offx.H offy.H > off.H 
	Graph < off.H  symbol='+' symbolsz=2 out=${RESDIR}/offset.v \
	label1="Offset" label2="Azimuth" title="Offset Geometry"  > /dev/null

${RESDIR}/cmp.v:  head.H
	Window <   head.H n2=1 f2=0 > cmpx.H
	Window <   head.H n2=1 f2=1 > cmpy.H
	Cmplx cmpx.H cmpy.H > cmp.H 
	Graph < cmp.H  symbol='+' symbolsz=2 out=${RESDIR}/cmp.v \
	label1="In-line midpoint" label2="Cross-line midpoint" title="CMP Geometry"  > /dev/null

${RESDIR}/cmp1.v:  head1.H
	Window <   head1.H n2=1 f2=0 > cmpx.H
	Window <   head1.H n2=1 f2=1 > cmpy.H
	Cmplx cmpx.H cmpy.H > cmp.H 
	Graph < cmp.H  symbol='+' symbolsz=2 out=${RESDIR}/cmp1.v \
	label1="In-line midpoint" label2="Cross-line midpoint" title="CMP Geometry"  > /dev/null

bin.H fold.H: grid.p dat.H head.H ${BINDIR}/Bin2.x
	${BINDIR}/Bin2.x \
	par=grid.p hff=head.H fold=fold.H interp=1 <dat.H > junk.H
	<junk.H Transp plane=31 > bin.H

bin-1.H fold-1.H: grid1.p dat1.H head1.H ${BINDIR}/Bin2.x
	${BINDIR}/Bin2.x \
	par=grid1.p hff=head1.H fold=fold-1.H interp=1 <dat1.H > junk.H
	<junk.H Transp plane=31 > bin-1.H

bin%.H: grid.p dat.H head.H ${BINDIR}/Bin2.x
	${BINDIR}/Bin2.x \
	par=grid.p hff=head.H interp=$* inv=1 <dat.H > bin$*.H

lin.H lfold.H: grid.p dat.H head.H ${BINDIR}/Bin2.x
	${BINDIR}/Bin2.x \
	par=grid.p hff=head.H fold=lfold.H interp=2 <dat.H > junk.H
	<junk.H Transp plane=31 > lin.H

elin.H: grid.p dat.H head.H ${BINDIR}/Bin2.x
	${BINDIR}/Bin2.x par=grid.p hff=head.H interp=2 inv=1 <dat.H > elin.H

inv%.H: grid.p dat.H head.H ${BINDIR}/Inv2.x
	${BINDIR}/Inv2.x \
	par=grid.p hff=head.H niter=5 interp=$* inv=2 <dat.H > junk.H
	<junk.H Transp plane=31 > inv$*.H

in1-%.H: grid1.p dat1.H head1.H ${BINDIR}/Inv2.x
	${BINDIR}/Inv2.x \
	par=grid1.p hff=head1.H niter=5 interp=$* inv=2 <dat1.H > junk.H
	<junk.H Transp plane=31 > in1-$*.H

in2-%.H: grid1.p dat1.H head1.H ${BINDIR}/Inv2.x
	${BINDIR}/Inv2.x \
	par=grid1.p hff=head1.H niter=5 interp=$* inv=3 <dat1.H > junk.H
	<junk.H Transp plane=31 > in2-$*.H

vbin.H: grid.p dat.H head.H ${BINDIR}/Inv2.x
	${BINDIR}/Inv2.x \
	par=grid.p hff=head.H niter=5 interp=1 inv=1 <dat.H > junk.H
	<junk.H Transp plane=31 > vbin.H

vbin-1.H: grid1.p dat1.H head1.H ${BINDIR}/Inv2.x
	${BINDIR}/Inv2.x \
	par=grid1.p hff=head1.H niter=5 interp=1 inv=1 <dat1.H > junk.H
	<junk.H Transp plane=31 > vbin-1.H

mres.H: grid.p dat.H head.H ${BINDIR}/Inv2.x
	${BINDIR}/Inv2.x \
	par=grid.p hff=head.H niter=5 interp=1 inv=0 <dat.H > mres.H

${RESDIR}/fold.v: fold.H
	<fold.H Taplot allpos=y pclip=100 > junk.A
	<junk.A Ta2vplot transp=n yreverse=n wantscalebar=1 \
	labelsz=10 title="Fold Map" out=${RESDIR}/fold.v > /dev/null

${RESDIR}/fold1.v: fold-1.H
	<fold-1.H Taplot allpos=y pclip=100 > junk.A
	<junk.A Ta2vplot transp=n yreverse=n wantscalebar=1 \
	labelsz=10 title="Fold Map" out=${RESDIR}/fold1.v > /dev/null

${RESDIR}/lfold.v: lfold.H
	<lfold.H Taplot allpos=y pclip=100 > junk.A
	<junk.A Ta2vplot transp=n yreverse=n wantscalebar=1 \
	labelsz=10 title="Linear Fold Map" out=${RESDIR}/lfold.v > /dev/null

${RESDIR}/hist.v: fold.H ${BINDIR}/Histogram.x
	<fold.H ${BINDIR}/Histogram.x > hist.H
	<hist.H Graph title="Fold Distribution" >/dev/null out=${RESDIR}/hist.v

${RESDIR}/hist1.v: fold-1.H ${BINDIR}/Histogram.x
	<fold-1.H ${BINDIR}/Histogram.x > hist.H
	<hist.H Graph title="Fold Distribution" >/dev/null \
	out=${RESDIR}/hist1.v

${RESDIR}/lhist.v: lfold.H ${BINDIR}/Histogram.x
	<lfold.H ${BINDIR}/Histogram.x > lhist$*.H
	<lhist.H Graph title="Linear Fold Distribution" >/dev/null out=${RESDIR}/lhist.v

${RESDIR}/bin.v: bin.H
	<bin.H Byte > bin.A
	<bin.A ${CUBE} title="Nearest Neighbor" out=${RESDIR}/bin.v

${RESDIR}/bin1.v: bin-1.H
	<bin-1.H Byte > bin.A
	<bin.A ${CUBE1} title=" " out=${RESDIR}/bin1.v

${RESDIR}/regbin.v: inv1.H
	<inv1.H Byte > inv1.A
	<inv1.A ${CUBE} title="Nearest Neighbor: Regularized Optimization" out=${RESDIR}/regbin.v

${RESDIR}/regbin1.v: in1-1.H
	<in1-1.H Byte > inv1.A
	<inv1.A ${CUBE1} title="Nearest Neighbor: Regularized Optimization" out=${RESDIR}/regbin1.v

${RESDIR}/invbin.v: vbin.H
	<vbin.H Byte > vbin.A
	<vbin.A ${CUBE} title="Nearest Neighbor: Optimization" out=${RESDIR}/invbin.v

${RESDIR}/invbin1.v: vbin-1.H
	<vbin-1.H Byte > vbin.A
	<vbin.A ${CUBE1} title=" " out=${RESDIR}/invbin1.v

${RESDIR}/lin.v: lin.H
	<lin.H Byte > lin.A
	<lin.A ${CUBE} title="Linear" out=${RESDIR}/lin$*.v

${RESDIR}/elin.v: elin.H
	<elin.H Transp plane=31 > junk.H
	<junk.H Byte > elin.A
	<elin.A ${CUBE} title="Linear: Data normalized" out=${RESDIR}/elin.v

${RESDIR}/reglin.v: inv2.H
	<inv2.H Byte > inv2.A
	<inv2.A ${CUBE} title="Linear: Regularized Optimization" out=${RESDIR}/reglin.v

${RESDIR}/reglin1.v: in1-2.H
	<in1-2.H Byte > inv2.A
	<inv2.A ${CUBE} title="Linear: Regularized Optimization" out=${RESDIR}/reglin1.v

${RESDIR}/regpef1.v: in2-2.H
	<in2-2.H Byte > inv3.A
	<inv3.A ${CUBE} title=" " out=${RESDIR}/regpef1.v

${RESDIR}/regmin.v: inv-4.H
	<inv-4.H Byte > inv-4.A
	<inv-4.A ${CUBE} title="Muir: Regularized Optimization" out=${RESDIR}/regmin.v

${RESDIR}/regtin.v: inv4.H
	<inv4.H Byte > inv4.A
	<inv4.A ${CUBE} title="Lagrange: Regularized Optimization" out=${RESDIR}/regtin.v

sall.% : 
	@-cd miss; $(MAKE) sall.$*;

mall.% : 
	@-cd miss; $(MAKE) mall.$*;

seabdat.% : 
	@-cd seab; $(MAKE) seabdat.$*;

seabold.% : 
	@-cd seab; $(MAKE) seabold.$*;

seabnew.% : 
	@-cd seab; $(MAKE) seabnew.$*;

${RESDIR}/auto%.v: \
	${BINDIR}/Helicon.x ${BINDIR}/Autocorr.x ${BINDIR}/Wilson.x %.H %.P
	Spike n1=10 n2=10 k1=5 k2=5 d1=1 d2=1 > spik.H
	<$*.H ${BINDIR}/Autocorr.x > auto$*.H
	<spik.H ${BINDIR}/Helicon.x adj=0 div=0 filt=$*.H > filt.H
	<filt.H ${BINDIR}/Helicon.x adj=1 div=0 filt=$*.H > auto.H
	<auto$*.H ${BINDIR}/Wilson.x par=$*.P >fact$*.H
	<spik.H ${BINDIR}/Helicon.x adj=0 div=0 filt=fact$*.H > fact.H
	<filt.H ${BINDIR}/Helicon.x adj=0 div=1 filt=fact$*.H > spik2.H
	<filt.H  ${GREY} out=filt.v  title=filter
	<auto.H  ${GREY} out=auto.v  title=autocorrelation
	<fact.H  ${GREY} out=fact.v  title="Wilson factor"
	<spik2.H ${GREY} out=spik.v  title="Wilson inversion"
	vp_SideBySideAniso filt.v auto.v > top.v 
	vp_SideBySideAniso fact.v spik.v > bot.v
	vp_OverUnderAniso top.v bot.v > junk.v
	<junk.v  vppen vpstyle=n txscale=2.4 > ${RESDIR}/auto$*.v

${RESDIR}/laplac.v : \
	${BINDIR}/Helicon.x ${BINDIR}/Autocorr.x ${BINDIR}/Wilson.x \
	thinfac42.H thinplate.H laplac.H
	Spike n1=40 n2=40 k1=20 k2=20 d1=1 d2=1 > spik.H
	<spik.H ${BINDIR}/Helicon.x adj=0 div=0 filt=laplac.H > lap1.H
	<spik.H ${BINDIR}/Helicon.x adj=1 div=0 filt=laplac.H > lap2.H
	Add lap1.H lap2.H > lap.H
	<spik.H ${BINDIR}/Helicon.x adj=0 div=0 filt=thinplate.H > bih1.H
	<spik.H ${BINDIR}/Helicon.x adj=1 div=0 filt=thinplate.H > bih2.H
	Add bih1.H bih2.H > bih.H
	<spik.H ${BINDIR}/Helicon.x adj=0 div=0 filt=thinfac42.H > fact.H
	<bih.H  ${BINDIR}/Helicon.x adj=0 div=1 filt=thinfac42.H > spk1.H
	<spk1.H ${BINDIR}/Helicon.x adj=1 div=1 filt=thinfac42.H > spk2.H
	<lap.H   Window n1=10 n2=10 f1=15 f2=15 | ${GREY} out=filt.v  title=laplacian
	<bih.H   Window n1=10 n2=10 f1=15 f2=15 | ${GREY} out=auto.v  title=autocorrelation
	<fact.H  Window n1=10 n2=10 f1=15 f2=15 | ${GREY} out=fact.v  title="Wilson factor"
	<spk2.H  Window n1=10 n2=10 f1=15 f2=15 | ${GREY} out=spik.v  title="Wilson inversion"
	vp_SideBySideAniso filt.v auto.v > top.v 
	vp_SideBySideAniso fact.v spik.v > bot.v
	vp_OverUnderAniso top.v bot.v > junk.v
	<junk.v  vppen vpstyle=n txscale=2.4 > $@

%.H: %.HH
	<$*.HH Dd esize=4 | Window >$*.H
%.P: %.H
	grep lag $*.H > $*.P
	grep n1 $*.H >> $*.P

# Mathematica creates a wrong Bounding Box, and we need to fix it with sed
${RESDIR}/helix1.ps: ${RESDIR}/helix.ps
	cat ${RESDIR}/helix.warning > ${RESDIR}/helix1.warning	
	cat helix.s  > $@	
	sed -e 1,2d ${RESDIR}/helix.ps >>$@

border.V: ${RAINDIR}/border.HH
	< ${RAINDIR}/border.HH Window n1=1 f1=0 n2=338 > x.H
	< ${RAINDIR}/border.HH Window n1=1 f1=1 n2=338 > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} title=" " plotcol=7 plotfat=2 > border1.V
	< ${RAINDIR}/border.HH Window n1=1 f1=0 f2=338 n2=234 > x.H
	< ${RAINDIR}/border.HH Window n1=1 f1=1 f2=338 n2=234 > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} title=" " plotcol=7 plotfat=2 > border2.V
	< ${RAINDIR}/border.HH Window n1=1 f1=0 f2=572  > x.H
	< ${RAINDIR}/border.HH Window n1=1 f1=1 f2=572 > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} title=" " plotcol=7 plotfat=2 > border3.V
	Vppen < border1.V border2.V border3.V \
	erase=once vpstyle=n big=n >$@

${RESDIR}/elev.v: ${RAINDIR}/elevation.HH ${BINDIR}/Igrad.x border.V
	< ${RAINDIR}/elevation.HH ${BINDIR}/Igrad.x > elev.H
	< elev.H Grey title=Elevation transp=n yreverse=n \
	wantaxis=n wantlabel=n wheretitle=t wherexlabel=b > elev.V
	Vppen < elev.V border.V  \
	erase=once vpstyle=n big=n >/dev/null out=$@

${RESDIR}/alldata.v: ${RAINDIR}/alldata.HH
	< ${RAINDIR}/alldata.HH Window n1=1 f1=0 > x.H
	< ${RAINDIR}/alldata.HH Window n1=1 f1=1  > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} symbol=x title="All data locations" \
	>/dev/null out=$@

obsdata.v: ${RAINDIR}/obsdata.HH border.V
	< ${RAINDIR}/obsdata.HH Window n1=1 f1=0  > x.H
	< ${RAINDIR}/obsdata.HH Window n1=1 f1=1  > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} symbol=o title="Observed data locations" \
	plotcol=7 >obs.V
	Vppen < obs.V border.V \
	erase=once vpstyle=n big=n >/dev/null out=$@

data.v: ${RAINDIR}/obsdata.HH ${RAINDIR}/coord.HH border.V
	< ${RAINDIR}/obsdata.HH Window n1=1 f1=0 > x.H
	< ${RAINDIR}/obsdata.HH Window n1=1 f1=1 > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} symbol=o title="Data locations" \
	plotcol=7 >obs.V
	< ${RAINDIR}/coord.HH Window n1=1 f1=0 > x.H
	< ${RAINDIR}/coord.HH Window n1=1 f1=1 > y.H
	Cmplx x.H y.H > xy.H
	< xy.H Graph ${MM} symbol=x title=" " >coord.V
	Vppen < obs.V coord.V border.V  \
	erase=once vpstyle=n big=n >/dev/null out=$@

${RESDIR}/raindata.v: obsdata.v data.v
	vp_SideBySideIso data.v obsdata.v > $@

${RESDIR}/wilson.txt: ${BINDIR}/Wilsontest.x
	${BINDIR}/Wilsontest.x > $@

lapinter%.H: ${RAINDIR}/elevation.HH ${RAINDIR}/obsdata.HH ${BINDIR}/Lapinter.x
	< ${RAINDIR}/obsdata.HH Window n1=2 | Transp > xy.H
	< ${RAINDIR}/obsdata.HH Window f1=2 > raindat.H
	time < ${RAINDIR}/elevation.HH \
	${BINDIR}/Lapinter.x coord=xy.H dat=raindat.H niter=$* > $@

lapinter%.v: lapinter%.H
	< lapinter$*.H Grey wantscalebar=y yreverse=n transp=n allpos=y \
	color=I pclip=100 title="Laplace $*" >/dev/null out=grey.v
	< lapinter$*.H Contour nc=4 c0=50 dc=100 yreverse=n transp=n \
	wanttitle=n plotcol=7 wantframe=y wantaxis=n \
	screenwd=12.17 >/dev/null out=cont.v
	vp_Overlay grey.v cont.v > $@

#lapinter%.v: lapinter%.H
#	< lapinter$*.H Grey wantscalebar=y yreverse=n transp=n \
#	color=j pclip=100 wheretitle=t wherexlabel=b \
#	title="Model-space: $* iterations" >/dev/null out=$@

${RESDIR}/lapinter.v: lapinter10.v lapinter100.v lapinter1000.v lapinter10000.v
	vppen < lapinter10.v lapinter100.v lapinter1000.v lapinter10000.v \
	vpstyle=n size=r gridnum=2,2 >$@

lapdat%.H: lapinter%.H ${RAINDIR}/coord.HH ${BINDIR}/Extract.x
	< ${RAINDIR}/coord.HH Window | Transp > xy.H
	< lapinter$*.H ${BINDIR}/Extract.x coord=xy.H > lapdat$*.H

${RESDIR}/lapstat.v: lapdat10000.H ${RAINDIR}/predict.HH 
	Cmplx ${RAINDIR}/predict.HH lapdat10000.H > xy.H
	Graph < xy.H symbol=* title="Laplace Correlation" >/dev/null out=$@

tridat.H: tridata.H ${RAINDIR}/coord.HH ${BINDIR}/Extract.x
	< ${RAINDIR}/coord.HH Window | Transp > xy.H
	< tridata.H ${BINDIR}/Extract.x coord=xy.H > $@

${RESDIR}/tristat.v: tridat.H ${RAINDIR}/predict.HH 
	Cmplx ${RAINDIR}/predict.HH tridat.H > xy.H
	Graph < xy.H symbol=* title="Triang Correlation" >/dev/null out=$@

lapstat.v: lapdat10000.H ${RAINDIR}/predict.HH 
	Cmplx ${RAINDIR}/predict.HH lapdat10000.H > xy.H
	Graph < xy.H symbol=* title="Correlation" >/dev/null out=$@

rand500.H: ${BINDIR}/Random.x
	${BINDIR}/Random.x tmax=80. n1=500 nd=50 > $@

rand50.H: ${BINDIR}/Random.x
	${BINDIR}/Random.x tmax=80. n1=50 nd=500 > $@

splbin%.H: rand%.H ${BINDIR}/Bin.x
	< rand$*.H ${BINDIR}/Bin.x n1=$* > $@

${RESDIR}/bin%.v: splbin%.H
	< splbin$*.H ${DOTS} out=$@

envint%.H emm%.H: rand%.H ${BINDIR}/Invinter.x ${BINDIR}/Spline1.x
	${BINDIR}/Spline1.x tension=0.01 > flt.H
	time < rand$*.H ${BINDIR}/Invinter.x niter=20 \
	n1=$* filt=flt.H mm=emm$*.H > envint$*.H

snvint%.H smm%.H: rand%.H ${BINDIR}/Invinter.x \
	${BINDIR}/Spline3.x ${BINDIR}/Splinefilter.x
	${BINDIR}/Spline3.x tension=0.01 > flt.H
	time < rand$*.H ${BINDIR}/Invinter.x niter=20 \
	n1=$* filt=flt.H mm=mm.H spline=y nw=4 > invint.H
	< mm.H ${BINDIR}/Splinefilter.x post=y nw=4 > smm$*.H
	< invint.H ${BINDIR}/Splinefilter.x post=y nw=4 > snvint$*.H

alias%.H: ${BINDIR}/Alias2.x
	${BINDIR}/Alias2.x tmax=80. n1=$* > $@

enorm%.H: emm%.H alias%.H ${BINDIR}/Modnorm1.x
	 < emm$*.H ${BINDIR}/Modnorm1.x ideal=alias$*.H > $@

snorm%.H: smm%.H alias%.H ${BINDIR}/Modnorm1.x
	 < smm$*.H ${BINDIR}/Modnorm1.x ideal=alias$*.H > $@

${RESDIR}/norm%.v: enorm%.H snorm%.H
	Cat axis=2 enorm$*.H snorm$*.H | \
	${GRAPH2} title="Model convergence $*" \
	legend=1 legendloc="tr" legendtitle=" " legendsz=10 \
	label1=Iterations label2="Model misfit" \
	curvelabel="Linear:B-Spline" out=$@

${RESDIR}/waves.v : ${BINDIR}/Helicon.x waves.H 
	Spike n1=20 n2=20 nsp=2 k1=3,7 k2=7,3 mag=1,-1  > spike.H
	Spike n1=20 n2=20 nsp=2 k1=19,15 k2=13,6 mag=1,-1 > inp2.H
	<spike.H ${BINDIR}/Helicon.x filt=waves.H > inp1.H
	Add inp1.H inp2.H > inp.H
	<inp.H ${BINDIR}/Helicon.x   filt=waves.H div=1 > div.H
	<inp.H ${GREY} out=inp.v  title=input
	<div.H ${GREY} out=div.v  title=input/filter
	vp_SideBySideAniso inp.v div.v > junk.v
	<junk.v  vppen vpstyle=n txscale=1.4  > $@

thinfac%.H: thinplate.H ${BINDIR}/Wilson.x fact%.p
	< thinplate.H ${BINDIR}/Wilson.x par=fact$*.p > $@

${RESDIR}/thin%.v : ${BINDIR}/Helicon.x thinfac%.H 
	Spike n1=40 n2=40 nsp=2 k1=11,16 k2=8,3 mag=1,-1  > spike.H
	Spike n1=40 n2=40 nsp=2 k1=31,28 k2=24,16 mag=1,-1 > inp2.H
	<spike.H ${BINDIR}/Helicon.x filt=thinfac$*.H > inp1.H
	Add inp1.H inp2.H > inp.H
	<inp.H ${BINDIR}/Helicon.x       div=1 filt=thinfac$*.H > div.H
	<div.H ${BINDIR}/Helicon.x adj=1 div=1 filt=thinfac$*.H > div2.H
	<inp.H  ${GREY} out=inp.v  clip=1 title=input
	<div.H  ${GREY} out=div.v  clip=1 title=input/filter
	<div2.H ${GREY} out=div2.v pclip=100 title="(input/filter)/filter'"
	vp_SideBySideAniso inp.v div.v div2.v > junk.v
	<junk.v  vppen vpstyle=n txscale=2.4 > $@

precinter%.H: ${RAINDIR}/elevation.HH thinfac42.H \
	      ${RAINDIR}/obsdata.HH ${BINDIR}/Pinter1.x
	< ${RAINDIR}/obsdata.HH Window n1=2 | Transp  > xy.H
	< ${RAINDIR}/obsdata.HH Window f1=2 > dat.H
	time < ${RAINDIR}/elevation.HH ${BINDIR}/Pinter1.x \
	niter=$* coord=xy.H dat=dat.H filt=thinfac42.H > $@

TRIFILES = list_struct delaunay heap predicates
TRIINCLUDS = ${addprefix ${SRCDIR}/, ${addsuffix .h, $(TRIFILES)}}
TRISOURCES = ${addprefix ${SRCDIR}/, ${addsuffix .c, $(TRIFILES)}}
TRIOBJECTS = ${addprefix ${OBJDIR}/, ${addsuffix .o, $(TRIFILES)}}

SEPLIB = -L$(SEPROOT)/lib -lsep
SEPROOT = /usr/local/SEP

# without graphics
${BINDIR}/Sept2a.x: ${SRCDIR}/sept2a.c $(TRISOURCES)
	$(MAKE) CFLAGS="-DWITH_TIME" ${OBJDIR}/sept2a.o $(TRIOBJECTS) 
	$(CC) -o $@ ${OBJDIR}/sept2a.o $(TRIOBJECTS) $(SEPLIB) -lm

# seplib i/o
${OBJDIR}/sept2a.o: ${SRCDIR}/sept2a.c 
	$(CC) -c $(CFLAGS) -I$(SRCDIR) -I$(SEPROOT)/include \
	-DSOURCE='"'$(SRCDIR)/sept2a.c'"' $(SRCDIR)/sept2a.c
	${MV} sept2a.o ${OBJDIR}

tridata.H edge.H: ${BINDIR}/Sept2a.x ${RAINDIR}/obsdata.HH
	< ${RAINDIR}/obsdata.HH Window | Transp > xyz.H
	< xyz.H ${BINDIR}/Sept2a.x n1=376 n2=253 \
	d1=1.009975 d2=1.009975 o1=-185.556375 o2=-127.261523 \
	> tridata.H edgeout=edge.H

${RESDIR}/trinter.v: tridata.H edge.H
	<edge.H Graph plotcol=1 title=" " wantaxis=n wantlabel=n \
	${MM} > edge.V
	< tridata.H Grey title=Triangulation transp=n yreverse=n \
	wheretitle=t wherexlabel=b color=I > data.V
	Vppen < data.V edge.V  \
	erase=once vpstyle=n big=n >/dev/null out=$@

precinter%.V: precinter%.H
	< precinter$*.H Grey wantscalebar=y yreverse=n transp=n \
	color=I pclip=100 title="Preconditioning $*" >$@

precinter%.v: precinter%.H
	< precinter$*.H Grey wantscalebar=y yreverse=n transp=n allpos=y \
	color=I pclip=100 title="Preconditioning $*" >/dev/null out=grey.v
	< precinter$*.H Contour nc=4 c0=50 dc=100 yreverse=n transp=n \
	wanttitle=n plotcol=7 wantframe=y wantaxis=n \
	screenwd=12.17 >/dev/null out=cont.v
	vp_Overlay grey.v cont.v > $@

#
#precinter%.v: precinter%.H
#	< precinter$*.H Grey wantscalebar=y yreverse=n transp=n \
#	color=j pclip=100 wheretitle=t wherexlabel=b \
#	title="Preconditioning: $* iterations" >/dev/null out=$@

${RESDIR}/precinter.v: precinter10.v precinter100.v
	vppen < precinter10.v precinter100.v \
	vpstyle=n size=r gridnum=2,1 >$@

clean: jclean
	@-${TOUCH} junk.P
	@-${RM}       *.P
	@-cd miss; $(MAKE) clean;
	@-cd seab; $(MAKE) clean;
	@-${RM} ${BINDIR}/*.x
	@-${RM} ${OBJDIR}/*.o

include ${SEPINC}/SEP.bottom
