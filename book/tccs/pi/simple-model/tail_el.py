from rsf.proj import *

def tail_el(path_int,happy,sad,
          nt,
          nsh,
          lrect1,
          lrect2,
          lniter):

    shifts_h = 'shifts_h_nsh%d'%nsh
    shifts_h = [happy]
    shifts_h_cat = happy+'-shifts'+'_nsh%d'%nsh
    for s in range(1,nsh):
        shift_h = 'shift_h-%d_nsh%d' % (s,nsh)
        Flow(shift_h,happy,'window f1=%d | pad end1=%d' % (s,s))
        shifts_h.append(shift_h)

        shift_h = 'shift_h+%d_nsh%d' % (s,nsh)
        Flow(shift_h,happy,'window n1=%d | pad beg1=%d' % (nt-s,s))
        shifts_h.append(shift_h)
    Flow(shifts_h_cat,shifts_h,'cat ${SOURCES[1:%d]}' % len(shifts_h))
    #Result('shifts_h_nsh%d'%nsh,shifts_h_cat,'grey pclip=100 gainpanel=each title="shifts happy"')
    print shifts_h
    print shifts_h_cat

    flt_h='flt_h_nsh%d'%nsh
    pre_h='pre_h_nsh%d'%nsh
    sig_h='sig_h_nsh%d'%nsh
    Flow([flt_h,pre_h],[shifts_h_cat,path_int],
                 '''
                 lpf match=${SOURCES[1]} pred=${TARGETS[1]} rect1=%d rect2=%d niter=%d
                 '''%(lrect1,lrect2,lniter))
    Flow(sig_h,[path_int,pre_h],'add scale=1,-1 ${SOURCES[1]}')

    #Result(pre_h,'grey gainpanel=each pclip=100 title="happy prediction"')
    #Result(sig_h,'grey gainpanel=each pclip=100 title="- happy pre"')

    shifts_s = 'shifts_s_nsh%d'%nsh
    shifts_s = [sad]
    shifts_s_cat = sad+'-shifts'+'_nsh%d'%nsh
    for s in range(1,nsh):
        shift_s = 'shift_s-%d_nsh%d' % (s,nsh)
        Flow(shift_s,sad,'window f1=%d | pad end1=%d' % (s,s))
        shifts_s.append(shift_s)

        shift_s = 'shift_s+%d_nsh%d' % (s,nsh)
        Flow(shift_s,sad,'window n1=%d | pad beg1=%d' % (nt-s,s))
        shifts_s.append(shift_s)
    Flow(shifts_s_cat,shifts_s,'cat ${SOURCES[1:%d]}' % len(shifts_s))
    #Result('shifts_s_nsh%d'%nsh,shifts_s_cat,'grey pclip=100 gainpanel=each title="shifts sad"')

    flt_s='flt_s_nsh%d'%nsh
    pre_s='pre_s_nsh%d'%nsh
    sig_s='sig_s_nsh%d'%nsh
    Flow([flt_s,pre_s],[shifts_s_cat,path_int],
                 '''
                 lpf match=${SOURCES[1]} pred=${TARGETS[1]} rect1=%d rect2=%d niter=%d
                 '''%(lrect1,lrect2,lniter))
    Flow(sig_s,[path_int,pre_s],'add scale=1,-1 ${SOURCES[1]}')

    #Result(pre_s,'grey gainpanel=each pclip=100 title="sad prediction"')
    #Result(sig_s,'grey gainpanel=each pclip=100 title="- sad prediction"')

    tails='tails_nsh%d'%nsh
    Flow(tails,[sig_s,sig_h],'add scale=1,1 ${SOURCES[1]}')
    #Result(tails,'grey gainpanel=each pclip=100 title="tails only"')
    
    path_int_te='path_int_te_nsh%d'%nsh
    Flow(path_int_te,[path_int, tails],'add scale=1,-1 ${SOURCES[1]}')
    #Result(path_int_te,'grey gainpanel=each pclip=100 title="path - tails"')

#    flt_h='flt_h_nsh%d'%nsh
#    flt_s='flt_s_nsh%d'%nsh
#    flt='flt_nsh%d'%nsh
#    Plot(flt_h,'window n1=1 n2=1 min1=1.5 f2=50 | dots')
#    Plot(flt_s,'window n1=1 n2=1 min1=1.5 f2=50 | dots')
#    Result(flt,[flt_h,flt_s],'OverUnderIso')






