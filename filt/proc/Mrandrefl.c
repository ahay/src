program Randrefl
  use randn
  use cumsum
  use synf
  use quicksort
  use sep

  implicit none
  integer                          :: nr, nt, it
  real                             :: t0, dt
  real, dimension (3)              :: fo
  real, dimension (:), allocatable :: t, pp, ps, ss, tpp, tss, tps
  real, dimension (:), allocatable :: ts, dtpp, tm, p2ss, p2ps, dtss, dtps, rs

  call sep_init ()
  call from_par ("nr",nr)
  call from_par ("n1",nt,3501)
  call from_par ("d1",dt,0.001)
  call from_par ("o1",t0,0.)
  call to_history ("n2",3)
  call to_history ("esize",4)
  call from_par ("fo",fo, (/20.,8.,5./))
  call to_history ("n1",nr,"vpvs")
  call to_history ("n2",2,"vpvs")
  call to_history ("esize",4,"vpvs")
  call sep_close ()

  allocate (t (nt), pp (nt), ps (nt), ss (nt), ts (nr))
  t = (/ (t0 + (it-1)*dt, it=1,nt) /)

  allocate (dtpp (nr), tm (nr), rs (nr), tpp (nr), tss (nr), tps (nr))
  allocate (p2ss (nr), p2ps (nr), dtss (nr), dtps (nr))

  ! ts - reflector positions
  call random_number (ts)
  ts = 0.1+0.9*ts
  call qsort_init (ts)
  call qsort ()

  ! dtpp - layer thickness in PP time
  dtpp(1) = ts(1)
  dtpp(2:nr) = ts(2:nr) - ts(1:nr-1)

  ! tm - time in the middle of the layer
  tm=ts-0.5*dtpp
  ! p2ss - Vp/Vs ratio as a function of time
  p2ss=1./(0.12+(0.5*tm))

  call sep_write (tm,"vpvs")
  call sep_write (p2ss,"vpvs")

  p2ps=0.5*(1+p2ss)

  dtss=p2ss*dtpp ! SS thickness
  dtps=p2ps*dtpp ! PS thickness

  ! rs - reflection coefficient
  call randn_number (nr, rs)
  rs = 0.1/nr + 0.05*rs

  call cum_sum(dtpp,tpp)
  call cum_sum(dtps,tps)
  call cum_sum(dtss,tss)

  call synf1(pp,t,fo(1),tpp,rs)
  call synf1(ps,t,fo(2),tps,rs)
  call synf1(ss,t,fo(3),tss,rs)

  call sep_write (pp)
  call sep_write (ps)
  call sep_write (ss)

  call exit (0)
end program Randrefl



