      program Test
      implicit none
      integer sf_gettype
      integer type, n1, i
      logical sf_histint
      integer*8 in, out, sf_input, sf_output
      character*100 sf_histstring, label1
      real trace(100)

      call sf_init()
      in = sf_input("in")
      out = sf_output("out")
      type = sf_gettype(in)
      call sf_putint(out,"n2",5)
      call sf_floatread(trace,100,in)
      do 10 i=1,5
         call sf_floatwrite(trace,100,out)
 10   continue
      stop 
      end

C	$Id: Testfile.f 982 2005-01-30 23:38:22Z shan $	
