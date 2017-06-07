from rsf.proj import *
import math

def vc3d(fft,out,vcname,cubesname,
          v0,
          nv,
          dv,
          eps):

	for i in range(nv):#i=0,1,...,99

    		vel=(v0+((i+1)*dv))
    		w1=1/(vel*vel)	
    		w12=0
    		w2=1/(vel*vel)
		
		# velocity continuation shift		

		shift='input*exp(I*%g*(x3*x3*%g-2*x2*x3*%g+x2*x2*%g)/(16.0*((x1/2.0)+%g)*(%g*%g-%g*%g)))'%(math.pi,w1,w12,w2,eps,w12,w12,w1,w2)

		# creating constant migration velocity image

		vc=vcname+'_%d'%(i)

		Flow(vc,fft,
			'''
			math output="%s" |
			fft3 axis=3 inv=y |
			fft3 axis=2 inv=y |
			fft1 inv=y |
			t2warp inv=y |
			spray axis=2 n=1 
			''' % (shift))
		
		cubesname.append(vc)

	Flow(out,cubesname,'rcat axis=2 ${SOURCES[1:%d]} | put o2=%g d2=%g label2="Velocity" unit2="km/s"'%(nv,v0+dv,dv))

def pi3d(epspi,erfi_input,fft,name,v_a,v_b):

	# 3D PI precompute

	eps_4pi=epspi

	const = '((-1)*sqrt(2*(x1+%g))*exp(I*0.25*%g)/( sqrt((x2+%g)*(x2+%g) + (x3+%g)*(x3+%g)) ))'%(eps_4pi,math.pi,eps_4pi,eps_4pi,eps_4pi,eps_4pi)

	erfi_a = '''
	         sfmath output="input*I*exp(I*0.75*%g)*%g*sqrt( 2*%g*((x2+%g)*(x2+%g)  +  (x3+%g)*(x3+%g)) )/(4*(sqrt(x1+%g)))" |
	         cerf | sfmath output="input*(-1)*I" | sfmath output="input*%s"
	         '''%(math.pi,v_a,math.pi,eps_4pi,eps_4pi,eps_4pi,eps_4pi,eps_4pi,const)

	erfi_b = '''
	         sfmath output="input*I*exp(I*0.75*%g)*%g*sqrt( 2*%g*((x2+%g)*(x2+%g)  +  (x3+%g)*(x3+%g)) )/(4*(sqrt(x1+%g)))" |
	         cerf | sfmath output="input*(-1)*I" | sfmath output="input*%s"
	         '''%(math.pi,v_b,math.pi,eps_4pi,eps_4pi,eps_4pi,eps_4pi,eps_4pi,const)

	erfi_a_f = 'erfi_a_4_'+name

	erfi_b_f = 'erfi_b_4_'+name

	erfi_ba_f = 'erfi_ba_4_'+name

	pi_erfi_fft = 'pi_erfi_fft_4_'+name

	pi = 'pi_4_'+name
	
	Flow(erfi_a_f,erfi_input,erfi_a)

	Flow(erfi_b_f,erfi_input,erfi_b)

	Flow(erfi_ba_f,[erfi_b_f,erfi_a_f],'add scale=1,-1 ${SOURCES[1]}')

	Flow(pi_erfi_fft,[fft,erfi_ba_f],'math  K=${SOURCES[1]} output="input*K"')

	Flow(pi,pi_erfi_fft,
                  '''
                  fft3 axis=3 inv=y |
                  fft3 axis=2 inv=y |
                  fft1 inv=y |
                  t2warp inv=y
                  ''')

def dpi3d(epspi,erfi_input,fft,name,v_a,v_b):

	# Double path integral

	eps_d=epspi

	k_sqr_sum = '((x2+%g)*(x2+%g) + (x3+%g)*(x3+%g))'%(eps_d,eps_d,eps_d,eps_d)

	shiftdpi_b='''
           	sfmath output="input*I*4*(x1+%g)*exp(-I*%g*%g*%s*%g/( 16*(x1+%g) ) )/(%s*%g)"
           	'''%(eps_d,v_b,v_b,k_sqr_sum,2*math.pi,eps_d,k_sqr_sum,math.pi)

	shiftdpi_a='''
           	sfmath output="input*I*4*(x1+%g)*exp(-I*%g*%g*%s*%g/( 16*(x1+%g) ) )/(%s*%g)"
           	'''%(eps_d,v_a,v_a,k_sqr_sum,2*math.pi,eps_d,k_sqr_sum,math.pi)

	derfi_a = 'derfi_a_4_'+name

	derfi_b = 'derfi_b_4_'+name

	derfi_ba = 'derfi_ba_4_'+name

	dpi_fft = 'dpi_fft_4_'+name

	dpi = 'dpi_4_'+name	

	Flow(derfi_a,erfi_input,shiftdpi_a)

	Flow(derfi_b,erfi_input,shiftdpi_b)

	Flow(derfi_ba,[derfi_b,derfi_a],'add scale=1,-1 ${SOURCES[1]}')

	Flow(dpi_fft,[fft,derfi_ba],'math  K=${SOURCES[1]} output="input*K"')

	Flow(dpi,dpi_fft,
                  '''
                  fft3 axis=3 inv=y |
                  fft3 axis=2 inv=y |
                  fft1 inv=y |
                  t2warp inv=y
                  ''')






