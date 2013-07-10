/** 
@page asgnotes Time Step Selection for Acoustic Staggered Grid Schemes

For lowest order staggered grid schemes, the CFL condition
<p>
dt < dxmin/(cmax*sqrt(dim))
<p>
(dxmin here is the min spatial step, cmax the max velocity, dim the spatial dimension) is <i>not</i>
sufficient to guarantee stability.  In fact stabiliy depends also on
the density field, markedly so when the density has
discontinuities. Gustafsson and Wahlund, SIAM J. Sci. Comp. 26 (2004),
pp.  272-293 derived a necessary and sufficient stability condition. 

IWAVE offers the user three ways to define the time step:
<ol>
<li>[Default] Use a prescribed fraction of the the maximum stable time step, computed from the Gustaffson criterion. This option is <i>adaptive</i>: it depends on the actual velocity and density (or equivalent) fields input to the simulation. The max stable step for the (2,2) scheme is scaled by a <i>CFL fraction</i> less than 1, to accommodate the stricter stability limits of higher-order schemes. 
<ul>
<li> cfl    = [float] [Default = CFL_DEF]</li>
</ul>
The default value of this CFL fraction (CFL_DEF, given in asg/include/defaults.h) is small enough to produce a stable step for all schemes implemented in IWAVE/asg. The user may optionally input a CFL fraction (parameter cfl), which must be less than or equal to CFL_DEF. Since higher order schemes have better dispersion characteristics with time steps well below the stability limit, specifying cfl to be well below CFL_DEF is usually a good idea. We have found that cfl=0.4 or cfl=0.5 give reasonably non-dispersive results for FD orders up to 10, for propagation distances of 50-100 wavelengths. </li>
<p>

<li>Compute the time step using the asserted max velocity (parameter cmax) and the standard CFL criterion (above), scaled by a prescribed fraction (parameter cfl). 
<ul>
<li> cmax   = [float] [Default = CMAX_DEF]</li>
<li> cfl    = [float] [Default = CFL_DEF]</li>
</ul>
(see \ref alternate for more on cmax).
This option is static (non-adaptive) and implicitly defines a feasible set of models on a family of spatial grids for which the chosen time-space step relation is stable. The step is checked against the Gustaffson step for the default CFL fraction, and rejected (runtime error) if it exceeds this maximum stable step. IWAVE provides this second option to enable simple accuracy tests for inversion applications based on IWAVE, with uniform time step over a class of models. For simulation using the IWAVE package itself, the first option is usually more efficient.
<p>
To enable this second (static) method for selecting time step, add the line
<p>
max_step=0
<p> 
to the parameter file.
</li>
<li>The user can specify the time step directly (parameter dt). <i> in this case, no stability checking is performed</i>. To enable this use case, add the line
<ul>
<li>dt   = [float] [Default: computed from hdrfile, cfl condition]</li>
</ul>
It is also possible to specify the total number of time steps, <i>provided</i> that dt has been set as just described (otherwise - an error!). To add this specification,
<ul>
<li>nt   = [int] [Default: computed from hdrfile, cfl condition]</li>
</ul>
If only dt is defined as described here, then nt is computed as in other use cases to approximate the time duration of the target traces - it is then substituted for the number of samples (SEGY keyword ns), as interpolation makes no sense in this case.
<p>
If both dt and nt are provided inthe parameter file, they override <i>both</i> the number of time samples 
 and time step computed from input data and stability criteria, for <i>both</i>
 internal simulator time stepping <i>and</i> output traces - i.e.
 if this option is chosen, then output traces will have
 nt samples and time step dt, overriding the values from  
 the input header file (though the source and receiver locations
 remain as determined by the input header file). Note that 
   either dt, or both nt and dt,  must be set - setting nt alone, without dt, is prohibitied. This option is 
 useful mostly for convergence studies, for precise
 control of time step. 
<p>
 In either case, if dt is provided, then <i>no</i> stability checking is
 performed by the code in this case, so ensuring stability
 is entirely up to the user - caveat emptor! 

<p>
dt=[specified time step]
<p> 
to the parameter file.
</li>
</ol>
*/
