/** \mainpage Time Stepping for Optimization

This package systematizes the construction of inversion via first-order optimization from a base simulator.
<p>
Main user-supplied or -decided requirements: 

<ul>
<li>  <a class="el" href="classTSOpt_1_1Time.html" title="The Base Time Class"> Time </a> type: Time provides an abstraction for storing time information.

\image html classTSOpt_1_1Time.png 

 The StdDiscreteTime subclass is the simplest specialization, storing an integer timestep index and nothing else.</li>

<li> State template parameter - required attributes are:
<ol>
<li> Time & getTime();</li>
<li> Time const & getTime() const;</li>
<li> void setTime(Time const &);
<li> void copy(State const &);
</ol>
giving access to an internally stored Time (subtype) and a copy facility. Otherwise, the internal structure of State is entirely at the user's disposal. Typically the State will encapsulate all data needed to describe the system state, and perhaps methods for timestepping as well.
</li>

<li>  <a class="el" href="classTSOpt_1_1TimeStep.html" title="The Base TimeStepping Class"> TimeStep </a>: simple interface for a \ref StateAlg, template on State type. On top of \c StateAlg's functionality, \c TimeStep also adds the functions \c setTime() and \c getTime() for reading and changing the simulation time. (Presumably, TimeStep's \c setTime() and \c getTime() functions use the corresponding methods defined in the \c State object.)

\image html classTSOpt_1_1TimeStep.png 

The user will need to provide several of these, corresponding to the base simulation, its linearization, and the adjoint linearization of the time step. The linearization should update only the perturbation fields, not the reference field. Similarly for the adjoint TimeStep. The RVLAlg::Algorithm::run() method for each of these classes implements a single time step for the relevant fields (encapsulated in the State type - the same type is presumed to be), either through delegation to update methods of the State or via independent functions or methods. The library supplies a standard construction (StdTimeStep), which may be further specialized.
</li>

<li> <a class="el" href="classTSOpt_1_1TimeTerm.html" title="The Base TimeTerm Class"> TimeTerm </a>: a subtype of Algorithm::Terminator that stores a TSOpt::Time object. Subclasses determine the return from Algorithm::Terminator::query by comparing the TSOpt::Time data member with the TSOpt::Time of a TSOpt::State object, typically a reference data member.

\image html classTSOpt_1_1TimeTerm.png

</li>

<li><a class="el" href="classTSOpt_1_1Sim.html" title="The Base Simulator Class"> Sim </a> encapsulates a time-stepping loop; it is a RVLAlg::StateAlg subtype. The RVLAlg::Algorithm::run method combines a TSOpt::TimeStep and a TSOpt::TimeTerm to create a Algorithm::LoopAlg. An additional feature is added to most subclasses, namely the ability to reconstruct TSOpt::State at an arbitrary time. This feature is essential for adjoint simulation.

\image html classTSOpt_1_1Sim.png

</ul>

The package provides a number of simple tests in \ref testsrc which illustrate the construction of various objects from library components, and which may be used to verify correct installation. \ref testsrc is set up as a regression set, so "make regress" will run all of them and output a report in regress.rpt.

<hr>
<a href="../../../doc/html/index.html">RVL Home Page</a>

*/
