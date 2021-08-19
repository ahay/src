# CRE Stacking package

This is a package with several C programs and SConstructs adapted to the Madagascar seismic processing package.
The mains objective of those programs is to obtain the zero-offset section using CRE stacking. The main advantage of the CRE
stacking is that allows to get the macrovelocity model and stacked section without the NMO velocity analisys.
And it can be used in a velocity inversion process, such as tomography algorithms, in order to get the depth velocity model.

The CRE stack process is done defining the seismic traces that belong to the CRE trajectories in a seismic
data cube (seismic data organized in CMP x Offset X Time coordinates) for each (m0, t0) pairs in the stacked section.
The set of seismic traces that belong to a CRE trajectory form a CRE Gather. The stacking throughout this CRE gathers using
the CRE traveltime approximation and the assignment of this amplitude summation to the (m0,t0) pairs in the stacked section
is the core of the CRE stacking process.

* For more theoretical details of CRE Gather stack please download and read this article:
  [The common reflecting element CRE method revisited (2000)](https://github.com/Dirack/creGatherStack/files/5213160/The_common_reflecting_element_CRE_method_revisited_cruz2000.pdf)
  
## CRE stack uses Very Fast Simulated Anneling (VFSA) global optimization

As explained above, a set of seismic traces that belong to a specific CRE trajectory form a CRE Gather.
Those trajectories are defined for each (m0,t0) pairs in the stacked section and for a given RNIP and BETA parameters
from the zero offset Common Reflection Surface (CRS) process. These parameters are extracted from the seismic data
cube using Very Fast Simulated Anneling (VFSA) global optimization algorithm.

For more details about VFSA, please look at the [VFSA package repository](https://github.com/Dirack/vfsa)

## CRE stack uses Predictive Adaptative Error Filters (PEF) interpolation:

CRE stack also needs increased CMP data sampling through Predictive Adaptative Error Filters (PEF) interpolation.
This happen because the CRE traveltime aproximation is derived from geometric considerations in a constant velocity model in the 
neighborhood of a normal ray and the sources (s_i) and receivers (r_i) in this geometry are distributed assimetrically along the aquisition
surface because of the reflector's curvature. Besides that, the reflection rays in the CRE gather _have the same reflection point_ and
are associated with the same normal ray defined by RNIP and BETA parameters (look at the Figure bellow, an schematic representation of a CRE gather geometry).
