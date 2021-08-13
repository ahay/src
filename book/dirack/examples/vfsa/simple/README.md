# Usage examples of sfvfsacrsnh and sfnhcrssurf

In this example we use a gaussian reflector model as input and kirchhoff modeling to obtain the
seismic data cube (Seismic data organized in CMP x Offset x Time coordinates) from it.

After that we use _sfvfsacrsnh_ to fit the non-hyperbolic CRS traveltime approximation in the reflection traveltime surface
exatracted from the data cube. Finally, we plot these surfaces and approximation error side by side such as in the [Fomel's original experiment](http://www.reproducibility.org/RSF/book/tccs/crs/paper_html/).

We use the _sfnhcrssurf_ to build the non-hyperbolic CRS traveltime approximation with RN, RNIP and BETA parameters as input.

## Expected outcomes

The CRS parameters will appear on your terminal screen. They are saved in the _crsParameters.rsf_ file.
You can see them running:
 
 ```sh
 ~$ sfdisfil < crsParameters.rsf
0:          3.23       0.9706    -0.003872       0.7791          0.5
5:            10          1.1            5
```

The first 3 numbers are RN, RNIP and BETA. The fourth one is the best semblance between the non-hyperbolic CRS approximation and data using this CRS parameters (it may change from one running to another). Finally, the fifth and sixth parameters are C0 and Temp0, arbitrary parameters used in vfsa that may change from problem to problem. t0 (Normal ray traveltime) and m0 (Central CMP of the approximation) are the last ones.

## The VFSA package

The VFSA programs use Very Fast Simulated Annealing (VFSA) global optimization to obtain the zero-offset
Common Reflection Surface (CRS) parameters (RN, RNIP, and BETA) from a reflection data cube (seismic data organized in CMP x half-offset X Time coordinates).
These parameters are obtained by fiting the non-hyperbolic CRS traveltime surface in the reflection data cube using the semblance between them
as the optimization criterion. The CRS paramters (RN, RNIP and BETA) that produced the minimum error will be the optimized ones. 

- For more details please check out ["What is VFSA?"](https://github.com/Dirack/vfsa/wiki/Very-Fast-Simulated-Anneling-(VFSA)) and ["What is non-hyperbolic CRS?"](https://github.com/Dirack/vfsa/wiki/Non-hyperbolic-CRS) in our Wiki.
