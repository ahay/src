#### Usage examples of sfvfsacrsnh for a set of (t0,m0) pairs

In this usage example we use the same gaussian reflector model as input and kirchhoff modeling. We use _sfvfsacrsnh_ to obtain the CRS parameters for a set of (t0,m0) pairs.

### How to run it?

- Run 'scons' in this directory (it may take some time):

```sh
~$ scons
```

#### Expected outcomes

Scons will generate a parameters cube in 'crsParameters.rsf' file. The optimized parameters
are organized in the current order: RN, RNIP, BETA, Semblance, C0, Temp0, t0, m0.
The RN, RNIP and BETA are the zero offset CRS parameters. Semblance is the semblance measure
between the non-hyperbolic CRS traveltime approximation and the seismic data cube in a CMP x
offset window around Central CMP m0. The parameters C0 and Temp0 are the VFSA parameters used
inside the optimization algorithm. And (t0,m0) pair is a time x CMP coordinate around the 
zero offset CRS approximation is calculated.

The parameters cube in 'crsParameters.rsf' file is organized as follows:

```sh
~$ sfin crsParameters.rsf

crsParameters.rsf:
    in="/home/dirack/rsfdata/vfsa/usage_examples/fullParametersOptimization/crsParameters.rsf@"
    esize=4 type=float form=native 
    n1=8           d1=1           o1=0          label1="parameters" unit1="" 
    n2=6           d2=1           o2=0          label2="(t0,m0) index" unit2="" 
    n3=1           d3=1           o3=0          label3="" unit3="" 
	48 elements 192 bytes
```

For each (t0,m0) pair (n2 coordinate) we get 8 parameters (n1 coordinate) as described above:

```sh
~$ sfdisfil < crsParameters.rsf col=8

   0:         2.677        2.228      -0.2981       0.4372          0.5           10          1.1            3
   8:         4.815        3.804      -0.2447        0.374          0.5           10          1.2            3
  16:         4.768        3.813      -0.2614       0.1921          0.5           10          1.1            4
  24:         3.127        1.081       -0.287       0.1964          0.5           10          1.2            4
  32:         3.105        1.129    -0.002743       0.9898          0.5           10          1.1            5
  40:          3.12        1.168     -0.01008       0.9028          0.5           10          1.2            5
```
