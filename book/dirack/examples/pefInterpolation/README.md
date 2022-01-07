#### Usage example of pefInterpolation.py recipe to double CMP data sampling

This Sconstruct shows how to use functions defined in pefInterpolation.py recipe.

As initial step, to build the input of pefInterpolation function (stacked section), it uses 
[kimodel.py](https://github.com/ahay/src/blob/master/book/Recipes/kimodel.py) and
[velocityAnalysis.py](https://github.com/ahay/src/blob/master/book/Recipes/velocityAnalysis.py) recipes to 
generates a multi layer model with 3 layers and 2 interfaces, apply Kirchhoff-Newton modeling
to generate a seismic data cube (Time x Offset x CMP data) and to do
velocity analysis with automatic velocity picking to generate the stacked section.

So, it generates another stacked section with doubled number of CMP samples through Predictive Adaptative Error Filters (PEF) interpolation
in two steps (it is done inside the function): First, it interleaves the original stacked section with zeroed traces and it calculates the
PEF coeficients using the program _sfapef_. Second, it interpolates zeroed traces using original traces in the neighborhood as a priori information.
Interpolation is done with the program _sfmiss4_.

```py
# Transpose to (time x offset x CMP)
Flow('transposedSection','stackedSection','transp plane=23')

# PEF coeficients
a1=10
a2=4

# PEF smooth radius
rect1=10
rect2=5

# PEF interpolation of zero offset section
pef('transposedSection',
	'section',
	nm=201,
	dm=0.025,
	nt=1001,
	dt=0.004,
	nhi=1,
	a1=a1,
	a2=a2,
	rect1=rect1,
	rect2=rect2
)
```

The same process is done with the first 3 offset gathers in the seismic data cube, stored in the file 'multiLayerDataCube.rsf' (data cube is the seismic data
organized in time x Offset x CMP coordinates). The number of offset gathers to interpolate extracted from the data cube is controled by _nhi_ parameter
passed to the function pefInterpolation (in this example called as pef function).

```py
# PEF interpolation of some offset gathers
# of the data cube (time x offset x CMP)
pef('multiLayerDataCube',
	'interpolatedDataCube',
	nm=201,
	dm=0.025,
	nt=1001,
	dt=0.004,
	nhi=3,
	a1=a1,
	a2=a2,
	rect1=rect1,
	rect2=rect2
)
```

The parameters _a1_ and _a2_ are the number of PEF coeficients to calculate for each sample, in time and space, respectively.
The _rect1_ and _rect2_ parameters are the smoothing radius in time and space, respectively. Those parameters were chosen through
trial and error for this problem (the parameters that build the better interpolated section were chosen).
