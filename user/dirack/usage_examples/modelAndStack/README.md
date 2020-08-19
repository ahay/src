#### Usage example of kimodel.py and velocityAnalysis.py recipes to Kirchhoff-Newton modeling and NMO stack

This Sconstruct shows how to use functions defined in kimodel.py and velocityAnalysis.py recipes. It generates
a multi layer model with 3 layers and 2 interfaces, apply Kirchhoff-Newton modeling to generate a seismic data cube
(Time x Offset x CMP data), it does velocity analysis with automatic velocity picking and generates the stacked section.

It also has predefined aliases to split the building in SConstruct.

```py
# Use aliases to split building
# Call it 'scons alias' to build
# target identified by alias
Alias('model','interfaces.rsf')
Alias('kimodel','multiLayerDataCube.rsf')
Alias('velan','stackedSection.rsf')
```

* Aliases, _model_, _kimodel_ and _velan_ can be used to do the building step-by-step. Example,
to generate the model with interfaces only:

```
~$ scons model
```

* To generate the data cube:

```shell
~$ scons kimodel
```

* To do velocity analysis, NMO correction and stack:

```sh
~$ scons velan
```
