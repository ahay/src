## Geoscience Hackathon 2024

Project directory for SoundFX team

### Overview

The primary objective of our project is to reproduce the influential geophysics paper "Trace Interpolation using the F-X Domain" by Simon Spitz (1991). Trace interpolation is a crucial step in enhancing the migration of seismic data. Traditional interpolation methods, such as spline or sinc functions, often struggle to prevent spatial aliasing during the interpolation process. Spitz proposed a novel approach utilizing the **Frequency-Time (F-X)** domain, allowing for the interpolation of regularly sampled data without requiring prior knowledge of event dips. This method involves estimating an interpolation operator through a spatial prediction filter derived from the recorded traces.

### Project structure

During this hackathon, our team successfully reproduced all synthetic results presented in Spitz's original paper. These results can be found in the directory:

* Synthetic data results : ```./synthspitz/```

Additionally, we have reproduced examples from the paper "Comparisons of Interpolation Methods" by Ray Abma and Nurul Kabir (2005). This study evaluates various interpolation techniques, including F-X interpolation (as proposed by Spitz, 1991) and T-X interpolation (introduced by Clearbout, 1992). The comparisons extend to real data from the North Sea, which can be accessed in the directory:

* Synthetic data results : ```./synthabma/```
* Real data results      : ```./real/```

In addition to traditional interpolation methods, our team has explored the implementation of a simple machine learning technique to interpolate the trace. We have utilized fully-connected network to predict the missing seismic traces by training them on a decimated version of the real data.

###  Key milestones

* **Synthetic Data Reproduction**: All synthetic results from Spitz's 1991 paper are meticulously recreated, allowing for validation and exploration of the F-X interpolation method.
* **Application on real data and comparison with alternate methods**: Reproduced examples from Abma and Kabir's 2005 paper which provide insights into the performance of  interpolation methods and comparison on real data examples, bridging the gap between theory and field application.
* 
