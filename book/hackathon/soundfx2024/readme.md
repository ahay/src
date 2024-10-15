<img src="https://github.com/sfomel/JSG-Open-Source-Hackathon-2024/blob/main/logo.png" style="background-color:white;">

# Project directory for the SoundFX team

<img src="SoundFX-logo.png" style="background-color:white;">

## Overview

The primary objective of our project is to reproduce computations from the influential geophysics paper "Trace Interpolation Using the F-X Domain" by Simon Spitz (1991). Trace interpolation is a crucial step in enhancing the processing of seismic data. Traditional interpolation methods, such as spline or sinc functions, often struggle to prevent spatial aliasing during the interpolation process. Spitz proposed a novel approach utilizing the **Frequency-Time (F-X)** domain, allowing for interpolating regularly sampled data without requiring estimation or prior knowledge of event dips. This method estimates an interpolation operator using a spatial prediction filter derived from the recorded traces.

Additionally, we have decided to reproduce examples from the paper "Comparisons of Interpolation Methods" by Ray Abma and Nurul Kabir (2005). 

Our project aims to reproduce the methodologies outlined in these papers and provides a comprehensive resource for researchers and practitioners interested in seismic data interpolation techniques. Integrating our project into the Madagascar software package should allow other researchers to verify our experiments and build upon these results.

##  Key milestones

- [x] **Synthetic Data Reproduction**: All synthetic results from Spitz's 1991 paper are meticulously recreated, allowing for validation and exploration of the F-X interpolation method.
- [X] **Application on real data and comparison with alternate methods**: Reproduced examples from Abma and Kabir's 2005 paper, which provide insights into the performance of  interpolation methods and comparison on real data examples, bridging the gap between theory and field application.
- [X] **Application of machine learning**: The trained fully connected Artificial Neural Network (ANN) can reconstruct the missing traces in our real data examples.

## Project structure

### Synthetic data results [synthspitz](synthspitz/SConstruct)

During this hackathon, our team successfully reproduced all synthetic results presented in Spitz's original paper. 

### Synthetic data results [synthabma](abma/SConstruct)

Synthetic examples from Abma and Kabir (2005) evaluate various interpolation techniques, including F-X interpolation (as proposed by Spitz, 1991) and T-X interpolation (introduced by Clearbout, 1992). The nonstationary T-X interpolation is reproduced from [Liu and Fomel (2011)](../../tccs/apefint/ray/SConstruct).

### Real data results [real](abma/SConstruct)   

The comparisons extend to a field dataset from the North Sea.

### Machine Learning  [SFX_ML.ipynb](SFX_ML.ipynb)

In addition to traditional interpolation methods, our team has implemented a straightforward machine-learning technique for trace interpolation. We utilized a fully connected network to predict the missing seismic traces by training the prediction on a decimated version of the input data in the FX domain.
