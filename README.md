# superkerrline
The superkerrline XSPEC model: for reflection modelling in the general Kerr metric.   

## Paper
This repository contains the XSPEC model that was developed in [Mummery and Ingram 2024](https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/stae140/7521315?utm_source=advanceaccess&utm_campaign=mnras&utm_medium=email). 

Please cite this paper if you use this model. 

## Code 
* The main program is a fortran90 file ``superkerrline.f90``
* The main program contains both of the models ``skine`` and ``skconv``
* Raytracing modules are included in the fortran90 file ``amodules.f90``, and associated ``X.mod`` files. 
* Files to load this into XSPEC ``load.xcm`` and ``lmodel.dat`` is also included. 


## Full model input parameter description

All of the model parameters are listed below, along with their units and limits of validity. Note that $`r_I`$ is the ISCO radius. 

| Parameter | Units | range | Interpretation |
 --- | --- | --- | ---
| $`a_\star `$ | Dimensionless | $`-\infty < a_\star < \infty`$ | Kerr metric spin parameter |
| $`i `$  | Degrees | $`0 < i < 90`$ | Source-observer inclination angle | 
| $`\gamma_i`$  | Dimensionless | $`0 < \gamma_i`$ | Inner emissivity profile $`\epsilon \propto r^{-\gamma_i}`$ |
| $`\gamma_o`$  | Dimensionless | $`0 < \gamma_o`$ | Outer emissivity profile $`\epsilon \propto r^{-\gamma_o}, \quad r \geq r_{\rm br}`$ |
| $`r_{\rm br}`$  | Gravitational radii  | $`r_I \leq r_{\rm br}`$ | Emissivity break radius |
| $`r_{i}`$  | Gravitational radii  | $`r_I \leq r_{i}`$ | Disc inner edge (if negative in units of ISCO) |
| $`r_{o}`$  | Gravitational radii  | $`r_I \leq r_{o}`$ | Disc outer edge |


In addition, the ``skline`` model takes as input a line energy $`E_{\rm line}`$, which has units of keV. 
The ``skconv'' model is an XSPEC convolution model, and takes as input a model spectrum.  

## Loading into XSPEC 
* For use in the XSPEC software (https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
* XSPEC requires the files: ``lmodel.dat``, ``load.xcm``, ``superkerrline.f90``, ``amodules.f90`` and all ``.mod`` files
* ``superkerrline.f90`` contains the actual fortran implementation of the model 
* See https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html for more information on how to load the model into XSPEC
* For any questions about the model, please email andrew.mummery@physics.ox.ac.uk 
