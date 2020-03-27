# Project Title

hs1D (hillslope 1D) solves groundwater flux and transport at the hillslope scale, assuming Boussinesq equations and taking into account seepage generations and vertical stratification. 

## Getting Started

```
https://github.com/jmarcais/hs1D_beta
```
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run hs1D, you need to have matlab installed and the image processing toolbox. For the topographic analysis and hillslope delineations, we rely heavily on Topotoolbox, a matlab toolbox for topographic analysis (Schwangart et al. 2010).
An older version of Topotoolbox is in the external tools folder of hs1D.

### Installing

To install, clone the master branch of the git repository:
```
git clone https://github.com/jmarcais/hs1D_beta.git
```

## Running the tests

You can run the example test file:
```
test_hs1D.m
```

This example file will:

1/ run the topographic analysis on a given watershed (ie a DEM with outlet coordinates).

2/ run the hillslope storage Boussinesq equations (ie groundwater flow) component that takes into account seepage and direct precipitations onto saturated areas dynamic (Marcais et al. 2017).

3/ run the Lagrangian transport equations than enables to retrieve dynamic transient transit time distributions in the river and transient residence time distributions in the aquifer (Marcais et al. in prep.).


## Authors

* **Jean Marcais** - *Initial work* - [hs1D](https://github.com/jmarcais/hs1D)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.


## Acknowledgments

I acknowledge the french Ministry for agriculture and the ANR through a MOPGA program who have funded this project.

