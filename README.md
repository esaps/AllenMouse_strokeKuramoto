# AllenMouse_strokeKuramoto
The repository contains scripts, functions and data used for the Brain Network Model of stroke and rehabilitation in mice in the article: "Experimental and computational study on motor control and recovery after stroke: towards a constructive loop between experimental and virtual embodied neuroscience", submitted in Frontiers of Systems Neuroscience in 2019.

_Author_: Spase Petkoski 12.01.2020, 

_Contact_: spase.petkoski[at]univ-amu.fr; spase.petkoski[at]gmail.com

### Brain Network Model of stroke and rehabilitation in mice, based on the Allen Mice Connectome and Kuramoto oscillators

## Datasets

- `res25minvox1minvol1.h5` contains a connectome (tract weights and lengths, region labels and coordinate sof the centres for each region) imported using the Allen Connectome Builder in TVB https://www.thevirtualbrain.org/tvb/zwei  on 26.01.2017, with the finest possible resolution (25 microns, minimum volume = 1, and minimum number of voxels = 1). The weights are defined as the ratio of projections versus the total number of injections (see http://docs.thevirtualbrain.org/_modules/tvb/adapters/creators/allen_creator.html for more details). 
There are 540 regions, of which 86 are cortical (43 at each hemisphere).
- `regsNamesAcr.mat` contains the accronims of the brain regions from the Interactive Allen Mouse Atlas. 
- `mask24acr.mat` contains the mapping between the number of the region in the mask, and the accronim in the list of accronims.
- `Im_MASKs24AM128.mat` contains masks of each of the 24 regions in the filed of view represented by 128x128 pixles. 
 
# scripts
- simKMcnctmAM.m --> is the scripts used for simulating the stroke versus the control/healthy state and printing the relative and absolute differences in the functional connectivity at both conditions. It gives the difference in the phase coherence for each pair of regions within the field of view (24 regions) between the control/healthy case and different strengths of stroke or recovery.
It contains the parameters sweeps for the global coupling, the stroke ratio, and rebound, as well as the definitions of 10 different sets of regions, of which the fifth is used in the results shown in the manuscript. 
 All the values are the same as the ones used for the results in the manuscript, but they can be further changed. 

# functions

- KMcnctmAM --> returns the simulated phases over time, within which the stroke has occured, so that the statistics will be obtained for the periods before and after the stroke. 
   
- KMcnctmHt0 --> performs the numerical integration of a network of instanteneously couled Kuramoto oscillators

- alloc_cnctmPlosCB --> allocates the initial values for the simulation of the Kuramoto oscillators

