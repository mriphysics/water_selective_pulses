## water_selective_pulses
MRI with parallel transmit for optimized water selective imaging.

The code provided here performs an example calculation of optimized spectral-spatial water selective excitation pulses using parallel transmission. An example dataset from a 3T MRI system with 8 channels is provided for testing. The method is fully described in [this paper](http://dx.doi.org/10.1002/mrm.22260). You are free to use this code under the terms of the MIT license, but please cite [this paper](http://dx.doi.org/10.1002/mrm.22260) in any subsequent publications.

### Usage
Run the script runme.m to generate an example for a 4 element spectral-spatial pulse, and a 1-3-3-1 binomial pulse for comparison.

### /lib
The library contains code for magnitude least squares optimization, for generating the system matrix, and for generating pulse and gradient waveforms that could be used on an MRI scanner. 

The functions for generating the RF and gradient waveforms are old and not thoroughly tested and are provided purely for demonstrative purposes. 