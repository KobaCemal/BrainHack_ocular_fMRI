# BrainHack_ocular_fMRI
Repository for scripts created during BrainHack2025 for ocular fMRI project


### Data types:

1- Raw data 
2- Denoised data
3- Smoothed data

### Masks

1- Whole-brain (no mask)
2- Big mask
3- Smaller mask
4- Rectangular
5- Optic nerve mask


### Isolating Strategies:
1- Miko’s U-net
2- “Defacer” search
3- Start from the other end of the image 


### Quantifying eye movements:
1- MVPA
2- Center of gravity 
3- Double motion correction on optic nerves
4- Slice-timing trick


### Example algortihms:
Double motion correction on optic nerves comply with eye movement tasks: https://doi.org/10.1002/hbm.10070 (motion correction and smoothing)

Variance of the time series from eye orbits match with eye movement tasks: https://doi.org/10.1002/mrm.10345 (no preprocessing?)

Semi-automated continous max-flow alghoritm to isolate eye orbits and follow the movements (this is the closest algorithm to our goal):  10.1186/s12868-016-0282-7 (only motion correction and smoothing)

Independent component analysis from eye regions relate with eye-tracker parameters:  https://doi.org/10.1101/2024.10.19.619187 (only motion correction is applied)

Data from the "large mask" reveals the gaze direction:  https://doi.org/10.1093/cercor/bhz157 (this model needs calibration and training)

DL-based model that was trained on datasets containing different eye movement tasks: https://doi.org/10.1038/s41593-021-00947-w (motion correction/ field bias correction)



Paper describing the methods: https://doi.org/10.1016/j.neuron.2015.02.027

