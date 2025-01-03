# ODF Fingerprinting scripts and example

Steven Baete [steven . baete (at) nyulangone . org], Patryk Filipiak [patryk . filipiak (at) nyulangone . org]

Matlab-scripts for ODF Fingerprinting as described in the paper:

* Steven H. Baete, Martijn A. Cloos, Ying-Chia Lin, Dimitris G. Placantonakis, Timothy Shepherd, Fernando E. Boada. Fingerprinting Orientation Distribution Functions in Diffusion MRI detects smaller crossing angles. [NeuroImage, 198:231-241m 2019](https://doi.org/10.1016/j.neuroimage.2019.05.024).

PLEASE NOTE: 

* *The software available on this page is provided free of charge and comes without any warranty. CAI2R and the NYU School of Medicine do not take any liability for problems or damage of any kind resulting from the use of the files provided. Operation of the software is solely at the user's own risk. The software developments provided are not medical products and must not be used for making diagnostic decisions.*

* *The software is provided for non-commercial, academic use only. Usage or distribution of the software for commercial purpose is prohibited. All rights belong to the author (Steven Baete) and the NYU School of Medicine. If you use the software for academic work, please give credit to the author in publications and cite the related publications.*

### Description

Diffusion tractography is routinely used to study white matter architecture and brain connectivity in vivo. A key step for successful tractography of neuronal tracts is the correct identification of tract directions in each voxel. Here we present a fingerprinting-based
methodology to identify these fiber directions in Orientation Distribution Functions, dubbed ODF-Fingerprinting (ODF-FP). In ODF-FP, fiber configurations are selected based on the similarity between measured ODFs and elements in a pre-computed library. In noisy ODFs, the library matching algorithm penalizes the more complex fiber configurations. 

The ODF-FP approach improves the detection of fiber pairs with small crossing angles while maintaining fiber direction precision, which leads to better tractography results. Rather than focusing on the ODF maxima, the ODF-FP approach uses the whole ODF shape to infer fiber directions to improve the detection of fiber bundles with small crossing angle. The resulting fiber directions aid tractography algorithms in accurately displaying neuronal tracts and calculating brain connectivity.

### Prerequisites for use of these scripts
* Matlab
* [MRtrix3](http://www.mrtrix.org). Some of the Matlab-scripts push commands to the system to MRtrix3. Make sure MRtrix3 has been set up and is on your path.
* [DSIStudio](http://dsi-studio.labsolver.org/dsi-studio-download) For display and tractography.

## External software packages used
* The package [Tools for NIfTI and ANALYZE image](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) as retrieved from [Mathworks.com](http://www.mathworks.com) in March 2018 and optimized for linux.
* Selected files from the [MIRT](http://web.eecs.umich.edu/~fessler/irt/irt/) (Michigan Image Reconstruction Toolbox) package by [Jeff Fessler](https://web.eecs.umich.edu/~fessler/) are included in this repository.

### Obtaining the code
Download from [cai2r.net -> Resources -> Software Downloads](https://www.cai2r.net/resources/software/odf-fingerprinting) or clone the repository from [bitbucket](https://bitbucket.org/sbaete/odffingerprinting)

        git clone https://bitbucket.org/sbaete/odffingerprinting.git

### Running demo-script
An example processing script can be found in *odffp_demo.m*. In this script the ODF Fingerprinting algorithm is used to identify fiber directions in simulated phantom data.

* Obtaining the simulated Phantomas phantom data: see *demodata/readme.md*. Basically, the simulated data can be easily generated with the [*rdsi_recon*-package](https://bitbucket.org/sbaete/rdsi_recon).
* Fill in the filenames of the *.src.gz* and *.fib.gz* phantom data in the *simulation Phantomas data*-section of the script *odffp_demo.m*.
* Run the script *odffp_demo.m*.
* You can open the resulting *.fp7.fib.gz* file with DSIStudio.

![Phantomas ODF-Fingerprinting tractography](https://bitbucket.org/sbaete/odffingerprinting/raw/65c83e4961aa6c24189774ff1545f92f610c718b/demodata/dwis.src.gz.fp7.fib.tracto.jpg)

## Questions, remarks and/or bugs
	
Steven Baete [steven . baete (at) nyulangone . org]
Patryk Filipiak [patryk . filipiak (at) nyulangone . org]

Department of Radiology, Center for Biomedical Imaging,  
[Center for Advanced Imaging Innovation and Research](https://www.cai2r.net),   
New York University School Of Medicine, New York, NY, USA.