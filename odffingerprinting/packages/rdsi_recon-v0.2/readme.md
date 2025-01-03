# Radial Diffusion Spectrum Imaging (RDSI) reconstruction

Matlab-scripts for Radial DSI reconstruction as described in the paper:

* Steven Baete, Stephen Yutzy and Fernando, E. Boada. Radial q-space sampling for DSI. [Magn Reson Med, 76(3):769-80, 2016](http://onlinelibrary.wiley.com/doi/10.1002/mrm.25917/abstract).

and expanded on in the paper:

* Steven H. Baete and Fernando E. Boada. Accelerated Radial Diffusion Spectrum Imaging using a Multi-Echo Stimulated Echo diffusion sequence. [Magn. Reson. Med., 79(1):306-316, 2018](http://onlinelibrary.wiley.com/doi/10.1002/mrm.26682/abstract).

PLEASE NOTE: 

* *The software available on this page is provided free of charge and comes without any warranty. CAI²R and the NYU School of Medicine do not take any liability for problems or damage of any kind resulting from the use of the files provided. Operation of the software is solely at the user's own risk. The software developments provided are not medical products and must not be used for making diagnostic decisions.*

* *The software is provided for non-commercial, academic use only. Usage or distribution of the software for commercial purpose is prohibited. All rights belong to the author (Steven Baete) and the NYU School of Medicine. If you use the software for academic work, please give credit to the author in publications and cite the related publications.*

## Data acquisition
A gradient table for the Siemens systems is included in the file *DiffusionVectors_RDSI.dvs*. This gradient scheme has 4 shells with each 59 directions according the the paper indicated above. It should be run with a b-value of 4000 s/mm2.

## Prerequisites for use of these scripts
* Matlab
* [DSIStudio](http://dsi-studio.labsolver.org/dsi-studio-download)
* The [SPAMS](http://spams-devel.gforge.inria.fr/downloads.html) (SPArse Modeling Software) package. An older version of this package is included in this repository. The necessary MacOS (MacOS X Maverick (10.9.5) and Linux-binaries are included. You might have to recompile the mex-files for your system otherwise.
	
## Processing steps

### Obtaining the code
Download from cai2r.net -> Resources -> Software Downloads or clone the repository from bitbucket

	git clone https://sbaete@bitbucket.org/sbaete/rdsi_recon.git

### Example processing
An example processing script can be found in *rdsi_example.m*. 

It is assumed that the diffusion data is preprocessed according to individual taste (denoising, gibbs ringing correction, eddy current correction, motion correction, etc) to a nii.gz-file named in the *inputnii*-variable and that the *inputnii* is accompanied by an FSL-style *data.bvec* and *data.bval*-file.

	inputnii = ['/Exp/rdsi_scripts/data.nii.gz'];
	
	ls /Exp/rdsi_scripts/data/
	data.bval  data.bvec  data.nii.gz
	
In *rdsi_example.m*, also make sure *dsistudiocmd*-variable points to the location where dsistudio is installed. In addition, set the the folder where you would like the output to be generated in the *outfold*-variable.

	dsistudiocmd ='/usr/local/dsistudio/dsi_studio_20170908/dsi_studio';
	outfold = '/Exp/rdsi_scripts/data/';

The script will automatically generate a dsistudio *src.gz*-file for further reconstruction. If this would not work on your system, you can perform this step directly using DSIStudio.

	/Exp/rdsi_scripts/data/data.src.gz

The actual RDSI reconstruction happens in *rdsi_recon.m*. The resulting *fib.gz*-file can be opened with DSIStudio for tractography. 

	/Exp/rdsi_scripts/data/data.src.gz.fib.gz

QA-maps of the first, second and third identified direction will be exported by the script to nifti-files with respective suffixes *.fa0.nii.gz*, *.fa1.nii.gz* and *.fa2.nii.gz*.

	QA1 -> /Exp/rdsi_scripts/data/data.src.gz.fib.gz.fa0.nii.gz
	QA2 -> /Exp/rdsi_scripts/data/data.src.gz.fib.gz.fa1.nii.gz
	QA3 -> /Exp/rdsi_scripts/data/data.src.gz.fib.gz.fa2.nii.gz
  
   
  
   	
	
Steven Baete and Fernando Boada

New York University School Of Medicine  
Department of Radiology, Center for Biomedical Imaging  
[Center for Advance Imaging Innovation and Research](cai2r.net)
