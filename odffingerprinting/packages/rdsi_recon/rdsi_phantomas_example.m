%% Steven Baete
%  NYU LMC CBI
%  March 2018

% Example appplication of the RDSI reconstruction on a simulated Phantomas
% RDSI dataset.

% -----------------------------------------------------------------------
% This method has been published in:
%  Baete, Yutzy and Boada. Radial q-space sampling for DSI. 
%       Magn Reson Med, 76(3):769-80, 2016
%  http://onlinelibrary.wiley.com/doi/10.1002/mrm.25917/abstract
%
%  and expanded on in:
%  Baete and Boada. Accelerated Radial Diffusion Spectrum Imaging using 
%       a Multi-Echo Stimulated Echo diffusion sequence. 
%	Magn. Reson. Med., 79(1):306-316, 2018.
%  http://onlinelibrary.wiley.com/doi/10.1002/mrm.26682/abstract
% -----------------------------------------------------------------------

close all;
clear all;

[fold,~,~] = fileparts(mfilename('fullpath'));fold = [fold filesep];
path(path,fold);
path(path,[fold filesep 'utils/']);
path(path,[fold filesep 'spams-matlab-v2.6/build/']);
cd(fold);

%% example usage of rdsi_recon

% command to call dsi_studio
dsistudiocmd ='/usr/local/dsistudio/dsi_studio_20170908/dsi_studio';
if (~(exist(dsistudiocmd,'file')==2))
    error('Could not find DSIStudio. Please adjust the variable "dsistudiocmd" to point to the executable of DSIStudio (http://dsi-studio.labsolver.org/).');
end;

% folder to save output
outfold = [fold 'demodata/'];

%% preprocess the diffusion data to a nii-file
%  dsistudio needs to be installed
inputnii = [outfold 'dwis.nii.gz'];
bvalfile = [outfold 'RDSI_withb0_252dirs.bval'];
bvecfile = [outfold 'RDSI_withb0_252dirs.bvec'];

%% create a dsistudio src-file

srcfile = [outfold filesep 'dwis.src.gz'];
if (~exist(srcfile,'file'))
    system([dsistudiocmd ' --action=src ' ...
            ' --source=' inputnii ...
            ' --bval=' bvalfile ...
            ' --bvec=' bvecfile ...
            ' --output=' srcfile]);
end;

%% run RDSI-reconstruction
FIB_FILE_OUT = rdsi_recon(srcfile,'phantommask',true,'flipBaxis',3,'phantom',true);

%% export data from the fib-file
% let's export the qa1, qa2 and qa3-maps
system([dsistudiocmd ' --action=exp --source=' FIB_FILE_OUT ...
            ' --export=fa0,fa1,fa2']);
