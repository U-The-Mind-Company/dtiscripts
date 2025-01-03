%% Steven Baete
%  NYU LMC CBI
%  September 2017

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

[fold,~,~] = fileparts(mfilename('fullpath'));
path(path,[fold filesep 'spams-matlab-v2.6/build/']);

%% example usage of dsistudio_recon_radial

% command to call dsi_studio
dsistudiocmd ='/usr/local/dsistudio/dsi_studio_20170908/dsi_studio';

% folder to save output
outfold = '/Exp/rdsi_scripts/data/';

%% preprocess the diffusion data to a nii-file
%  dsistudio needs to be installed

inputnii = [outfold 'data.nii.gz'];

%% create a dsistudio src-file

srcfile = [outfold filesep 'data.src.gz'];
if (~exist(srcfile,'file'))
    system([dsistudiocmd ' --action=src --source=' inputnii ...
            ' --output=' srcfile]);
end;

%% run RDSI-reconstruction
FIB_FILE_OUT = rdsi_recon(srcfile,'flipBaxis',2);


%% export data from the fib-file
% let's export the qa1, qa2 and qa3-maps
system([dsistudiocmd ' --action=exp --source=' FIB_FILE_OUT ...
            ' --export=fa0,fa1,fa2']);
