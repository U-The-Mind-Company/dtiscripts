%% Steven Baete
%  NYU SOM CBI
%  March 2018

% Demonstration of ODF Fingerprinting on simulated Phantomas DWI data.

% -----------------------------------------------------------------------
% This method has been published in:
%  Baete, Cloos, Lin, Placantonakis, Shepherd, Boada. Fingerprinting Orientation
%     Distribution Functions in Diffusion MRI detects smaller
%     crossing angles. Submitted to NeuroImage, 2019.
% -----------------------------------------------------------------------

clear all;
close all;

restoredefaultpath
path(path,'utils/');
path(path,'utils/NIfTI_20140122/');

%% simulation Phantomas data
%  simulation data can be easily generated with the simulated DWI MRI data 
%  and the script rdsi_phantomas_example.m in the rdsi_recon repository at
%       https://sbaete@bitbucket.org/sbaete/rdsi_recon.git
%   or  cai2r.net -> Resources -> Software Downloads

%  When done, copy rdsi_recon/demodata/dwis.src.gz.fib.gz to the subfolder
%  demodata/ of this folder.
srcfile = 'demodata/dwis.src.gz';
fibfile = 'demodata/dwis.src.gz.rdsi.f3.fib.gz';

%% create library

generate_FP_library

display(['generated ODF FP library ' libname]);

%% identify fiber directions with ODF FP
tic;
opt.libname = libname;
opt.matchmethod = 7;
opt.srcfile = srcfile;
method = 'FP';
fibfileout = run_on_fib(fibfile,method,opt);

display([' ran FP with noise penalty  in  ' num2str(toc) ' s']);
display([' result saved in ' fibfileout]);
display([' You can open this file with DSIStudio.']);
