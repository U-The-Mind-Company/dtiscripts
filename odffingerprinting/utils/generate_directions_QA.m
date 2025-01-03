%% Create a library of voxel configurations for the ODF FP library

% Steven Baete
% NYU SOM CBI
% October 2016

%   Input:
%       anglefile: tesselation file to use, (default odf8.mat)
%       max_fiber: maximum number of fibers per voxel (default 3)
%       water: water fraction in the voxels (default 0.1)
%       adclim: min and max adc values (default 1.0e-9 m^2/s)
%       microlim (n x 2): min max values for microstructure 
%               (default, DTI with FA 0.2 - 0.6)
%   Options (in 'keyword','value' -format)
%       ratiosteps: number of steps in the fiber volume fraction between
%		water fraction and 1 (default 9)
%       microsteps (n x 1): number of microstructure param steps between 
%               min and max microstructure values (default 5)
%       adcsteps: number of adc steps between min and max adc (default 1)
%       alignZ: align the first/main fiber direction with the Z-direction 
%           (default true)
%       maxmicrodiff (n x 1): maximum difference between the microstructure
%           parameters of fibers in one voxel (default 0.3)
%       maxratdiff: maximum difference between the volume fractions of the 
%           fibers (default 1.0 = no limitation)
%   Output:
%       vox: structure with the voxel configurations
%       opt: structure of configuration library options


function [vox,opt] = generate_directions_QA(anglefile,max_fiber,water,adclim,microlim,opt)

% options in opt will overwrite other inputs!

if (nargin < 1 || isempty(anglefile))  anglefile = 'odf8.mat'; end;
if (nargin < 2 || isempty(max_fiber))  max_fiber = 3; end;
if (nargin < 3 || isempty(water))  water = 0.1; end;
if (nargin < 4 || isempty(adclim))  adclim = 1.0; end; %m^2/s
if (nargin < 5 || isempty(microlim))  microlim = [0.2,0.6]; end; %DTI FA by default if not specified
if (nargin < 6 || isempty(opt)) opt = struct; end;
if (~isfield(opt,'ratiosteps'))  ratiosteps = 9; else ratiosteps = opt.ratiosteps; end; %m^2/s
if (~isfield(opt,'microsteps'))  microsteps = 5; else microsteps = opt.microsteps; end; %m^2/s
if (~isfield(opt,'adcsteps'))  adcsteps = 1; else adcsteps = opt.adcsteps; end; %m^2/s
if (~isfield(opt,'alignZ')) alignZ = true; else alignZ = opt.alignZ; end;
if (~isfield(opt,'maxmicrodiff')) maxmicrodiff = 0.3; else maxmicrodiff = opt.maxmicrodiff; end;
if (~isfield(opt,'maxratdiff')) maxratdiff = 1.0; else maxratdiff = opt.maxratdiff; end;
if (~isfield(opt,'dispersionsteps')) dispersionsteps = 1; else dispersionsteps = opt.dispersionsteps; end;
if (~isfield(opt,'maxdispersion')) maxdispersion = 0.0; else maxdispersion = opt.maxdispersion; end;
if (~isfield(opt,'maxdispersiondiff')) maxdispersiondiff = 90.0; else maxdispersiondiff = opt.maxdispersiondiff; end;
if (~isfield(opt,'microstruct')) microstruct = 'DTI'; else microstruct = opt.microstruct; end;
if (~isfield(opt,'minangle')) minangle = []; else minangle = opt.minangle; end;
% check if some parameters in the opt-structure overrule regular inputs
optvars = {'anglefile';'max_fiber';'water';'microlim';'adclim';};
for i = 1:length(optvars)
    if (isfield(opt,optvars{i}))
        eval([optvars{i} ' = opt.' optvars{i} ';']);
    end;
end;
% check # of parameters for microstructure model
display(['   generate_directions_QA: ']);
display(['      generate voxel directions for microstructure model ' microstruct]);
switch (microstruct)
    case 'DTI'
        % FA, incombination with the adc
        nmicropar = 1;
        microlimdefault = [0.2,0.6];
    case 'twocomp'
        % two compartment, see Jelescu2017 eq (14)
        % (lambda_i, lambda_e, lambda_e_p, v1)
        % the adc now only counts for free water
        nmicropar = 4;
        microlimdefault = [1.5,2.5;1.5,2.5;0.5,1.5;0.2,0.6];
end;
if (size(microlim,1) > nmicropar) microlim = microlim(:,1:nmicropar); end;
if (size(microlim,1) < nmicropar) microlim = microlimdefault; end;
if (length(microsteps) > nmicropar) microsteps = microsteps(nmicropar); end;
if (length(microsteps) < nmicropar) microsteps = microsteps(1)*ones(nmicropar,1); end;
if (length(maxmicrodiff) > nmicropar) maxmicrodiff = maxmicrodiff(nmicropar); end;
if (length(maxmicrodiff) < nmicropar) maxmicrodiff = 0.3*ones(nmicropar,1); end;

if (maxdispersion > 50), maxdispersion = 50; end;
if (length(adclim)==1) adclim = [adclim,adclim]; end;
if (size(microlim,2)==1) microlim = [microlim,microlim]; end;

%% library options

opt_adcs = adclim(1) + (adclim(2)-adclim(1))*(0:1/(adcsteps-1):1);
% number of fibers
opt_ndir = max_fiber;
% ratios
for i = 1:length(water)
    opt_ratios{i} = [];
    ratios{i} = flip(0:((1-water(i))/ratiosteps):(1-water(i)));
    opt_ratios{i} = ratio_options(ratios{i},max_fiber,true,true,maxratdiff);
end;
% microstructure parameters 
for i = 1:nmicropar
    microranges{i} = flip(microlim(i,1) + (microlim(i,2)-microlim(i,1))*(0:1/(microsteps(i)-1):1));
end;
micros = micro_options(microranges(1:end));
for i = 1:max_fiber
    opt_micros{i} = ratio_options_micro(micros,i,maxmicrodiff);
end;
% falim = [0.2,0.6];
% fasteps = 5;
% maxfadiff = 0.3;
% fas = falim(1) + (falim(2)-falim(1))*(0:1/(fasteps-1):1);
% for i = 1:max_fiber
%     opt_fas{i} = ratio_options(flip(fas),i,false,false,maxfadiff);
% end;
% dispersion
dispersions = maxdispersion*(0:1/(dispersionsteps-1):1);
for i = 1:max_fiber
    opt_dispersions{i} = ratio_options(dispersions,i,false,false,maxdispersiondiff);
end;
% directions
load(anglefile)
angles = odf_vertices(:,1:(end/2));
[azimuth,elevation,r] = cart2sph(angles(1,:),angles(2,:),angles(3,:));
opt_angles = [azimuth;elevation]';
if (alignZ)
    opt_ang{1} = 1;    
    for i = 2:max_fiber
        tmp = angle_options(2:size(opt_angles,1),i-1);
        opt_ang{i} = [ones(size(tmp,1),1) tmp];
    end;
else        
    for i = 1:max_fiber
        opt_ang{i} = angle_options(1:size(opt_angles,1),i);
    end;
end;
% clear out too small crossing angles
if (~isempty(minangle))
    for i = 2:max_fiber    
        for ii = 1:(size(opt_ang{i},2)-1)
            for j = (ii+1):size(opt_ang{i},2)
                sel = sum(angles(:,opt_ang{i}(:,ii)) ...
                    .*angles(:,opt_ang{i}(:,j)),1) ...
                    < cos(minangle/180*pi);
                opt_ang{i} = opt_ang{i}(sel,:);
            end;
        end;
    end;
end;

%% estimate of number of directions
display(['      # water              ' num2str(length(water))]);
display(['      # adcs               ' num2str(length(opt_adcs))]);
display(['      # ratio combinations ' num2str(size(opt_ratios{1},1))]);
display(['      # dispersion options ' num2str(size(opt_dispersions{end},1))]);
display(['      # microstructure conf ' num2str(size(opt_micros{end},1))]);
display(['      # orientations        ' num2str(size(opt_ang{end},1))]);
estim = 0;
nf = sum(opt_ratios{1}>0,2);
for i = 1:max(nf)
    nnf = sum(nf==i);
    estim = estim + nnf*size(opt_dispersions{i},1) ...
        *size(opt_micros{i},1)*size(opt_ang{i},1);
end;
estim = length(water)*length(opt_adcs)*(1+estim);
display(['      large estim # library elements ' num2str(estim)]);

%% prepare list of all voxels

vox = struct;

n = 1;
for w = 1:length(water)
    for a = 1:length(opt_adcs)
        diff = opt_adcs(a);
        diffw = diff;

        fiber_water = [water(w),diffw,diffw,diffw,pi/2,pi/2,0];

        vox(n).ndir = 0;
        vox(n).directions = [];
        vox(n).FA = [];   
        vox(n).ADC = diffw;
        vox(n).fibers = [1,diffw,diffw,diffw,pi/2,pi/2,0];
        vox(n).dirs = [];
        vox(n).disp = [];
        vox(n).microstruct = microstruct;
        switch microstruct
            case 'DTI'
                micro_water = 0;
            case 'twocomp'
                micro_water = [0,diffw,diffw,0];
        end;
        vox(n).micro = micro_water;
        n = n+1;

        if (a == 1), [vox(2:estim).ADC] = deal(0); end;

        for i = 1:size(opt_ratios{w},1)
            fiberratio = opt_ratios{w}(i,:);
            nfiber = sum(fiberratio > 0);
            fiberratio = fiberratio(1:nfiber)';

            microsopt = opt_micros{nfiber};
            angopt = opt_ang{nfiber};
            dispopt = opt_dispersions{nfiber};
            for d = 1:size(dispopt,1)
                disp = dispopt(d,:)';

                for j = 1:size(microsopt,1)
                    micros = microsopt(j,:)';
                    switch microstruct
                        case 'DTI'
                            fa = micros;
                            micro = [micro_water;fa];
                            d2 = sqrt(3)*diff*fa./(sqrt(9-6*fa.^2));
                            l1 = diff + 2*d2;
                            l2 = diff - d2;

                            fibers = [fiber_water;fiberratio,l2,l2,l1,zeros(nfiber,2),disp];
                            ADC = (l1+2*l2)/3;
                        case 'twocomp'
                            % two compartment, see Jelescu2017 eq (14)
                            % (lambda_i, lambda_e, lambda_e_p, v1)
                            micro = reshape(micros,[4,nfiber])';
                            % fa and ADC of external compartment (it's a choice)
                            [fa,ADC] = FAcalc(micro(:,2),micro(:,3),micro(:,3));
                            fibers = [fiber_water;fiberratio,micro(:,3),micro(:,3),micro(:,2),...
                                zeros(nfiber,2),disp];                        
                            micro = [micro_water;micro];
                    end;
                    fiber_angle(1,1:2) = [pi;pi/2];
                    for k = 1:size(angopt,1)
                        dirs = angopt(k,:);
                        fiber_angle = opt_angles(dirs,:);

                        % fibers and water component
                        fibers(1+(1:nfiber),5:6) = fiber_angle;

                        vox(n).ndir = nfiber;
                        vox(n).directions = fiber_angle;
                        vox(n).FA = [0;fa];   
                        vox(n).ADC = [diff;ADC];
                        vox(n).fibers = fibers;
                        vox(n).dirs = dirs;
                        vox(n).disp = [0;disp];
                        vox(n).microstruct = microstruct;
                        vox(n).micro = micro;
                        n = n+1;
                    end;
                end;
            end;
        end;
    end;
end;

optvars = {'anglefile';'max_fiber';'water';'diff';'microlim';'ratiosteps';...
    'microsteps';'adcsteps';'adclim';'alignZ';'maxmicrodiff';'maxratdiff';...
    'maxdispersion';'dispersionsteps';'nmicropar';'maxdispersiondiff';...
    'minangle'};
for i = 1:length(optvars)
    eval(['opt.' optvars{i} ' = ' optvars{i} ';']);
end;

end

function opt_ratios = ratio_options_micro(ratios,nfiber,maxdiff)
    if (nargin < 3) maxdiff = []; end;
    opt_ratios = [];
    if nfiber > 1
        npar = size(ratios,2);
        for ra = 1:size(ratios,1)
            ratios_small = ratio_options_micro(ratios((ra):end,:),nfiber-1,maxdiff);%ratios((ra):end,:);
            opt_ratio_t = [repmat(ratios(ra,:),[size(ratios_small,1),1]) ratios_small];
            opt_ratios = [opt_ratios;opt_ratio_t];
        end;        
        % clean out unwanted combinations
        if (~isempty(ratios))
            sel = ones(size(opt_ratios,1),1) == 1;
            for i = 1:npar
                sel = sel & (max(abs(diff(opt_ratios(:,i:npar:end),[],2)),[],2) < maxdiff(i));
            end;
            opt_ratios = opt_ratios(sel,:);
        end;
    else
        opt_ratios = ratios;
    end;
end

function opt_ratios = ratio_options(ratios,nfiber,clean,sum1,maxdiff)
    if (nargin < 3) clean = true; end;
    if (nargin < 4) sum1 = true; end;
    if (nargin < 5) maxdiff = 0; end;
    opt_ratios = [];
    if nfiber > 1
        for ra = ratios
            ratios_small = ratios((ratios <= ra));
            if (sum1)
                ratios_small = ratios_small((ratios_small <= (ratios(1)*1.001 - ra)));
            end;
            opt_ratio_t = ratio_options(ratios_small,nfiber-1,false,sum1);
            opt_ratio_t = [ones(size(opt_ratio_t,1),1)*ra opt_ratio_t];
            opt_ratios = [opt_ratios;opt_ratio_t];
        end;        
        if (maxdiff > 0)
            opt_ratios = opt_ratios(-min(diff(opt_ratios,[],2),[],2) < maxdiff,:);
        end;
    else
        opt_ratios = ratios';
    end;
    if (clean)
        opt_ratios = opt_ratios(abs(sum(opt_ratios,2) - ratios(1))<0.001,:);
    end;
end

function micro_ratios = micro_options(micros)
    if (length(micros) > 1)
        micro_ratios = [];
        tmp = micro_options(micros(2:end));
        for l = 1:length(micros{1})
            micro_ratios = [micro_ratios;...
                micros{1}(l)*ones(size(tmp,1),1),tmp];
        end;
    else
        micro_ratios = micros{1}';
    end;
end

function opt_ang = angle_options(angles,nfiber)
    opt_ang = [];
    if nfiber > 1
        for ang = angles
            angles_small = angles(angles~=ang);
            opt_ang_t = angle_options(angles_small,nfiber-1);
            opt_ang_t = [ones(size(opt_ang_t,1),1)*ang opt_ang_t];
            opt_ang = [opt_ang;opt_ang_t];
        end;
    else
        opt_ang = angles';
    end;
end
