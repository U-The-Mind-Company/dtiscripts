%% Run fiber identification on a fib.gz-file

% Steven Baete
% November 2016
% NYU SOM CBI

%   Input:
%       fibfile: location of fib.gz-file
%       method: fiber identification method
%           'FP', 'Regular', 'MrTrix', 'FSL', 'FSLDTI','FSLQBOOT','msmtCSD'
%       opt: option structure
%           For 'FP':
%               opt.libname: filename of ODF FP library
%               opt.matchmethod: ODF Fingerprinting option
%                   0: without fiber complexity penalty
%                   7: with fiber complexity penalty (default), 
%                       needs srcfile-input
%               opt.srcfile: filename of src.gz-file (for noise estimation)
%   Output:
%       fibfileout: name of the created fib.gz-file

function fibfileout = run_on_fib(fibfile,method,opt)
fstr = 'run_on_fib';

MAX_FIBERS = 5;
matchmethod = -1;

%% check input
if (~(exist(fibfile,'file')==2))
    error([fstr ' : not a valid fibfile: ' fibfile '.']);
end;
if (nargin < 2)
    error([fstr ' : did not receive which method to run.']);
end;
if (nargin < 3)
    opt = [];
end;
switch(method)
    case 'FP'
        if ((nargin < 3) || isempty(opt))
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with fields libname.']);
        end;
        if ((~isfield(opt,'libname')) || (isempty(opt.libname)) ...
                || (~(exist(opt.libname,'file')==2)) )
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with a field libname pointing to a valid library.']);
        end;
        if ((~isfield(opt,'matchmethod')) || (isempty(opt.libname)) )
            matchmethod = 7;
        else
            matchmethod = opt.matchmethod;
        end;
        if (matchmethod == 7 & ( (~isfield(opt,'srcfile')) || (isempty(opt.srcfile)) ...
                || (~(exist(opt.srcfile,'file')==2)) ) )
            if (~(exist(opt.srcfile,'file')==2))
                display([' Could not find srcfile : ' opt.srcfile]);
            end;
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with a field srcfile pointing to a valid srcfile.']);
        end;
        methodtag = ['fp' num2str(matchmethod)];        
    case 'Regular'
        methodtag = 'reg';    
    case 'MrTrix'
        methodtag = 'trx';  
    case 'FSL'
        if ((nargin < 3) || isempty(opt))
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with fields srcfile.']);
        end;
        if ((~isfield(opt,'srcfile')) || (isempty(opt.srcfile)) ...
                || (~(exist(opt.srcfile,'file')==2)) )
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with a field srcfile pointing to a valid srcfile.']);
        end;
        methodtag = 'fsl'; 
    case 'FSLDTI'
        if ((nargin < 3) || isempty(opt))
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with fields srcfile.']);
        end;
        if ((~isfield(opt,'srcfile')) || (isempty(opt.srcfile)) ...
                || (~(exist(opt.srcfile,'file')==2)) )
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with a field srcfile pointing to a valid srcfile.']);
        end;
        methodtag = 'fsldti';
    case 'FSLQBOOT'
        if ((nargin < 3) || isempty(opt))
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with fields srcfile.']);
        end;
        if ((~isfield(opt,'srcfile')) || (isempty(opt.srcfile)) ...
                || (~(exist(opt.srcfile,'file')==2)) )
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with a field srcfile pointing to a valid srcfile.']);
        end;
        methodtag = 'fslqboot';
    case 'msmtCSD'
        if ((nargin < 3) || isempty(opt))
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with fields srcfile.']);
        end;
        if ((~isfield(opt,'srcfile')) || (isempty(opt.srcfile)) ...
                || (~(exist(opt.srcfile,'file')==2)) )
            error([fstr ' : method ''' method ''' expects as third input '...
                'a structure with a field srcfile pointing to a valid srcfile.']);
        end;
        methodtag = 'csd';
    otherwise
        error([fstr ' : ''' method ''' is not a valid method. Options are FP, Regular']);
end;
%% load fib-file

fib = fib_obj(fibfile);

%% copy fib-file

fibfileout = strrep(fibfile,'.fib.gz',['.' methodtag '.fibt']);
unix(['rm ' fibfileout]);
fibcp = fib.copy(fibfileout);
clear fib;

%% load the FP library

if (strcmp(method,'FP'))
    lib = load(opt.libname);
    odflib = lib.odfrot;
    dirslib = lib.dirrot;
    adclib = lib.adc;
    if (~isfield(lib,'micro'))
        microlib = lib.fa;
    else
        microlib = lib.micro;
    end;
    if (isfield(lib,'disp'))
        displib = lib.disp;
    else
        displib = zeros(size(adclib));
    end;
    ratlib = lib.rat;
    clear lib;
end;

%% mask

fa0 = fibcp.getvar('fa0');
mask = (fa0~=0);
% mask(:,:,1:23) = 0;
% mask(:,:,25:end) = 0;

if (isfield(opt,'slices') & ~isempty(opt.slices))
    if (~isempty(opt.slices{1}))
        sel = (ones(1,size(mask,1))==true);
        sel(opt.slices{1}) = false;
        mask(sel,:,:) = 0;
    end;
    if (~isempty(opt.slices{2}))
        sel = (ones(1,size(mask,2))==true);
        sel(opt.slices{2}) = false;
        mask(:,sel,:) = 0;
    end;
    if (~isempty(opt.slices{3}))
        sel = (ones(1,size(mask,3))==true);
        sel(opt.slices{3}) = false;
        mask(:,:,sel) = 0;
    end;
end;

%% load the src-file and bvals

if (strcmp(method,'FSL') | strcmp(method,'FSLDTI') | ...
        strcmp(method,'FSLQBOOT') | strcmp(method,'msmtCSD') | ...
        (strcmp(method,'FP') & matchmethod == 7))
    src = src_obj(opt.srcfile);
    b_table = src.b_table;
    dwi = src.getvar('image');
    sdwi = size(dwi);
    if (strcmp(method,'msmtCSD'))
        dwiall = dwi;
        frfnames{1} = [opt.srcfile(1:(end-7)) '.wm.txt'];
        frfnames{2} = [opt.srcfile(1:(end-7)) '.gm.txt'];
        frfnames{3} = [opt.srcfile(1:(end-7)) '.csf.txt'];
    end;
    dwi = reshape(dwi,[prod(sdwi(1:3)),sdwi(4)]);
    dwi = dwi(mask(:),:);
    if ( strcmp(method,'FSL') | strcmp(method,'FSLDTI') | ...
        strcmp(method,'FSLQBOOT')  | strcmp(method,'msmtCSD'))
        if (regexp(fibfile,'\.f\d'))
            ind = regexp(fibfile,'\.f\d')+1;
            dotind = strfind(fibfile,'.');
            flipstr = fibfile((dotind(find(dotind < ind,1,'last'))+2):(dotind(find(dotind > ind,1,'first'))-1));
            %flipstr = flipstr(flipstr < 3);
            for ii = 1:length(flipstr)
                if (str2num(flipstr(ii)) < 3)
                    b_table(1+str2num(flipstr(ii)),:) = -b_table(1+str2num(flipstr(ii)),:);
                end;
            end;
        end;
    end;
end;

%% read odfs

odf_vertices = fibcp.odf_vertices;
odf_faces = fibcp.odf_faces;
odf = fibcp.getvar('odf');
s = size(odf,1:4);
odf = reshape(odf,[s(1),prod(s(2:4))]);
odfm = odf(:,mask(:));
sm = size(odfm);
clear odf;

%% initialize matrices
noisestdm = zeros(sm(2),1);
for i = 1:(MAX_FIBERS+1) % index 1 is the free water
    fam{i} = zeros(sm(2),1);
    indexm{i} = zeros(sm(2),1);
    adcm{i} = zeros(sm(2),1);
    dispm{i} = zeros(sm(2),1);
    ratm{i} = zeros(sm(2),1);
    if (strcmp(method,'FP'))
        microm{i} = zeros(sm(2),size(microlib,3));
    end;
end;
if (strcmp(method,'FP'))
    rmsem = zeros(sm(2),1);
end;

%% loop over the odfs in blocks

blocksize = 1e4;
tic;
for il = 1:blocksize:sm(2)
    %display(['    identifying fiber directions : [' num2str(il) '/' num2str(sm(2)) ']']);
    blockind = (il -1) + (1:blocksize);
    blockind = blockind(blockind <= sm(2));
    bsize = length(blockind);
    
    %% match the odfs with the best fit and find their FP-directions
    odfmt = odfm(:,blockind)';
    odfmt(odfmt < 0 ) = 0;

    switch (method)
        case 'FP'
            if (matchmethod == 6 | matchmethod == 7)
                dwimt = dwi(blockind,:);
                var = dwi_variance(dwimt,b_table);
                [fpdirst,fpmicro,fpadc,fprat,fpdisp,fprmse] = ...
                 find_ODF_peak_FP(odfmt,odf_vertices,[],odflib,...
                    dirslib,microlib,MAX_FIBERS,adclib,ratlib,...
                    matchmethod,var,displib);
                noisestdm(blockind) = var;
            else         
                [fpdirst,fpmicro,fpadc,fprat,fpdisp,fprmse] = ...
                 find_ODF_peak_FP(odfmt,odf_vertices,[],odflib,...
                    dirslib,microlib,MAX_FIBERS,adclib,ratlib,...
                    matchmethod,[],displib);
            end;
        case 'Regular'
            [fpdirst] = find_ODF_peak(odfmt,odf_faces,odf_vertices);
        case 'MrTrix'
            [fpdirst] = find_ODF_peak_mrtrix(odfmt,odf_vertices);
        case 'FSL'
            dwimt = dwi(blockind,:);
            [fpdirst] = find_ODF_peak_bedpostx(dwimt,b_table);
        case 'FSLDTI'
            dwimt = dwi(blockind,:);
            [fpdirst] = find_ODF_peak_dtifit(dwimt,b_table);
        case 'FSLQBOOT'
            dwimt = dwi(blockind,:);
            [fpdirst] = find_ODF_peak_qboot(dwimt,b_table);
        case 'msmtCSD'
            dwimt = dwi(blockind,:);
            [fpdirst] = find_ODF_peak_msmtCSD(dwimt,b_table,dwiall,frfnames);
    end;
    
    %% calculate the QA and index
    sfpd = size(fpdirst);
    fpdirst = reshape(fpdirst,[sfpd(1)*sfpd(2),sfpd(3)]);
    [x,y,z] = sph2cart(fpdirst(:,1),fpdirst(:,2),ones(sfpd(1)*sfpd(2),1));
    [val,fpindt] = min((repmat(odf_vertices(1,:),[sfpd(1)*sfpd(2),1])-repmat(x,[1,size(odf_vertices,2)])).^2 ...
        + (repmat(odf_vertices(2,:),[sfpd(1)*sfpd(2),1])-repmat(y,[1,size(odf_vertices,2)])).^2 ...
        +(repmat(odf_vertices(3,:),[sfpd(1)*sfpd(2),1])-repmat(z,[1,size(odf_vertices,2)])).^2,[],2);
    fpindt(isnan(val)) = -1;
    fpindt(fpindt > size(odf_vertices,2)/2) = fpindt(fpindt > size(odf_vertices,2)/2) - size(odf_vertices,2)/2;
    fpindt = reshape(fpindt,[sfpd(1),sfpd(2)]);
        
    sodfmt = size(odfmt);
    for i = 1:MAX_FIBERS
        indexm{i}(blockind) = max(fpindt(:,i),1) - 1;
        fatmp = odfmt(sub2ind(sodfmt,1:sodfmt(1),max(fpindt(:,i)',1))) - min(odfmt,[],2)';
        fam{i}(blockind) = fatmp.*(fpindt(:,i)>0)';
    end;
    if (strcmp(method,'FP'))
        for i = 1:(MAX_FIBERS+1) % 1 is the free water
            if (exist('fpadc','var')), adcm{i}(blockind) = fpadc(:,i);  end;
            if (exist('fprat','var')), ratm{i}(blockind) = fprat(:,i);  end;
            if (exist('fpdisp','var')), dispm{i}(blockind) = fpdisp(:,i);  end;
            if (exist('fpmicro','var')), microm{i}(blockind,:) = squeeze(fpmicro(:,i,:));  end;
        end;
        if (exist('fprmse','var')), rmsem(blockind) = fprmse(:);  end;
    end;
    display(['  run ' num2str(blockind(end)) '/' num2str(sm(2)) '  ' datestr(toc/1000/60,'dd HH:MM:SS')]);
end;
clear odflib; clear odfm; clear dwimt;

%% write to fib-file
setstr = [];clearstr = [];
noisestd = zeros(size(mask));
noisestd(mask == 1) = noisestdm;
setstr = strcat(setstr,',''noisestd'',noisestd');
clearstr = strcat(clearstr,' noisestd');

for ff = 0:(MAX_FIBERS-1)
    eval(sprintf('fa%d = zeros(size(mask));', ff));
    eval(sprintf('fam{%d}(fam{%d} <= 0) = 1e-4;', ff+1, ff+1));
    eval(sprintf('fa%d(mask == 1) = fam{%d};', ff, ff+1));
    eval(['setstr = strcat(setstr,'',''''fa' num2str(ff) ''''',fa' num2str(ff) ''');']);
    eval(['clearstr = strcat(clearstr,'' fa' num2str(ff) ''');']);
%     eval(sprintf('fibcp.setvolume(''fa%d'',fa%d);', ff,ff));
%     eval(sprintf('clear fa%d', ff));
    eval(sprintf('index%d = zeros(size(mask));', ff));
    eval(sprintf('index%d(mask == 1) = indexm{%d};', ff, ff+1));
    eval(['setstr = strcat(setstr,'',''''index' num2str(ff) ''''',index' num2str(ff) ''');']);
    eval(['clearstr = strcat(clearstr,'' index' num2str(ff) ''');']);
%     eval(sprintf('fibcp.setvolume(''index%d'',index%d);', ff,ff));
%     eval(sprintf('clear index%d', ff));
end;
% save microstructure parameters
if (strcmp(method,'FP'))
    pars ={'adc';'rat';'disp';'micro'};
    for ff = 0:(MAX_FIBERS)
        if (ff == 0), comp = 'fw_';
        else, comp = ['fib' num2str(ff-1)];  end;
        for j = 1:length(pars)
            eval(['np = size(' pars{j} 'm{' num2str(ff+1) '},2);']);
            for k = 1:np
                if (np > 1), kstr = ['_' num2str(k)]; else, kstr = ''; end;
                eval([comp pars{j} kstr ' = zeros(size(mask));']);
                eval(sprintf([comp pars{j} kstr '(mask == 1) = ' pars{j} 'm{%d}(:,k);'], ff+1));
                eval(['setstr = strcat(setstr,'',''''' comp pars{j} kstr ''''',' ...
                    comp pars{j} kstr ''');']);
                eval(['clearstr = strcat(clearstr,'' ' comp pars{j} kstr ''');']);
            end;
        end;
    end;
	rmse = zeros(size(mask));
    rmse(mask == 1) = rmsem;
    eval(['setstr = strcat(setstr,'',''''rmse'''',rmse'');']);
    eval(['clearstr = strcat(clearstr,'' rmse'');']);

end;
i = 0;
while ~isempty(whos(fibcp.m,['odf' num2str(i)]))
    setstr = [setstr ',''odf' num2str(i) ''',[]'];
    i = i+1;
end;
eval(['fibcp.setvolume(' setstr(2:end) ');']);
%figure;imagesc(montageSB(sqrt(noisestd)))
eval(['clear' clearstr ]);

%% close the new fib-file
fibcp.delete;
clear fibcp;

fibfileout = strrep(fibfileout,'fibt','fib.gz');
