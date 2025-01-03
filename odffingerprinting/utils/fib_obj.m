% Steven Baete
% NYU LMC CBI
% Februari 2015

classdef fib_obj < handle
    % class to hold information about and read from fib.gz-files
    % (DSIStudio)

    properties
        output;
    end

    properties(GetAccess='public', SetAccess='protected')
        m;
        varnamesall;
        varnamesfull;
        varnamesequal;
        varnamesmasked;
        varnamesmaskedsize;

        dimension;
        voxel_size;
        odf_faces;
        odf_vertices;
        file;
        fileunzip;
        filezipped;
        changed;
    end

    methods
        % Constructor
        function this = fib_obj(file)
            % the user used a wildcard in the folder-name
            if (length(strfind(file,'*')) > 0)
                file = filestar(file);
            end;
            file = strrep(file,[filesep filesep],[filesep]);
            
            if ~(exist(file)==2)
                display(['    fib_obj - did not find ' file]);
                error(['fib_obj: did not find requested file ' file]);
            end;

            this.output = true;
            if (this.output);display(['    fib_obj - reading ' file]);end;
            % unzip
            if (~isempty(strfind(file((end-5):end),'.gz')))
                this.filezipped = true;
                if (this.output);display(['     unzipping gz-file']);end;
                if (isunix == 1)
                    fileunzip = [file(1:(end-3)) 't'];
                    if (~(exist(fileunzip)==2))
                        unix(['gunzip -c ' file ' > ' fileunzip]);
                    else % check if the gz file is not newer that the fibt-file
                        FileInfo = dir([file]);
                        TimeStamp = FileInfo.datenum;
                        FileInfo = dir([fileunzip]);
                        TimeStampUnzip = FileInfo.datenum;
                        if (TimeStamp > TimeStampUnzip)
                            unix(['gunzip -c ' file ' > ' fileunzip]);
                        end;
                    end;
                else
                    fileunzip = cell2mat(gunzip(file));
                end;
            else
                this.filezipped = false;
                fileunzip = file;
            end;

            this.file = file;
            this.fileunzip = fileunzip;

            % read contents of file and index
            this = this.readandindex;

            this.disp;
            this.changed = false;
        end;
            
        % Read contents of file and index
        function this = readandindex(this)
            % what is in our file
            this.m = matfile(this.fileunzip,'Writable',true);

            % read in some parameters
            this.dimension = this.m.dimension;
            this.voxel_size = this.m.voxel_size;
            this.odf_faces = this.m.odf_faces;
            this.odf_vertices = this.m.odf_vertices;

            % make a list of the varnames
            this.varnamesall = whos(this.m);
            % make a list of the available varnames for reading
            leq = 1;lfull = 1;this.varnamesmasked = {};this.varnamesmaskedsize = {};
            verts = size(this.odf_vertices,2);
            for i = 1:length(this.varnamesall)
                varname = this.varnamesall(i).name;
                varsize = this.varnamesall(i).size;
                % fields which need to be copied
                if ((varsize(1) < 4 && varsize(2) < 3*verts) || (varsize(2) == 3) || (varsize(2) == 1)) % such as: dimension, idxb0,idcb0compl,idxb0echo,max_odf,odf_faces,odf_vertices,voxel_size
                    this.varnamesequal{leq} = varname;leq = leq+1;
                    continue;
                end;
                % fields which are saved fully
                if ((varsize(1) == 1 && varsize(2) == prod(this.dimension)) || (varsize(1) > 2*verts))
                    this.varnamesfull{lfull} = varname;lfull = lfull+1;
                    continue;
                end;
                % fields which are saved masked
                if (varsize(1) > 1)% && varsize(2) > 2000)
                    % do we already have a list of these?
                    found = false;
                    for j = 1:length(this.varnamesmasked)
                        ll = min(3,length(varname));
                        ll2 = min(3,length(this.varnamesmasked{j}{1}));
                        if (~found && strcmp(this.varnamesmasked{j}{1}(1:ll2),varname(1:ll)))
                            this.varnamesmasked{j}{length(this.varnamesmasked{j})+1} = varname;
                            this.varnamesmaskedsize{j}(size(this.varnamesmaskedsize{j},1)+1,:) = varsize;
                            found = true;
                        end;
                    end;
                    if (~found)
                        this.varnamesmasked{length(this.varnamesmasked)+1}{1} = varname;
                        this.varnamesmaskedsize{length(this.varnamesmaskedsize)+1}(1,:) = varsize;
                    end;
                    continue;
                end;
                if (this.output);display(['    fib_obj - I could not classify variable ' varname ' (size ' num2str(varsize) ') in file ' fileunzip '. Skipping from read!']);end;
            end;
            % sort the varnamesmasked
            for i = 1:length(this.varnamesmasked)
                temp = this.varnamesmasked{i};
                tempsize = this.varnamesmaskedsize{i};
                for j = 1:length(temp)
                    num = str2num(temp{j}(4:end));
                    if (isempty(num))
                        num = 0;
                    end;
                    this.varnamesmasked{i}{num+1} = temp{j};                    
                    this.varnamesmaskedsize{i}(num+1,:) = tempsize(j,:);
                end;
            end;
        end;

        % Destructor
        function delete(this)
            if (this.filezipped || this.changed)
                if (this.output);display(['    fib_obj - deleting original unzipped file ' this.file]);end;
                if (this.changed) % rezip the fib-file
                    if (isunix == 1)
                        o = unix(['pigz -c ' this.fileunzip ' > ' this.file]);
                        if (o==127) % no pigz on this system
                            unix(['gzip -c ' this.fileunzip ' > ' this.file]);
                        end;
                    else
                        gzip(file);
                    end;
                end;
                delete(this.fileunzip);
                if (this.output);display(['              ... done']);end;
            end;
        end;

        function disp(this)
            display(['  Fib_obj from ' this.file]);
            display(['     dimension [' num2str(this.dimension,'% 4.0f') '], voxel_size [' num2str(this.voxel_size,'% 5.2f') ']']);
            fprintf(['     available ']);fprintf('%s, ',this.varnamesequal{1:(end-1)});fprintf([this.varnamesequal{end} '\n']);
            fprintf(['               ']);fprintf('%s, ',this.varnamesfull{1:(end-1)});fprintf([this.varnamesfull{end} '\n']);
            if (length(this.varnamesmasked) > 1)
                fprintf(['               ']);fprintf('%s, ',this.varnamesmasked{1:(end-1)}{1}(1:(end-1)));fprintf([this.varnamesmasked{end}{1}(1:(end-1)) '\n']);
            else
                if (length(this.varnamesmasked) > 0)
                    fprintf(['               ']);fprintf([this.varnamesmasked{end}{1}(1:(end-1)) '\n']);
                end;
            end;
        end;

        % read the fields from the fib-file and return the reshaped vars
        function varargout = subsref(this,S)
            switch S(1).type
                case '.'
                    % is this one of the fib.gz variables?
                    varloc = this.findvar(S(1).subs);
                    if (~isempty(varloc))
                        varargout{:} = this.getvar(S(1).subs);
                        return;
                    end
                    if (strcmp(S(1).subs,'montage'))
                        for i = (length(S(2).subs)+1):4
                            S(2).subs{i} = [];
                        end
                        builtin('subsref',this,S);
                    end;
                    if nargout == 0
                        builtin('subsref', this, S);
                    else
                        varargout      = cell(1, nargout);
                        [varargout{:}] = builtin('subsref', this, S);
                    end
                    return;
                case {'()','{}'}
                otherwise
                    error('operator not supported');
            end;
        end;
        
        function isf = isvar(this,var)
            isf = ~isempty(findvar(this,var));
        end;
        
        function [varloc,i] = findvar(this,var)
            % variable is on the varnamesequal list
            for i = 1:length(this.varnamesequal)
                if (strcmp(var,this.varnamesequal{i}))
                    varloc = 'equal';
                    return;
                end;
            end;
            % variable is on the varnamesfull list
            for i = 1:length(this.varnamesfull)
                if (strcmp(var,this.varnamesfull{i}))
                    varloc = 'full';
                    return;
                end;
            end;
            % variable is on the varnamesmasked list
            for i = 1:length(this.varnamesmasked)
                if (strcmp(var,this.varnamesmasked{i}{1}(1:(end-1))))
                    varloc = 'masked';
                    return;
                end;
            end;
            varloc = '';
            return;
        end;

        function outvar = getvar(this,var)
            [varloc,fullind] = this.findvar(var);

            if (strcmp(varloc,'equal') || strcmp(varloc,'full'))
                eval(['outvar = this.m.' var ';']);
                if (strcmp(varloc,'full'))
                    outvar = reshape(outvar,this.dimension);
                end;
                return
            end;
            if (strcmp(varloc,'masked'))
                str = load(this.fileunzip,'-mat');
                tempvar = zeros(this.varnamesmaskedsize{fullind}(1,1),...
                    sum(this.varnamesmaskedsize{fullind}(:,2)));
%                 tempvar = [];
                varsize(1) = 0;
                for j = 1:length(this.varnamesmasked{fullind})
                    varname = this.varnamesmasked{fullind}{j};
%                     eval(['tempvar = [tempvar this.m.' varname '];']);
%                     eval(['temptemp = this.m.' varname ';']);
                    eval(['temptemp = str.' varname ';']);
                    tempvar(:,sum(this.varnamesmaskedsize{fullind}(1:(j-1),2)) + ...
                        (1:this.varnamesmaskedsize{fullind}(j,2))) = temptemp;
                    varsize(j+1) = sum(this.varnamesmaskedsize{fullind}(1:j,2)) ...
                                -sum(varsize(1:j));
                    if (varsize(j+1) < varsize(2))
                        break;
                    end;
                end;
                clear str;
                nvol = size(tempvar,1);
                outvar = zeros([nvol,prod(this.dimension)]);
                %outvar = zeros([nvol,this.dimension]);
                try
                    fa0 = this.m.fa0;
                    non_zero_index = (fa0(:) ~= 0.00);
                    clear fa0;
                    if (sum(non_zero_index) == (length(tempvar) - varsize(2)))
                        tempvar = tempvar(:,(varsize(2)+1):end);
                    end;
                    if (this.changed && sum(non_zero_index) ~= length(tempvar))
                        fa0 = this.m.fa0orig;
                        non_zero_index = (fa0(:) ~= 0.00);
                        clear fa0;
                    end;
                    %outvar(:,non_zero_index) = tempvar;
                    tt = find(non_zero_index);
                    tempvar(:,find(sum(tempvar,1)==0)) = [];
                    outvar(:,tt) = tempvar;
                catch
                    try
                        try
                            fa0 = this.m.FAq;
                            non_zero_index = (fa0(:) ~= 0.00);
                            clear fa0;
                            outvar(:,non_zero_index) = tempvar;
                        catch
                            try
                                fa0 = this.m.indexdti;
                                non_zero_index = (fa0(:) ~= 0.00);
                                clear fa0;
                            catch
                                fa0 = this.m.fa0orig;
                                non_zero_index = (fa0(:) ~= 0.00);
                                clear fa0;
                                tempvar = tempvar(:,sum(tempvar,1)>0);
                            end;
                            outvar(:,non_zero_index) = tempvar;
                        end;
                    catch
                        fa0 = this.m.fa0;fa1 = this.m.fa1;fa2 = this.m.fa2;
                        non_zero_index = (fa0(:) ~= 0.00 | fa1(:) ~= 0.00 | ...
                            fa2(:) ~= 0.00);
                        clear fa0 fa1 fa2;
                        tempvar = tempvar(:,sum(tempvar,1)>0);
                        outvar(:,non_zero_index) = tempvar;
                    end;
                end;
                clear fa0;
                outvar = reshape(outvar,[nvol,this.dimension]);
            end;
        end;

%             See how MATLAB calls subsref for the expression:
%             A.field
%             The syntax A.field calls B = subsref(A,S) where S.type='.' and S.subs='field'.
        function imagine(this,volume)
            try
                voldat = this.getvar(volume);
                imagine(permute(voldat,[2,1,3]));
            catch
                display([' Error: Please choose an available volume!']);
                this.disp();
            end;
            return
        end;

        function montage(this,volume,slices,msize,sel,varargin)
            cmax = 0;
            for i=1:2:size(varargin,2)
                if (~ischar(varargin{i}))
                    error('fib_obj, montage: Please use string-value pairs for input');
                end;
                switch varargin{i}
                    case 'cmax'
                        cmax = varargin{i+1};
                    otherwise
                        error(['EPIreconMB: I did not recognize the input-string ' varargin{i}]);
                end;
            end;

            try
                voldat = this.getvar(volume);
            catch
                display([' Error: Please choose an available volume!']);
                this.disp();
                return
            end;
            X = zeros([this.dimension([2,1]),1,this.dimension(3)]);
            X(:,:,1,:) = permute(voldat,[2,1,3]);
            X = flipdim(X,4);
            if (isempty(slices))
                slices = 1:size(X,4);
            end;
            if (length(slices) ==3)
                s = size(slices);
                if (isempty(slices{3}))
                    if (isempty(slices{1}))
                        X = permute(X,[1,4,3,2]);
                        slicest = slices{2};
                    else
                        X = permute(X,[4,2,3,1]);
                        slicest = slices{1};
                    end;
                else
                    slicest = slices{3};
                end;
                clear slices;slices = slicest;
            end;
            slices(slices < 1) = 1;slices(slices > size(X,4)) = size(X,4);
            slices = size(X,4) - slices+1;
            if (isempty(msize))
                msize = ceil(sqrt(length(slices)))*[1,1];
            else
                msize(msize < 1) = 1;
                msize(2) = ceil(length(slices)/msize(1));
            end;
            if (isempty(sel))
                sel(1:2) = [1,size(X,1)];
                sel(3:4) = [1,size(X,2)];
            else
                sel = max(sel,1);
                sel([1,2]) = min(sel([1,2]),size(X,1));
                sel([3,4]) = min(sel([3,4]),size(X,2));
            end;
            X = X(sel(1):sel(2),sel(3):sel(4),:,:);
            if (cmax == 0)
                cmax = max(X(:));
            end;
            montage(X,'Indices',slices,'Size',msize);caxis([0 cmax]);
            ylim([1 msize(1)*size(X,1)]);xlim([1 msize(2)*size(X,2)]);
            return
        end;
        function plot_odf(this,pos)
            try
                odf = this.getvar('odf');
            catch
                display([' Error: Could not load odfs!']);
                this.disp();
                return
            end;
            s = size(odf);
            pos = max(pos,1);pos = min(pos,s(2:4));
            load('odf8.mat');
            odf = odf(:,pos(1),pos(2),pos(3));
            odf_points = odf_vertices.*repmat(odf-min(odf),[2 3])';
            polygon.vertices = odf_points';
            polygon.faces = odf_faces'+1;
            polygon.facevertexcdata = abs(odf_vertices)';
            %figure;
            p = patch(polygon);
            set(p,'FaceColor','flat','EdgeColor','none');
            daspect([1 1 1])
            view(3); axis tight
            camlight
            lighting gouraud
            hold on;
            lax = max(odf_points,[],2)*1.20;
            lax(lax < max(lax)) = (lax(lax < max(lax)) + max(lax))/2;
            plot3([0,-lax(1)],[0,0],[0,0],'-k','Color','w');text(-lax(1)*1.05,0,0,'L','Color','w');
            plot3([0,0],[0,lax(2)],[0,0],'-k','Color','w');text(0,lax(2)*1.05,0,'A','Color','w');
            plot3([0,0],[0,0],[0,lax(3)],'-k','Color','w');text(0,0,lax(3)*1.05,'S','Color','w');
            return;
        end;
        function plotmicro(this,slice,famaskvalue)
            if (nargin < 2 | isempty(slice)), slice = floor(this.dimension(3)/2); end;
            if (nargin < 3 | isempty(famaskvalue)), famaskvalue = 0; end;
            % max number of fibers
            for i = 0:4
                tmp = this.getvar(['fa' num2str(i)]);
                mask{i+1} = tmp(:,slice,:) > famaskvalue;
                if (sum(tmp(:)>1e-4))
                    nfib = i+1;
                end;
            end;
            mask{6} = mask{1};
            rmseplot = 0;
            % number of micro-params
            nmicro = 0;
            for i = 1:length(this.varnamesfull)
                if (~isempty(strfind(this.varnamesfull{i},'fib0'))) 
                    nmicro = nmicro+1; 
                    microvar{nmicro} = this.varnamesfull{i};
                end;
                if (~isempty(strfind(this.varnamesfull{i},'rmse')))
                    rmseplot = 1;
                end;
            end;
            figure;
            l = 0;
            for i = 1:(nfib+1)
                if (i <= nfib), fibstr = ['fib' num2str(i-1)]; else, fibstr = 'fw_'; end;
                for j = 1:(nmicro+1)
                    l = l+1;
                    subplot(nfib+1,nmicro+1,l);
                    if (j > 1), parstr = strrep(microvar{j-1},'fib0',fibstr); 
                        else parstr = ['fa' num2str(i-1)]; end;
                    if ((i == nfib+1) && (j == 1))
                        parstr = 'rmse';
                    end;         
                    if (i <= nfib | j > 1 | rmseplot)
                        tmp = this.getvar(parstr);
                        imagesc(flipdim(squeeze(tmp(:,slice,:).*mask{i})',1));
                        title(strrep(parstr,'_','\_'));colormap bone;
                        axis equal;colorbar;
                    end;                    
                    axis off;
                end;
            end;        
        end;            
        function obj = copy(this,newfilename)
            obj = this;

            % copy unzipped fibfile
            if (obj.output);display(['     copying fib-file']);end;
            if (isunix == 1)
                if (~(exist(newfilename)==2))
                    unix(['cp ' obj.fileunzip ' ' newfilename]);
                else
                    display(['     fib-file already exists ...']);
                end;
            end;
            if (this.filezipped)
                delete(this.fileunzip);
                if (this.output);display(['              ... done']);end;
            end;
            % change the fibfilename
            obj.fileunzip = newfilename;
            obj.file = [obj.fileunzip(1:(end-1)) '.gz'];
            obj.changed = true;
            obj.m = matfile(obj.fileunzip,'Writable',true);
        end;

        function setvolume(this,varargin)
            this.changed = true;
            m2 = load(this.fileunzip,'-mat');
            for vi=1:2:size(varargin,2)
                var = varargin{vi};
                vardata = varargin{vi+1};
                if (this.output);display(['     setting volume ' var]);end;
                
                [varloc,fullind] = this.findvar(var);

                % remove what was there before
                if (~isempty(varloc))
                    if (strcmp(varloc,'equal') || strcmp(varloc,'full'))
                        % rmfield(this.m,var);
                        % eval(['this.m.' var ' = [];']);
                        if strcmp(varloc,'equal')
                            this.varnamesequal{fullind} = [];
                        end;
                        if strcmp(varloc,'full')
                            this.varnamesfull{fullind} = [];
                        end;
                    end;
                end;

                if (isempty(varloc))
                    s = size(vardata);
                    if (s(1) == this.dimension(1) || length(vardata) == prod(this.dimension))
                        varloc = 'full';
                        this.varnamesfull{length(this.varnamesfull)+1} = var;
                        fullind = length(this.varnamesfull);
                    else
                        varloc = 'masked';
                        this.varnamesmasked{length(this.varnamesmasked)+1}{1} = var;
                        fullind = length(this.varnamesmasked);
                    end;
                end;

                if (strcmp(varloc,'equal') || strcmp(varloc,'full'))
                    if strcmp(varloc,'equal')
                        this.varnamesequal{fullind} = var;
                    end;
                    if strcmp(varloc,'full')
                        this.varnamesfull{fullind} = var;
                    end;
                    vardata = reshape(vardata,[1,prod(size(vardata))]);
                    if (strcmp(var,'fa0') & ~isfield(m2,'fa0orig'))
                        m2.fa0orig = m2.fa0;
                        this.varnamesfull{length(this.varnamesfull)+1} = 'fa0orig';
                    end
                    eval(['m2.' var ' = vardata;']);
    %                 save(this.fileunzip,'-struct','m2','-mat', '-v4');
    %                 this = this.readandindex;
    %                 this = clearemptyfields(this);
    %                 return
                end;
                if (strcmp(varloc,'masked'))
                    % reshape our variable
                    nvol = size(vardata,1);
                    vardata = reshape(vardata,[nvol,prod(this.dimension)]);
                    % mask our variable
                    try
                        fa0 = m2.fa0;
                        non_zero_index = (fa0(:) ~= 0.00);
                        vardata = vardata(:,non_zero_index);
                    catch
                        try
                            fa0 = m2.indexdti;
                            non_zero_index = (fa0(:) ~= 0.00);
                        catch
                            fa0 = m2.fa0;
                            non_zero_index = (fa0(:) ~= 0.00);
                        end;
                        vardata = vardata(:,non_zero_index);
                    end;
                    % size of the saved variables
                    varname = this.varnamesmasked{fullind}{1};
                    try
                        eval(['tempvar = size(m2.' varname ');']);
                    catch
                        varnametmp = this.varnamesmasked{fullind}{1};
                        eval(['tempvar = size(m2.' varnametmp ');']);
                    end;
                    blocksize = tempvar(2);
                    % remove the old variable
                    for i = 2:length(this.varnamesmasked{fullind})
                        % rmfield(m2,this.varnamesmasked{fullind}{i});
                        eval(['m2.' this.varnamesmasked{fullind}{i} '=[];']);
                    end;
                    % save our variable
                    this.varnamesmasked{fullind}(2:end) = [];
                    l = 0;
                    for j = 1:blocksize:size(vardata,2)
                        varname = strrep(this.varnamesmasked{fullind}{1},'0',num2str(l));
                        sel = j-1+(1:blocksize);
                        sel = sel(sel <= size(vardata,2));
                        eval(['m2.' varname ' = vardata(:,sel);']);
                        this.varnamesmasked{fullind}{l+1} = ...
                            varname;
                        l = l+1;
                    end;
                end;
            end;
            
            fieldsrm = {};
            % which fields to delete
            fie = fields(m2);
            for i = 1:length(fie)
              eval(['val=isempty(m2.' fie{i} ');']);
              if (val)
                fieldsrm{end+1} = fie{i};
              end;
            end;

            if ~isempty(fieldsrm)
              % the actual deletion
              %str = load(this.fileunzip,'-mat');
              for i = 1:length(fieldsrm)
                m2 = rmfield(m2,fieldsrm{i});
              end;
            end;       
            
            save(this.fileunzip,'-struct','m2','-mat', '-v4');
            this = this.readandindex;
            %this = clearemptyfields(this);
            if (this.output);display(['              ... done']);end;
        end;
        
        function [mask,non_zero_index] = getmask(this)
            try
                fa0 = this.m.fa0;
                non_zero_index = (fa0(:) ~= 0.00);
                clear fa0;                
            catch
                try
                    try
                        fa0 = this.m.FAq;
                        non_zero_index = (fa0(:) ~= 0.00);
                        clear fa0;
                    catch
                        try
                            fa0 = this.m.indexdti;
                            non_zero_index = (fa0(:) ~= 0.00);
                            clear fa0;
                        catch
                            fa0 = this.m.fa0orig;
                            non_zero_index = (fa0(:) ~= 0.00);
                            clear fa0;
                        end;
                    end;
                catch
                    fa0 = this.m.fa0;fa1 = this.m.fa1;fa2 = this.m.fa2;
                    non_zero_index = (fa0(:) ~= 0.00 | fa1(:) ~= 0.00 | ...
                        fa2(:) ~= 0.00);
                    clear fa0 fa1 fa2;
                end;
            end;
            mask = zeros(this.dimension)==1;
            mask(non_zero_index) = 1;
        end;

            
        % clear out no longer needed fields of the mat-file
        % when writing new volumes, it might be necessary to clear out some of
        % the fields in the mat-file. Since matlab has no function for this,
        % the rmfield-option no longer works, the code sets the to be deleted
        % fields to empty and this function will clear them
        function this = clearemptyfields(this)
            str = load(this.fileunzip,'-mat');
            fieldsrm = {};
            % which fields to delete
            fie = fields(str);
            for i = 1:length(fie)
              eval(['val=isempty(str.' fie{i} ');']);
              if (val)
                fieldsrm{end+1} = fie{i};
              end;
            end;

            if ~isempty(fieldsrm)
              % the actual deletion
              %str = load(this.fileunzip,'-mat');
              for i = 1:length(fieldsrm)
                str = rmfield(str,fieldsrm{i});
              end;
              save(this.fileunzip,'-struct','str','-mat', '-v4');
              % reread file
              this = readandindex(this);
            else
              clear str;
            end;            
        end;
    end;
end
