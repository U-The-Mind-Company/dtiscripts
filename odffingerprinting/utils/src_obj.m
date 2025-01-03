% Steven Baete
% NYU LMC CBI
% Februari 2015

classdef src_obj < handle
    % class to hold information about and read from src.gz-files
    % (DSIStudio)
        
    properties
        output;
    end
    
    properties(GetAccess='public', SetAccess='protected')
        m;
        varnamesall;
        varnamesfull;
        varnamesequal;
        
        dimension;
        voxel_size;
        echonum;
        nvol;
        file;
        fileunzip;
    end
    
    methods
        % Constructor
        function this = src_obj(file)
            % the user used a wildcard in the folder-name
            if (length(strfind(file,'*')) > 0)
                file = filestar(file);
            end;
            file = strrep(file,[filesep filesep],[filesep]);

            this.output = true;
            if (this.output);display(['    src_obj - reading ' file]);end;
            % unzip
            if (this.output);display(['     unzipping gz-file']);end;
            if (isunix == 1)
                fileunzip = file(1:(end-3));
                if (~(exist(fileunzip)==2))
                    unix(['gunzip -c ' file ' > ' fileunzip]);
                end
            else
                fileunzip = cell2mat(gunzip(file));
            end;
            
            % what is in our file
            this.m = matfile(fileunzip);
            this.file = file;
            this.fileunzip = fileunzip;
                       
            % read in some parameters
            this.dimension = this.m.dimension;
            this.voxel_size = this.m.voxel_size;
            try
                this.echonum = this.m.echonum;
                this.nvol = length(this.echonum);
            catch
                this.echonum = [];
                this.nvol = size(this.m.b_table,2);
            end;
            
            % make a list of the varnames
            this.varnamesall = whos(this.m);
            % make a list of the available varnames for reading
            leq = 1;lfull = 1;
            for i = 1:length(this.varnamesall)
                varname = this.varnamesall(i).name;
                varsize = this.varnamesall(i).size;
                % fields which need to be copied
                if (varsize(1)*varsize(2) < 7*this.nvol) % such as: b_matrix, b_table, dimension, echonum, q, t+table.voxel_size
                    this.varnamesequal{leq} = varname;leq = leq+1;
                    continue;
                end;
                % fields which are saved fully
                if ((varsize(1) == 1 && varsize(2) == prod(this.dimension)) | (prod(varsize) == prod(this.dimension)))
                    this.varnamesfull{lfull} = varname;lfull = lfull+1;
                    continue;
                end;
                if (this.output);display(['    src_obj - I could not classify variable ' varname ' (size ' num2str(varsize) ') in file ' fileunzip '. Skipping from read!']);end;
            end;
            this.disp;
        end;
        
        % Destructor
        function delete(this)
            if (this.output);display(['    src_obj - deleting unzipped file ' this.file]);end;
            delete(this.fileunzip);
            if (this.output);display(['              ... done']);end;
        end;
        
        function disp(this)
            display(['  src_obj from ' this.file]);
            display(['     dimension [' num2str(this.dimension,'% 4.0f') '], voxel_size [' num2str(this.voxel_size,'% 5.2f') '], nvol [' num2str(this.nvol,'% 5.2f') ']']);
            fprintf(['     available ']);fprintf('%s, ',this.varnamesequal{1:(end-1)});fprintf([this.varnamesequal{end} '\n']);
            fprintf(['               ']);fprintf('%s, ',this.varnamesfull{1:(end-1)});fprintf([this.varnamesfull{end} '\n']);
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
                        %this.montage(S(2).subs{1},S(2).subs{2},S(2).subs{3});
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
            if (strcmp(var,'image'))
                str = load(this.fileunzip,'-mat');
                outvar = zeros([this.dimension,this.nvol]);
                for j = 1:this.nvol
                    varname = [var num2str(j-1)];
                    eval(['outvar(:,:,:,j) = reshape(str.' varname ',this.dimension);']);
                end;
                clear str;
            end;
        end;      
        
        function imagine(this,volume)
            if (volume > 0 && volume < this.nvol)
                imagine(permute(this.getvar(['image' num2str(volume-1)]),[2,1,3]));
            else
                display([' please choose a volume in the range [1,' num2str(this.nvol) ']']);
            end;            
            return
        end; 
                
        function outvar = rawdata(this,ind1,ind2,ind3,ind4)
            if (exist('ind1') && ~strcmp(ind1,':')) 
                ind1 = ind1((ind1 > 0) & (ind1 <= this.dimension(1)));
            else ind1 = 1:this.dimension(1); end;
            if (exist('ind2') && ~strcmp(ind2,':')) 
                ind2 = ind2((ind2 > 0) & (ind2 <= this.dimension(2)));
            else ind2 = 1:this.dimension(2); end;
            if (exist('ind3') && ~strcmp(ind3,':')) 
                ind3 = ind3((ind3 > 0) & (ind3 <= this.dimension(3)));
            else ind3 = 1:this.dimension(3); end;
            if (exist('ind4') && ~strcmp(ind4,':')) 
                ind4 = ind4((ind4 > 0) & (ind4 <= this.nvol));
            else ind4 = 1:this.nvol; end;            
            
            outvar = zeros(length(ind1),length(ind2),length(ind3),length(ind4));
            for i = ind4
                temp = this.getvar(['image' num2str(i-1)]);
                outvar(1:length(ind1),1:length(ind2),1:length(ind3),i) = temp(ind1,ind2,ind3);
            end;            
            return;
        end;
        
        
        function montage(this,volume,slices,msize,sel)
            volume = volume(1);
            if (volume > 0 && volume < this.nvol)
                X = zeros([this.dimension([2,1]),1,this.dimension(3)]);
                X(:,:,1,:) = permute(this.getvar(['image' num2str(volume-1)]),[2,1,3]);
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
                    msize = [NaN,NaN];
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
                montage(X,'Indices',slices,'Size',msize);caxis([0 max(X(:))]);
                ylim([1 msize(1)*size(X,1)]);xlim([1 msize(2)*size(X,2)]);
            else
                display([' please choose a volume in the range [1,' num2str(this.nvol) ']']);
            end;            
            return
        end; 
        
    end;
end
