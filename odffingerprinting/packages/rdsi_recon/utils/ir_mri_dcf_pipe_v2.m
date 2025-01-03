% ir_mri_dcf_pipe()
% Pipe&Menon, based on equalty 1 = A A' w so w = w ./ (A A' w)
%
% Steven Baete
% NYU SOM CBI
% March 2019
% based on ir_mri_dcf_pipe from the irt toolbox of Jeff Fessler
% SB: added filter function Psi as in the original Pipe&Menon 1999 paper
% (Pipe199b.pdf)
% Pipe et al, MRM 41:179-186, 1999
% the function used for psi is an approximation, not sure why the iteration
% doesn't converge faster
function wi = ir_mri_dcf_pipe_v2(kspace, G, varargin)
arg.isave = [];
arg.niter = 50;
arg.thresh = 0.01;
arg.fov = 1;
arg.unitv = [];
arg.wi = []; % see below
arg = vararg_pair(arg, varargin, 'allow_new',true);

% P = G.arg.st;
if (isstruct(G))
    P = G.arg.st.p;
    Dims = length(G.arg.st.Kd);
    N = G.arg.st.Kd(1);
    if (Dims > 2)
        NZ = G.arg.st.Kd(3);
    end;
else
    P = G;
    Dims = length(size(arg.Nsize));
    N = arg.Nsize(1);
    if (Dims > 2)
        NZ = arg.Nsize(3);
    end;
end;
goal = inf;
iter = 0;
%saver = zeros(arg.niter,1);

beta = 16;
order = 8;
Npsi = 1024;
Psi1D = mkaiserwin(Npsi,beta,order);   %k-space
% psi1D = fftshift((fft(Psi1D))); %image
% the /9 is closer to the Pipe-paper, but the line below gives a better
% balanced result
%kind = ((-Npsi/2):(Npsi/2-1))/(Npsi/9);
kind = ((-Npsi/2):(Npsi/2-1))/(Npsi/32);
% rind = (kind(2)-kind(1))*((-Npsi/2):(Npsi/2-1))*4;

% Psi1D, |k| > 2 = 0; psi1D, |r| > 1 = 0
% Psi1D(abs(kind)>2) = 0;

% for i = 1:10  
%     Psi1Dold = Psi1D;
%     psi1D(abs(rind)>1) = 0;
%     Psi1D = ifft((fftshift(psi1D)));
%     Psi1D(abs(kind)>2) = 0;
%     psi1D = fftshift((fft(Psi1D)));
%     figure;
%     subplot(3,1,1);plot(kind,Psi1D);xlim([-2 2]);
%     subplot(3,1,2);plot(rind,abs(psi1D));gg = gca;gg.YScale = 'log';xlim([-2 2]);
%     subplot(3,1,3);plot(kind,Psi1D-Psi1Dold);xlim([-2 2]);
%     title(num2str(i));
% end;

% N = arg.Nsize;%
% N = G.arg.st.Kd(1);
if (Dims == 2)
    [X,Y] = meshgrid(-N/2:(N/2-1),-N/2:(N/2-1));
    kdist = sqrt(X.^2+Y.^2);
    Psi = reshape(interp1(kind((end/2+1):end),...
                    Psi1D((end/2+1):end),...
                    kdist(:)),[N,N]);
    %Psi(abs(kdist)>2) = 0;
    Psi(isnan(Psi(:))) = 0;
    psi = abs(ifft2c(Psi));
else
    [X,Y,Z] = meshgrid(-N/2:(N/2-1),-N/2:(N/2-1),-NZ/2:(NZ/2-1));
    kdist = sqrt(X.^2+Y.^2+Z.^2);
    Psi = reshape(interp1(kind((end/2+1):end),...
                    Psi1D((end/2+1):end),...
                    kdist(:)),[N,N,NZ]);
    %Psi(abs(kdist)>2) = 0;
    Psi(isnan(Psi(:))) = 0;
    psi = abs(ifft3c(Psi));
end;
psi = psi/sum(abs(psi(:)))*numel(psi);

%figure;plot(Psi1D);
%figure;plot(psi1D((end/2+1):end));gg = gca;gg.YScale = 'log';
%figure;imagesc(Psi)
%figure;plot(abs(Psi(:,end/2+1)));
% figure;plot(abs(psi(:,end/2+1)));gg = gca;gg.YScale = 'log';

%{
%removed 2017-04-11 in favor
scale = G.arg.st.sn(end/2,end/2)^(-2) / arg.fov^2 ...
	/ prod(G.arg.st.Kd) * prod(G.arg.st.Nd);
scale = reale(scale);
%}

if isempty(arg.wi)
	wi = ones(size(kspace,1), 1); % default initial
    wi = wi/sum(wi(:));
else
%	minmax(arg.wi)
%	wi = arg.wi / scale; % todo: why?
	wi = arg.wi;
end
wi_save = wi;
while(max(abs(goal-1)) > arg.thresh)
	iter = iter + 1;
	if (isstruct(G) & isfield(G.arg.st,'interp_table'))
		goal = G.arg.st.interp_table(G.arg.st, ...
			G.arg.st.interp_table_adj(G.arg.st, wi) );
    else
        wpsf = P' * wi; %image space (+- psf)
        if (~strcmp(class(P),'NUFFT'))
            wpsi = wpsf.*psi(:); % image space: conv -> mult
        else
            wpsi = wpsf.*psi;
        end;
		goal = P * wpsi;
	end
% 	wi = wi ./ real(goal);
	%wi = wi ./ abs(goal);
    wi = wi ./ goal;
	if iter > arg.niter
		%warn 'iteration stuck?'
		break
	end
%	saver(iter) = max(abs(goal-1));
	if ~isempty(arg.isave)
		wi_save = [wi_save, wi];
    end
    %display([' iter ' num2str(iter) ' goal ' num2str(max(abs(goal-1)))]);
end
%printm('pipe ended at iteration %d with %g', iter, max(abs(goal-1)))
%plot(saver(2:end))

wi = real(wi);

% scale
if isempty(arg.unitv)
	%warn 'no unitv so no scaling - user must fix scaling'
else
	fail 'todo: this needs tested1'
	last = wi_save(:,1);
	psf_raw = G' * last;
	e0 = ig.unitv % todo
	scale = dot(e0, psf_raw) / norm(psf_raw(:)).^2;
	wi_save = scale * wi_save;
end

if ~isempty(arg.isave)
	wi = wi_save;
end

%wi = wi / scale; % todo: why?
