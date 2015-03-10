function Afull = STA_system_matrix(xlist,klist,t,tx,b0,mask,varargin)

%%% function Afull = STA_system_matrix(x,k,t,Tx,b0,mask,varargin)
%%% Generates system matrix for STA RF pulse calculation using the Grissom
%%% spatial domain design method (http://dx.doi.org/10.1002/mrm.20978).
%%% This version expects 2D sp-sp designs (i.e. 'space' consists of
%%% x,y,omega and 'k-space' is kx,ky,time

gamma_mT = 267522.119;      % units rad s^-1 mT^-1

% quad = false;
% st_corr = false;
% two_d=false;
% loopcalc=false; % 9-2-11: allow S*A operation to be done in loop
% quad3d=false;
% freq_norev = false; %30-1-13: do not reverse sign of frequency term
% 
% for ii=1:length(varargin)
%     % 7-11-09, quad mode => Nc=1
%     if strcmpi(varargin{ii},'quad')
%         quad = true;
%     end
%     % 3-12-09, staircase corr
%     if strcmpi(varargin{ii},'cor')
%         st_corr = true;
%     end
%     % 28-1-11: force 2d
%     if strcmpi(varargin{ii},'2d')
%         two_d = true;
%     end
%     % 9-2-11: allow S*A operation to be done in loop
%     if strcmpi(varargin{ii},'loop')||strcmpi(varargin{ii},'loopcalc')
%         loopcalc = true;
%     end
%     % 10-1-13: 3D quad mode
%     if strcmpi(varargin{ii},'quad3d')
%         quad3d = true;
%     end
%     % 30-1-13: opt out of frequency axis negation
%     if strcmpi(varargin{ii},'freq_norev')
%         freq_norev = true;
%     end
% end
% No=1;
% switch ndims(tx)
%     case 2
%         [Nx Ny] = size(tx);Nc=1;
%         Nz=1;No=1;
%     case 3
%         if quad
%             [Nx Ny No] = size(tx);
%             Nc=1;
%         else
%             if quad3d
%                 [Nx Ny Nz] = size(tx);
%                 Nc=1;
%             else
%                 [Nx Ny Nc] = size(tx);
%                 Nz=1;No=1;
%             end
%         end
%     case 4
%         if ~two_d %edit 28-1-11
%         [Nx Ny Nz Nc] = size(tx);
%         if (Nc>10)||quad
%             No=Nc;Nc=1;
%         end
% %         % for 3d shim
% %         if (Ny/Nx)>1.5
% %             Nz=1;No=Nc;Nc=1;
% %         end
%         % 24-11-10: edit this, z dim removed for 3d shim
%         if (Ny/Nx)>1.5
%             No=Nz;Nz=1;
%         end
%         else
%             [Nx Ny No Nc] = size(tx);
%             Nz=1;
%         end
%     case 5
%         [Nx Ny Nz No Nc] = size(tx);
% end

[Nx Ny No Nc] = size(tx);


idx = find(mask(:));
Nidx = length(idx);
% generate indices for multicoil variables
N = Nx*Ny*No;

idx_mc = logical(repmat(mask,[1 1 1 Nc]));

% Time axis
dt = t(2)-t(1);
T = t(end);

%% set up b0 - assume in uT
dB0mat = 1e-3 * squeeze(b0(idx))*(t-T); % convert from uT to mT
M0 = 1;
%%% by default use minus sign for frequency: this depends on defition of
%%% phase used in MRI images used to define B0 maps
dB0mat = 1i * gamma_mT * M0 * dt * exp(-1i * gamma_mT * dB0mat);


%%% Fourier Matrix 
F = xlist(idx,:) * klist'; % this is a matrix of dot products
F = exp(1i*F);
A = dB0mat .* F; % element wise multiplication

%%% Finally use Sensitivity information to build system matrix

M = size(A,2);
Afull=zeros([Nidx Nc*M]);
for ii=1:Nc
    %idx1 = (1:Nidx) + (ii-1)*Nidx;
    idx2 = (1:M) + (ii-1)*M;
    tmp=squeeze(tx(idx+(ii-1)*N));
    S = repmat(tmp(:),[1 M]);
    Afull(:,idx2) = S.*A;
end


end