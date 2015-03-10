%%% This script can be used for optimizing binomial pulses with parallel
%%% transmit MRI as outlined in Malik et al 2009 (http://dx.doi.org/10.1002/mrm.22260)

%%% In order to understand this code please read the paper. The code is
%%% provided along with some example data and running this script will
%%% generate all of the results.

%%% Copyright Shaihan Malik, 2015. shaihan.malik@gmail.com
%%% This code is provided under the terms of the MIT license. For details
%%% please see license.txt. If you use this code please cite the reference
%%% listed above.



%% constant definitions
gamma_mT = 267522.119;      % units rad s^-1 mT^-1
dt_system = 6.4e-6;         % System RF dwell time
f_fat = -434;               % Chemical shift of fat in Hz @3T
t_wf = 1/(2*abs(f_fat));    % Time for water and fat to be out of phase

%%% Design parameters
Nsubpulse = 4;              % Number of subpulses in SpSp pulse (default=4)
tgap = t_wf;                % Gap between subpulses
flip = 20;                  % total flip angle, degrees
weighting_factor = 1e-3;    % weighting of water wrt fat

%%% Pulse specific parameters
slab_thick = 0.2;           % slab thickness, metres
slab_offc = 0.0;            % slab offcentre
TBP = 4;                    % time bandwidth product of subpulse
gauss_factor = 1.5;         % gauss weighting in s-g pulse
slew = 100e3;               % slew rate, mT/m/s
% gradient_orientation = [0 0 1]; % gradient orientation vector


addpath(genpath('lib'))
addpath(genpath('bin'))


%% Load in B1 and B0 maps and define FoX

%%% FOX
x.dx = 8e-3;
x.dy = 8e-3;
x.FOV = [-0.1895 0.1895 -0.1895 0.1895];

% spectral domain
omega = 2*pi* (-600:20:600);

% Overall spatial domain is x,y,omega
[x.x,x.y,x.o] = ndgrid(x.FOV(1):x.dx:x.FOV(2),x.FOV(3):x.dy:x.FOV(4),omega);
[Nx,Ny,No] = size(x.x);N = numel(x.x);
xlist = cat(2,x.x(:),x.y(:),x.o(:)); % Nx3 column matrix


%%% --- b0, b1 maps ---
Nc = 8;

% Load in data:
% contains tx=transmit sensitivity (relative units), b0 (uT) and a mask
load b1_b0_data 

%%% Define specific bands for optimization and masks
maskf=zeros([Nx Ny No]);
water_band = find(((omega/(2*pi))>-30)&(omega/(2*pi))<30);
fat_band = find(((omega/(2*pi))<-390)&(omega/(2*pi))>-450);
maskf(:,:,water_band)=1;
maskf(:,:,fat_band)=1;


%%% get indices of mask
mask = maskf .* m;
w = find(mask(:));
wxy = find(m(:));
Nw=length(w);
wc = find(repmat(mask,[1 1 1 Nc]));

%%% weighting matrix 
W = mask;
W(:,:,water_band)=W(:,:,water_band) * weighting_factor;


%% k-space 

% timing term
t = tgap*(0:Nsubpulse-1);
T=t(end);

kx=zeros([Nsubpulse 1]);
ky=kx;
[k.x,k.y,k.o] = ndgrid(0,0,t-T);


klist = cat(2,kx(:),ky(:),k.o(:)); % Mx3 column matrix


%% Based on these grids, k-spaces and field maps, construct system matrices 
% for optimization

%%% Full system matrix for all frequencies
Afull = STA_system_matrix(xlist,klist,t,tx,b0,m);

%%% Reduced matrix for only frequencies of interest
Ared = STA_system_matrix(xlist,klist,t,tx,b0,mask);

%% Forward model of Binomial sp-sp pulse

%%% Generate Binomial sequence
pulse_amps = zeros([Nsubpulse 1]);
for ii=1:Nsubpulse
    pulse_amps(ii) = nchoosek(Nsubpulse-1,ii-1);
end
phase_shift = -(1-(tgap/t_wf))*pi;
pulse_amps = pulse_amps .* exp(1i*(0:(Nsubpulse-1))'*phase_shift);
b = repmat(pulse_amps(:),[1 Nc]);

%%% scaling: we want the nominal flip angle of total to be as specified
%%% above (in radians)
b = b * (flip*pi/180)*((1/sum(abs(pulse_amps))) / (gamma_mT * tgap));


P_w = Afull*b(:);
P = zeros([Nx Ny No]);
P(wxy)=P_w;


%% Now make an optimized excitation

%%% Specify target
P_target = zeros(size(mask));
P_target(:,:,water_band)=flip*pi/180; 
P_target(:,:,fat_band)=0;


%%% MLS solution
% apply weighting
ww = diag(W(w));
ww=sqrt(ww); 
wA = ww*Ared;
P_target_w = ww*P_target(w);
clear ww

[b_mls,z2,res] = MLS(wA,P_target_w,0,1,1000,1e-10);

%  time domain samples of RF pulse 
bt_mls = reshape(b_mls,[Nsubpulse Nc]);

%%% Map this back to whole frequency domain
pt_mls_w = Afull * b_mls;
Pmls = zeros([Nx Ny No]);
Pmls(wxy) = pt_mls_w;

%% Display some figures

figure;
nr=2;nc=2;
subplot(nr,nc,1);
imagesc(abs(P(:,:,water_band(1)))*180/pi,[0.2 1.2]*flip);
title('Water excitation (binomial), degrees')
colorbar
axis off

subplot(nr,nc,2);
imagesc(abs(P(:,:,fat_band(1)))*180/pi,[0 0.2]*flip);
title('Fat excitation (binomial), degrees')
colorbar
axis off
hold on

subplot(nr,nc,nc+1);
imagesc(abs(Pmls(:,:,water_band(1)))*180/pi,[0.2 1.2]*flip);
title('Water excitation (binomial), degrees')
colorbar
axis off

subplot(nr,nc,nc+2);
imagesc(abs(Pmls(:,:,fat_band(1)))*180/pi,[0 0.2]*flip);
title('Fat excitation (binomial), degrees')
colorbar
axis off

colormap gray
set(gcf,'position',[100 100 750 550])

%%% Pulse amplitudes

figure
subplot(2,1,1)
plot(abs(b)*gamma_mT*tgap*180/pi,'k','linewidth',2)
hold on
grid on
plot(abs(bt_mls)*gamma_mT*tgap*180/pi)
set(gca,'xtick',1:Nsubpulse)
xlim([0.5 Nsubpulse+0.5])
xlabel('Subpulse number')
ylabel('Nominal Flip angle, deg')
title('Optimized amplitudes (binomial in heavy black)')

subplot(2,1,2)
plot(angle(b)*180/pi,'k','linewidth',2)
hold on
grid on
plot(angle(bt_mls)*180/pi)
set(gca,'xtick',1:Nsubpulse)
xlim([0.5 Nsubpulse+0.5])
xlabel('Subpulse number')
ylabel('Phase, deg')
title('Optimized Phases (binomial in heavy black)')

set(gcf,'position',[100 100 450 800])


%% Interactive explore frequency response. Ctrl-C to exit
vox_x = [10 40];
vox_y = [30 15];
figure;
set(gcf,'position',[100 200 800 400])
while 1
subplot(1,2,1)
imagesc(b0(:,:,1) * 1e-3*gamma_mT /(2*pi));
axis off
colormap gray
title('B0/Hz Click to see frequency response')
drawnow

[gx,gy]=ginput(1);
gx=fix(gx);gy=fix(gy);

subplot(1,2,2)
cla
hold on
grid on
plot(omega/(2*pi),squeeze(abs(P(gy,gx,:)))*180/pi);
plot(omega/(2*pi),squeeze(abs(Pmls(gy,gx,:)))*180/pi);
xlim([min(omega/(2*pi)) max(omega/(2*pi))])
legend('Binomial','Optimized')
title('Frequency responses. Green=fat Red=water')
xlabel('f/Hz')
ylabel('flip/degrees')
ylim([0 flip])

p1 = patch(omega([water_band fliplr(water_band)])/(2*pi),[zeros(size(water_band)) flip*ones(size(water_band))],[1 0 0]);
p2 = patch(omega([fat_band fliplr(fat_band)])/(2*pi),[zeros(size(fat_band)) flip*ones(size(fat_band))],[0 1 0]);
set([p1 p2],'facealpha',0.2,'edgealpha',0)

end

%% Generate pulse sequence

flips_mls = (bt_mls) * gamma_mT * tgap;

%%% The function gen_binomial_sequence is old (circa 2009) and has not been
%%% tested extensively. It does however generate fly-back sp-sp pulses in a
%%% form that may be useful to someone so please feel free to use it
[pulse_out,g,timing] = gen_binomial_sequence(flips_mls,'thk',slab_thick,'pl',...
    'ncycle',TBP,'slew',slew,'ori',[0 0 1],'dz',slab_offc,...
    'b1',10e-3,'flip',flip*pi/180);
