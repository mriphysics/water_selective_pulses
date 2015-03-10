%%=========================================================================
% 4-2-09: SJM Generate subpulse waveform for binomial sequence

function [pulse,t] = sinc_gauss_subpulse(dur,varargin)

%% set constants and defaults
t0 = dur/2;     % time of pulse centre (?)
ncycles = 5;    % scanner sets this to 2 if 2D
gauss_sd_factor = 1.5;  % scanner sets this to 1.0 for 2D (18-2-09, removed e-3)
dt = 6.4e-6;
t = 0:dt:dur;
t = t - t0;

%% check varargin
for ii=1:nargin-1,
    % ncycles
    if strcmp(varargin{ii},'ncycles')
        ncycles = varargin{ii+1};
        continue;
    end
    % gauss_sd_factor
    if strcmp(varargin{ii},'gauss_sd_factor')||strcmp(varargin{ii},'gf')
        gauss_sd_factor = varargin{ii+1};
        continue;
    end   
    % dt
    if strcmp(varargin{ii},'dt')
        dt = varargin{ii+1};
        t = 0:dt:dur;
        t = t - t0;
        continue;
    end
    % specify time for ramp sampling
    if strcmp(varargin{ii},'t')
        t = varargin{ii+1};
        continue;
    end
end

%% calculate RF samples
K = 2 * (t0/gauss_sd_factor)^2; % gaussian time constant
wt = 2*pi*ncycles*t/dur;

pulse = ( sin(wt)./wt ) .* exp(-t.^2/K);

% if t==0 get NaN
pulse(isnan(pulse))=1;
end