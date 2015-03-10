%%=========================================================================
% 4-2-09: SJM Generate pulse and gradient waveforms for binomial sequence
% 23-3-09: allow sequence with no flyback gradient

function [rfpulse,g,timing] = gen_binomial_sequence(weights,varargin)

%% --- Define Constants and parameters ---
gamma_mT = 267522.1;    % rad s^-1 mT^-1
dt_system = 6.4e-6;     % sample rate
dt = 0.1e-6;            % higher sample rate, interpolate down afterwards
max_gr = 10;            % mT m^-1
max_slew = 16.667e3; % mT m^-1 s^-1

[Npulse Nc] = size(weights);    % weights are Npulse * Ncoil
% weights are scaled s.t. each is nominal
% flip in radians
ori = [0 0 1];      % orientation vector
slab_thick = 0.3;   % slab thickness, metres
z_offset = 0.0;     % slice offset

flip_max_sp = max(abs(weights(:)));% maximum flip angle from one sub pulse
define_flip = false;
max_b1 = 10e-3;     % 10uT
ncycle = 5;         % # sinc pulse cycles, scanner sets this to 2 if 2D
t_wf = 1.1513e-3;   %water-fat precession time (sec)
t_gap = t_wf;
gauss_sd_factor = 1.5; % Gaussian envelope factor for subpulse pulse

g_times_dur = 4*pi*ncycle / (gamma_mT * slab_thick); % gradient time product

balanced = false;   % balanced => replicate rephase grad at start
gen_plots = false;  % produce a plot of the gradient and RF waveforms
maxdur = false;     % aim to make a maximum duration solution by increasing the slew rate
flyback = true;  % use flyback gradients as default 

%% get extra inputs
for ii=1:nargin-1,
    % slab thickness
    if strcmp(varargin{ii},'thk')
        slab_thick = varargin{ii+1};
        continue;
    end
    % slab thickness
    if strcmp(varargin{ii},'dz')||strcmp(varargin{ii},'z_offset')...
            ||strcmp(varargin{ii},'dx')||strcmp(varargin{ii},'dy')
        z_offset = varargin{ii+1};
        continue;
    end
    % number of sinc cycles
    if strcmp(varargin{ii},'ncycle')
        ncycle = varargin{ii+1};
        continue;
    end
    % gaussian envelope for sinc-gauss pulses
    if strcmp(varargin{ii},'gauss_sd_factor')||strcmp(varargin{ii},'gf')
        gauss_sd_factor = varargin{ii+1};
        continue;
    end   
    % orientation
    if strcmp(varargin{ii},'ori')
        ori = varargin{ii+1};
        ori = ori/norm(ori); % normalise orientation vector
        continue;
    end   
    % dt
    if strcmp(varargin{ii},'dt')
        dt_system = varargin{ii+1};
        continue;
    end
    % balanced
    if strcmp(varargin{ii},'balanced')||strcmp(varargin{ii},'bal')
        balanced=true;
        continue;
    end
    % plots
    if strcmp(varargin{ii},'pl')||strcmp(varargin{ii},'plots')
        gen_plots=true;
        continue;
    end
    % max B1 (mT)
    if strcmp(varargin{ii},'b1')||strcmp(varargin{ii},'B1')
        max_b1 = varargin{ii+1};
        continue;
    end
    % user specified max flip angle (rad)
    if strcmp(varargin{ii},'flip')||strcmp(varargin{ii},'flipmaxflip')
        flip_max_sp = varargin{ii+1};
        fprintf(1,'gen_binomial_seq: \nMax flip specified is %2.1f deg; ref B1 corresponds to %2.1f deg\n',...
            max(abs(weights(:)))*180/pi,flip_max_sp*180/pi);
        define_flip=true;
        continue;
    end
    % max G (mT/m)
    if strcmp(varargin{ii},'G')||strcmp(varargin{ii},'gr')
        max_gr = varargin{ii+1};
        continue;
    end
    % max slew (mT/m/s)
    if strcmp(varargin{ii},'S')||strcmp(varargin{ii},'slew')
        max_slew = varargin{ii+1};
        continue;
    end
    % gap between pulses
    if strcmp(varargin{ii},'tgap')||strcmp(varargin{ii},'gap')
        t_gap = varargin{ii+1};
        continue;
    end
    % drives to a maximum slew solution which gives maximum sub pulse duration
    if strcmp(varargin{ii},'maxdur')
        maxdur=true;
        continue;
    end
    % transmit on negative gradient lobes
    if strcmp(varargin{ii},'nofly')||strcmp(varargin{ii},'noflyback')
        flyback=false;
        fprintf(1,'gen_binomial_seq: not using flyback gradients\n');
        continue;
    end
end
% may have changed:
g_times_dur = 4*pi*ncycle / (gamma_mT * slab_thick);

%% initial definition of sub-pulses
% start with max gradient => min duration but max peak B1
grad = max_gr;
% calculate sub pulse duration based on the slab width
dur = g_times_dur / grad;
subpulse = sinc_gauss_subpulse(dur,'dt',dt,'ncycles',ncycle,'gf',gauss_sd_factor);
% B1  integral of pulse (this is the effective B1 integral, no units
% i.e. flip (rad) = gamma * sp_b1_int * sp_b1 where sp_b1 is peak b1 of
% strongest subpulse
sp_b1_int = sum(subpulse) * dt;
sp_b1 = (flip_max_sp) / (gamma_mT * sp_b1_int);


% first check B1
if sp_b1 > max_b1
    % increase duration by factor sp_b1/max_b1 
    dur = (sp_b1/max_b1) * dur;
    grad = g_times_dur / dur; % this is now the maximum acceptable gradient
    subpulse = sinc_gauss_subpulse(dur,'dt',dt,'ncycles',ncycle,'gf',gauss_sd_factor);
    sp_b1_int = sum(subpulse) * dt;
    sp_b1 = (flip_max_sp) / (gamma_mT * sp_b1_int);
end

% now calculate minimum required slew rate
if flyback
    slew = 4*grad*t_gap/(t_gap-dur)^2;
else
    slew = 2*grad/(t_gap-dur);
end

%% 23-3-09: following calculation depends on whether we use flyback gradients or not
if flyback
    
% check that the peak rephase gradient is OK
dur_s_r = 0.5*((t_gap-dur) - 2*grad/slew);
grad_r = slew * dur_s_r/2;

if (slew > max_slew)||maxdur||(dur_s_r<0) 
    % go to max slew
    slew=max_slew;
    % can we reduce the gradient to make this acceptable?
    a = 4*t_gap/max_slew;
    b = -t_gap^2;
    c = 2*g_times_dur*t_gap;
    d = -g_times_dur^2;
    grad_poss = roots([a b c d]);
    % is there a real root?
%     discr = -4*b^3*d + b^2*c^2 - 4*a*c^3 + 18*a*b*c*d - 27*a^2*d^2;

    % check all roots
    found_solution=false;
    for ii=3:-1:1,
%         disp(ii)
        if ~isreal(grad_poss(ii))
            continue; % cplx root
        end
        % check that the resulting rephase duration/strength are acceptable
        dur_poss = g_times_dur / grad_poss(ii);
        dur_s_r = 0.5*((t_gap-dur_poss) - 2*grad_poss(ii)/slew);
        grad_r = slew * dur_s_r/2;
        dur_r=0;

        if (dur_s_r<0)||(grad_r<0)||(grad_r>max_gr)||(grad_poss(ii)>grad)
            continue;
        else
            grad = grad_poss(ii);
            dur = dur_poss;
            subpulse = sinc_gauss_subpulse(dur,'dt',dt,'ncycles',ncycle,'gf',gauss_sd_factor);
            sp_b1_int = sum(subpulse) * dt;
            sp_b1 = (flip_max_sp) / (gamma_mT * sp_b1_int);
            slew = 4*grad*t_gap/(t_gap-dur)^2;
            fprintf(1,'Found slew limited solution\n');
            fprintf(1,'Max B1 %3.2f uT, grad %3.2f mT/m, slew %5.3f mT/m/s, dur %3.3f ms, grad_r %3.2f mT/m, dur_s_r %3.3f ms\n',...
                sp_b1*1e3,grad,slew,dur*1e3,grad_r,dur_s_r*1e3);
            found_solution=true;
            ramp_sample = false;
            break
        end
    end
    if ~found_solution
        fprintf(1,'No solution is possible without ramp sampling ...\n');
%         rfpulse=[];g=[];timing=[];
%         return
        ramp_sample = true;
    end
else

    fprintf(1,'Max B1 %3.2f uT, grad %3.2f mT/m, slew %5.3f mT/m/s, dur %3.3f ms, grad_r %3.2f mT/m, dur_s_r %3.3f ms\n',...
                sp_b1*1e3,grad,slew,dur*1e3,grad_r,dur_s_r*1e3);
    ramp_sample = false;
end
   
% another step might be to now try increasing the slew rate so that the
% negative lobes can be trapezoidal and not triangular...

%==============================
% ramp sampling bit for flyback
%==============================
if ramp_sample
    % do above again but allow f to increase. f = dur_rf/dur
    found_solution=false;
    f=1;fmax=2;
    while (~found_solution)&&(f<fmax)
        slew=max_slew;
        % calculte possible gradients
        a = 4*t_gap/max_slew;
        b = -t_gap^2;
        c = 2*(g_times_dur/f)*t_gap;
        d = -(g_times_dur/f)^2;
        grad_poss = roots([a b c d]);


        for ii=3:-1:1,
            if ~isreal(grad_poss(ii))
                continue; % cplx root
            end
            % check that the resulting rephase duration/strength are acceptable
            dur_poss = (g_times_dur/f) / grad_poss(ii);
            dur_s_r = 0.5*((t_gap-dur_poss) - 2*grad_poss(ii)/slew);
            grad_r = slew * dur_s_r/2;
            dur_r=0;

            if (dur_s_r<0)||(grad_r<0)||(grad_r>max_gr)||(grad_poss(ii)>grad)
                continue;
            else
                grad = grad_poss(ii);
                dur = dur_poss;
                subpulse = sinc_gauss_subpulse(dur*f,'dt',dt,'ncycles',ncycle,'gf',gauss_sd_factor);
                sp_b1_int = sum(subpulse) * dt;
                sp_b1 = (flip_max_sp) / (gamma_mT * sp_b1_int);
                slew = 4*grad*t_gap/(t_gap-dur)^2;
                fprintf(1,'Found ramp sample solution with f = %1.2f\n',f);
                fprintf(1,'Max B1 %3.2f uT, grad %3.2f mT/m, slew %5.3f mT/m/s, dur %3.3f ms, grad_r %3.2f mT/m, dur_s_r %3.3f ms\n',...
                    sp_b1*1e3,grad,slew,dur*1e3,grad_r,dur_s_r*1e3);
                found_solution=true;
                break
            end
        end
        
        if ~found_solution
            f = f + 0.01;
        end
    end
    if ~found_solution
        fprintf(1,'Ramp sampling will not work, try increasing slew rate or slab thickness\n')
        rfpulse=[];g=[];timing=[];
        return
    end

    if f>(1+2*grad/(slew*dur))
        fprintf(1,'Ramp sample solution will not work: sub pulses extend beyond positive gradient lobes\n')
        fprintf(1,'Try again with increased slab thickness or higher maximum slew rate\n')
        rfpulse=[];g=[];timing=[];
        return
    end

else
    f=1;
end

end

%% no flyback solution
if ~flyback
   
if (slew > max_slew)||maxdur||(dur>t_gap/2)
    % go to max slew solution
    slew=max_slew;
    % calculate acceptable gradient value
    a = 2;
    b = -slew*t_gap;
    c = g_times_dur*slew;
    grad_poss = roots([a b c]);
    
    found_solution=false;
    for ii=2:-1:1,
%         disp(ii)
        if ~isreal(grad_poss(ii))
            % cplx root => Not possible with any gradient
            break;
        end
        if grad_poss(ii)<0
            continue; % -ve root
        end
        if grad_poss(ii)>max_gr
            continue; % exceeds maximum
        end
        % check that the resulting duration is acceptable (no ramp
        % sampling)
        dur_poss = g_times_dur / grad_poss(ii);
        if dur_poss <= (t_gap - 2*grad_poss(ii)/slew)
            %  accept the solution (I don't know why this condition is
            %  violated but it is sometimes)
            grad = grad_poss(ii);
            dur = dur_poss;
            subpulse = sinc_gauss_subpulse(dur,'dt',dt,'ncycles',ncycle,'gf',gauss_sd_factor);
            sp_b1_int = sum(subpulse) * dt;
            sp_b1 = (flip_max_sp) / (gamma_mT * sp_b1_int);
            slew = 2*grad/(t_gap-dur);
            fprintf(1,'Found slew limited solution w/o flyback\n');
            fprintf(1,'Max B1 %3.2f uT, grad %3.2f mT/m, slew %5.3f mT/m/s, dur %3.3f ms\n',...
                sp_b1*1e3,grad,slew,dur*1e3);
            found_solution=true;
            ramp_sample = false;
            break
        end
    end
    
    if ~found_solution
        fprintf(1,'No solution is possible without ramp sampling, this is not implemented for NO FLYBACK ...\n');
%         rfpulse=[];g=[];timing=[];
        ramp_sample=true;
%         return
    end
else
        fprintf(1,'Max B1 %3.2f uT, grad %3.2f mT/m, slew %5.3f mT/m/s, dur %3.3f ms\n',...
                sp_b1*1e3,grad,slew,dur*1e3);
            ramp_sample = false;
end



%==================================
% ramp sampling bit for no flyback 
%==================================
if ramp_sample
    % do above again but allow f to increase. f = dur_rf/dur
    found_solution=false;
    f=0.5;fmax=2;
    while (~found_solution)&&(f<fmax)
        slew=max_slew;
        % calculte possible gradients
        a = 2;
        b = -slew*t_gap;
        c = (g_times_dur/f)*slew;
        grad_poss = roots([a b c]);

        for ii=2:-1:1,
            if ~isreal(grad_poss(ii))
                break; % cplx root
            end
            if grad_poss(ii)<0
                continue; % -ve root
            end
            if grad_poss(ii)>max_gr
                continue; % exceeds maximum
            end
            
            % check that the resulting duration is acceptable
            dur_poss = (g_times_dur/f) / grad_poss(ii);
            if dur_poss <= (t_gap - 2*grad_poss(ii)/slew)
                %  accept the solution (I don't know why this condition is
                %  violated but it is sometimes)
                grad = grad_poss(ii);
                dur = dur_poss;
                subpulse = sinc_gauss_subpulse(dur*f,'dt',dt,'ncycles',ncycle,'gf',gauss_sd_factor);
                sp_b1_int = sum(subpulse) * dt;
                sp_b1 = (flip_max_sp) / (gamma_mT * sp_b1_int);
                slew = 2*grad/(t_gap-dur);
                fprintf(1,'Found slew limited solution w/ flyback, f = %1.2f\n',f);
                fprintf(1,'Max B1 %3.2f uT, grad %3.2f mT/m, slew %5.3f mT/m/s, dur %3.3f ms\n',...
                    sp_b1*1e3,grad,slew,dur*1e3);
                found_solution=true;
                break
            end
        
        end
        
        if ~found_solution
            f = f + 0.01;
        end
    end
    if ~found_solution
        fprintf(1,'Ramp sampling will not work, try increasing slew rate or slab thickness\n')
        rfpulse=[];g=[];timing=[];
        return
    end

    if f>(1+2*grad/(slew*dur))
        fprintf(1,'Ramp sample solution will not work: sub pulses extend beyond positive gradient lobes\n')
        fprintf(1,'Try again with increased slab thickness or higher maximum slew rate\n')
        rfpulse=[];g=[];timing=[];
        return
    end

else
    f=1;
end


end

%% slab offset (19-2-09)
Nsp= length(subpulse);
if abs(z_offset) > 0
    dw = gamma_mT * grad * z_offset;
    % phase ramp for each subpulse
    phase_ramp = exp(-i*(1:Nsp)*dt*dw);
else
    phase_ramp = ones([1 Nsp]);
end
   
%% --- Build gradient waveform ----
% positive lobes
Nramp = ceil((grad/slew)/dt);
Nflat = ceil(dur/dt);
g_ramp = (0:(Nramp-1))*slew*dt;
g_flat = g_ramp(end) * ones([1 Nflat]);
g_positive = [g_ramp,g_flat,fliplr(g_ramp)];

if flyback
    % negative lobes ASSUME TRIANGULAR FOR NOW
    Nramp_r = ceil(dur_s_r/dt);
    g_ramp_r1 = (1:(Nramp_r-1))*slew*dt; % start at 1 so concatenation is smooth
    g_ramp_r2 = (Nramp_r-1:-1:1)*slew*dt;
    g_negative = -[g_ramp_r1,g_ramp_r2];
else
    g_negative = -g_positive;
end

% build into train
if flyback
    g_core_samples = [g_positive repmat([g_negative g_positive],[1 Npulse-1])];
else
    if mod(Npulse,2)==0
        g_core_samples = repmat([g_positive g_negative],[1 Npulse/2]);
    else
        g_core_samples = [repmat([g_positive g_negative],[1 (Npulse-1)/2]) g_positive];
    end
end
    


%=============================================
% Gradient rewinder, do this @max slew if poss
rewind_area = 0.5*sum(g_positive)*dt;%0.5*(grad^2 / slew + grad*dur);
% max_slew=40e3;
grad_r = sqrt(rewind_area * max_slew);
dur_r = grad_r/max_slew;
slew_r=max_slew;

% check if Grad is too high
if grad_r <= max_gr
    % OK, use triangular blip
    Nramp_r = ceil(dur_r/dt);
    g_ramp_r1 = (1:(Nramp_r-1))*slew_r*dt;
    g_ramp_r2 = (Nramp_r:-1:1)*slew_r*dt;
    g_rewind_lobe = [g_ramp_r1,g_ramp_r2];
else
    % Not OK, use max grad instead with trapezoidal lobe
    grad_r = max_gr;
    slew_r = max_slew;
    dur_r = (rewind_area - grad_r^2/slew_r)/grad_r;
    dur_r_ramp = grad_r/slew_r;
    Nramp_r = ceil(dur_r_ramp/dt);
    Nflat_r = ceil(dur_r/dt);
    g_ramp_r = (0:(Nramp_r-1))*slew_r*dt;
    g_flat_r = g_ramp_r(end) * ones([1 Nflat_r]);
    g_rewind_lobe = [g_ramp_r,g_flat_r,fliplr(g_ramp_r)];
end

% Polarity of rewind lobe
if flyback
    % always negative
    g_trail_samples = -g_rewind_lobe;
else
    if mod(Npulse,2)==0
        g_trail_samples = g_rewind_lobe;
    else
        g_trail_samples = -g_rewind_lobe;
    end
end

% if balanced add a (p)rephase gradient
if balanced
    g_lead_samples = g_trail_samples;
else
    g_lead_samples = 0; % simple leading zeros
end


%% --- Build RF waveform ----
Np = length(g_positive);
Nn = length(g_negative);
Nsp= length(subpulse);
% make an array of indices to contain RF samples
if flyback
    for ii=1:Npulse,
        pulse_centre(ii) = fix(Np/2) + (ii-1) * (Np+Nn);
        pulse_sample{ii} = (1:Nsp) - fix(Nsp/2) + pulse_centre(ii);
    end
else
    for ii=1:Npulse,
        pulse_centre(ii) = fix(Np/2) + (ii-1) * (Np);
        pulse_sample{ii} = (1:Nsp) - fix(Nsp/2) + pulse_centre(ii);
    end
end
if ramp_sample
    % get gradient samples during which the pulse is on
    g_rf_ramp = g_positive((1:Nsp) - fix(Nsp/2) + pulse_centre(1));

    for ii=1:length(g_rf_ramp),tt(ii) = -sum(g_rf_ramp(ii:end))*dt;end
    tp = (tt)/grad;
%     tp = (1:length(g_rf_ramp))*dt;
    tp = tp - tp(fix(length(tp)/2));
    % get new pulse
    subpulse = sinc_gauss_subpulse(dur*f,'dt',dt,'ncycles',ncycle,'gf',...
                        gauss_sd_factor,'t',tp);
    % density compemsation ...
    subpulse = subpulse .* (g_rf_ramp/grad);
end    

% add phase ramp to subpulse
subpulse = subpulse .* phase_ramp;

% now make RF samples array
fullpulse = zeros([length(g_core_samples) Nc]);
% SJM11-3-09: replace line below with one after
if ~define_flip
    rfweights = sp_b1 * weights / max(abs(weights(:)));
else
    rfweights = sp_b1 * weights / flip_max_sp;
    fprintf(1,'gen_binomial_seq: \nActual peak B1 is %2.2f uT, not %2.2f uT\n',...
        max(abs(rfweights(:)))*1e3,sp_b1*1e3);
end
for ii=1:Npulse,
    for jj=1:Nc,
        fullpulse(pulse_sample{ii},jj) = rfweights(ii,jj) .* subpulse;
    end
end

%% resample all RF and gradients for output
% do core samples first
ti = (0:(length(g_core_samples)-1)) .* dt;
to = 0:dt_system:ti(end);
rfpulse = interp1(ti,fullpulse,to);
if Nc==1
    rfpulse = rfpulse.';
end
g_core_out = interp1(ti,g_core_samples,to);

% do same for lead and trail just for gradients
if balanced
    ti = (0:(length(g_lead_samples)-1)) .* dt;
    to = 0:dt_system:ti(end);
    g_lead_out = interp1(ti,g_lead_samples,to);
else
    g_lead_out = g_lead_samples; % just a zero
end
ti = (0:(length(g_trail_samples)-1)) .* dt;
to = 0:dt_system:ti(end);
g_trail_out = interp1(ti,g_trail_samples,to);

%% ---- variables for output ---------------
g.Nlead = numel(g_lead_out);
g.Ncore = numel(g_core_out);
g.Ntrail = numel(g_trail_out);

g.x = ori(1) * [g_core_out g_trail_out]; % leading samples are not used here
g.all.x = ori(1) * [g_lead_out g_core_out g_trail_out];
g.y = ori(2) * [g_core_out g_trail_out]; 
g.all.y = ori(2) * [g_lead_out g_core_out g_trail_out];
g.z = ori(3) * [g_core_out g_trail_out]; 
g.all.z = ori(3) * [g_lead_out g_core_out g_trail_out];

timing.Ncore=g.Ncore;timing.Ntrail=g.Ntrail;timing.Nlead=g.Nlead;

%%
if gen_plots
    nrow=3;ncol=1;
    tt = (0:(g.Nlead+g.Ncore+g.Ntrail-1)) * 1e3*dt_system;

    figure;

    subplot(nrow,ncol,1);
    plot(tt,g.all.x,'r')
    hold;grid on;
    plot(tt,g.all.y,'g')
    plot(tt,g.all.z,'b')
    xlim(tt([1 end]))
    title('gradients')
    legend({'x','y','z'},'location','eastoutside');
    
    subplot(nrow,ncol,2);
    plot(tt,abs( cat(1,zeros([g.Nlead Nc]),rfpulse,zeros([g.Ntrail Nc]))))
    grid
    xlim(tt([1 end]))
    title('RF amplitude')

    subplot(nrow,ncol,3);
    plot(tt,angle( cat(1,zeros([g.Nlead Nc]),rfpulse,zeros([g.Ntrail Nc]))))
    grid
    xlim(tt([1 end]))
    title('RF Phase')

end
%%
end