function [x,z,residual,xall] = MLS(A,b,z_initial,lambda,maxIT,tol,varargin)
% Solves Magnitude Least Squares problem, using local variable 
% exchange method.
% USAGE:    [x,z] = MLS(A,b,z_initial,maxIT,tol,lambda). 
%           [x,z] = MLS(A,b,z_initial,lambda,maxIT,tol,varargin)
% A is the system 
% matrix. b is a positive real valued target; A and x are both complex
% maxIT and tol are self explanatory. lambda is a tikhonov regularization
% term.
% See Setsompop et al, DOI:10.1002/mrm.21513 for more info

if ~exist('maxIT','var')==1
    maxIT = 50;
end
if ~exist('tol','var')==1
    tol = 1e-4;
end

% look in extra ARGS - 18-6-09
verbose=true;
manstop=true;
saveall=false;
inner_its=60;
for ii=1:numel(varargin)
    if strcmp(varargin{ii},'quiet')
        verbose=false;
    end
    %...
    if strcmp(varargin{ii},'nobox')
        manstop=false;
    end
    %...
    if strcmp(varargin{ii},'saveall')
        saveall=true;
    end
    %...
    if strcmp(varargin{ii},'innerIT')
        inner_its=varargin{ii+1};
    end
end
        
% if no z specified, use random
if numel(z_initial)>1
    z = exp(complex(0,1)*z_initial);
else
   z = exp(1i*2*pi*rand(size(A,1),1));
end
% put in regularizer
if numel(lambda)==1
    if lambda>0
        idx = size(A,1)+(1:size(A,2));
        A = [(A);lambda*(eye(size(A,2)))];
        b = [b;zeros([size(A,2) 1])];
        z = [z;zeros([size(A,2) 1])];
    end
else
    idx = size(A,1)+(1:size(A,2));
    A = [(A);diag(lambda)];
    b = [b;zeros([size(A,2) 1])];
    z = [z;zeros([size(A,2) 1])];
end

% set initial conditions
norm_target = norm(b);
err = [1 1];err_change=1;
it = 0;
b = abs(b);

% 8-12-08, interruptble
if manstop
    hf = figure('name','continue...?');
    set(hf,'position',[420 420 200 100])
    h = uicontrol('Style', 'checkbox', 'String', 'End after current LSQR',...
        'position',[20 20 200 100]);
    drawnow;
end
keepgoing=true;

% allocate cell array to save all solutions
xall={};

while (it<maxIT&&err_change>tol&&keepgoing)
    tic;
    if it>0
        % on all iterations after first, calculate z by setting it equal to phase from Ax
        z_i = A*x;
        z = exp(complex(0,1)*angle(z_i));
        
        % enforce last zeros for regularization are zero
        if lambda>0
            z(idx)=0;
        end            
    end
    
    % now minimise ||Ax - bz||^2 using z above
    if verbose
        x = lsqr((A),b.*z,1e-6,inner_its);
    else
        [x,flag] = lsqr((A),b.*z,1e-6,inner_its);
    end
    % calculate error
    e = norm(abs(A*x)-b) ./ norm_target;
    err(1) = err(2);err(2) = e;
    err_change = err(1)-err(2);
    
    it = it+1;
    
    % save this solution
    if saveall
        xall{it}=x;
    end
    
    residual(it) = e;
    if verbose
        disp(['Iteration ' num2str(it) ' ,  Relative error: ' num2str(e),' change: ' num2str(err_change)]);
    end
    
    it_time=toc;
    if manstop
        pause(it_time/100);
        boxtick=get(h,'Value');
        
        if boxtick==1
            keepgoing=false;
            fprintf(1,'====================\nInterrupted by User\n====================\n')
        end
    end
end

if manstop
    close(hf)
end
end