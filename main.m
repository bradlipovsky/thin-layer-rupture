function out = main(sigma,dc,tides,d,varargin)
tic; 


M.Lup  =  150; %  km
M.Ldown = 100; %  km
M.Vup  =  1e-2; % sliding rate units:  mm/s

M.a    =  0.050;  % Value applied in margins
M.ars  =  0.020;  % Value applied in RS region
M.arw  =  0.020;  % Value applied in RW region

M.b    =  0.005;  % Value applied in margins
M.brs  =  0.015;  % Value applied in RS region
M.brw  =  0.025;  % Value applied in RW region

M.V0   =  M.Vup;  
M.f0   =  0.4;

M.G    =  3.6;    
M.poi  =  1/3;
M.rho  =  0.916;  
M.c    =  sqrt(M.G/M.rho);

% Parameters read as input arguments
W_domain     =  400;  
M.H    =  0.800;
M.tmax =  84600*d;
res = 4;
block=0;
M.W = 120;
M.Wrs = 0;
M.sigma=sigma;
M.Dc = dc;
M.Tides=tides;
M.TidePeriod = 86400;
M.TideHeight = 1e-3;
M.BigLeft=0;
M.prs=1;
M.drs=0;
M.floating=[];
saveit=0;
fig=1;
outfunc=@odeOutFcn;
atol = 1e-9;
rtol = 1e-9;
quiet=0;
M.MaxTimeSteps = 1e20;  
histo=0;
M.TidePhase = 2*pi * 1.28333/24;
bc='zero';
% Process input arguments
if numel(varargin)>0
    for i = 1:numel(varargin)
        if strcmpi(varargin{i},'res')
            res = varargin{i+1};
        end    
        if strcmpi(varargin{i},'W')
            M.W = varargin{i+1};
        end    
        if strcmpi(varargin{i},'Wrs')
            M.Wrs = varargin{i+1};
        end 
        if strcmpi(varargin{i},'block')
            block = varargin{i+1};
        end 
        if strcmpi(varargin{i},'TideHeight')
            M.TideHeight = varargin{i+1};
        end 
        if strcmpi(varargin{i},'TidePeriod')
            M.TidePeriod = varargin{i+1};
        end 
        if strcmpi(varargin{i},'BigLeft')
            M.BigLeft = varargin{i+1};
        end
        if strcmpi(varargin{i},'brs')
            M.brs = varargin{i+1};
        end
        if strcmpi(varargin{i},'brw')
            M.brw = varargin{i+1};
        end      
        if strcmpi(varargin{i},'arw')
            M.arw = varargin{i+1};
        end    
        if strcmpi(varargin{i},'prs')
            M.prs = varargin{i+1};
        end
        if strcmpi(varargin{i},'drs')
            M.drs = varargin{i+1};
        end
        if strcmpi(varargin{i},'a')
            M.a = varargin{i+1};
        end
        if strcmpi(varargin{i},'vup')
            M.Vup = varargin{i+1};
            M.V0 = M.Vup;
        end
        if strcmpi(varargin{i},'floating')
            M.floating = varargin{i+1};
        end
        if strcmpi(varargin{i},'lup')
            M.Lup  =  varargin{i+1};
        end
        if strcmpi(varargin{i},'ldown')
            M.Ldown  =  varargin{i+1};
        end
        if strcmpi(varargin{i},'save')
            saveit  =  varargin{i+1};
        end
        if strcmpi(varargin{i},'fig')
            fig  =  varargin{i+1};
        end
        if strcmpi(varargin{i},'maxtime')
            outfunc = @odeOutFcnTmax;
        end
        if strcmpi(varargin{i},'lowres')
            atol = 1e-6;
            rtol = 1e-5;
        end
        if strcmpi(varargin{i},'hires')
            atol = 1e-10;
            rtol = 1e-10;
        end
        if strcmpi(varargin{i},'quiet');
            quiet=1;
        end
        if strcmpi(varargin{i},'hist');
            histo=varargin{i+1};
        end
        if strcmpi(varargin{i},'L');
            W_domain=varargin{i+1};
        end
        if strcmpi(varargin{i},'bc');
            bc=varargin{i+1};
        end
        if strcmpi(varargin{i},'H');
            M.H=varargin{i+1};
        end
        if strcmpi(varargin{i},'W_domain');
            W_domain=varargin{i+1};
        end
    end
end

% Finite difference stuff
M.Lc = 4*sqrt(M.H*M.Dc*M.G/M.sigma/abs(M.brw-M.arw));
M.Lcrs = 4*sqrt(M.H*M.Dc*M.G/M.prs/M.sigma/abs(M.brs-M.ars));
M.tfric= 1/((M.brw-M.arw)/M.arw * M.V0/M.Dc);
if ~block
    % Minimum 50 nodes in the frictional evolution zone
    L0 = min([M.Lc M.Lcrs]);
    dx = L0/50/res; M.N = numel(0:dx:W_domain);  
    if M.N < 51, M.N = 51; dx = W_domain/(M.N-1); end
    if mod(M.N,2)==0, M.N=M.N+1; end
    [~, M.D2,~] = createSBPOps_mc(M.N,dx,ones(M.N,1),4);
    x = linspace(0,W_domain,M.N);
else
    % Block-slider, no horizontal shear, spring stiffness is G/Lup=G/W
    M.N=1;
    M.Lup = M.W;
end




% Rate weakening region
RwRegion = abs(x-mean(x)) < M.W/2;
M.a = M.a*ones(M.N,1);
M.b = M.b*ones(M.N,1); 
M.b( RwRegion )   = M.brw;
M.a( RwRegion )   = M.arw;

% Rate strengthening region
middle = mean(x);
if M.Wrs > 0
    W0 = (M.W-M.Wrs)/2;
    switch M.BigLeft 
        case 0
            RsRegionOffset = middle;
        case 1
            RsRegionOffset = middle - 0.05*M.W;
            if middle-RsRegionOffset < 4*dx, error('asymmetry is poorly resolved'); end
        case 2
        	RsRegionOffset = middle - 0.02*M.W;
            if middle-RsRegionOffset < 4*dx, error('asymmetry is poorly resolved'); end
    end
    RsRegion = abs(x-RsRegionOffset) < M.Wrs/2;
    M.b( RsRegion ) = M.brs;
else
    W0 = M.W;
    RsRegion = [];
end
if M.N == 1, M.b=M.brw; end % Slider block

% Floating ice
if M.floating
    M.floating = (abs(x-mean(x)) > M.W/2);
end

% RS patch is also a patch of HIGH pressure
if M.prs ~= 1
    prs=M.sigma * M.prs;
    M.sigma = M.sigma * ones(M.N,1);
    M.sigma( RsRegion ) = prs;
    M.sigma = smooth(M.sigma,10,'lowess');
end

% Patch of high Dc
if M.drs ~= 0
    M.Dc = M.Dc * ones(M.N,1);
    M.Dc( RsRegion ) = M.drs;
    M.Dc = smooth(M.Dc,10,'lowess');
end
    
% Initial conditions
if strcmpi(bc,'zero')
    V0 = M.V0*ones(M.N,1);
    Q0 =  M.f0*ones(M.N,1);
    U0=zeros(M.N,1);
    InitialCond = [U0; V0; Q0];
    M.Tau0 = 0;
elseif strcmpi(bc,'fast')
    V0 = 25*M.V0*ones(M.N,1);
    Q0 = M.f0*ones(M.N,1) + 5*(M.b-M.a);
    U0=zeros(M.N,1);
    InitialCond = [U0; V0; Q0];
    M.Tau0 = 0;
end




if quiet
    opt = odeset('abstol',atol,'reltol',rtol);
else
    if ~block
        disp(['dx    = ' num2str(dx)]);
        disp(['Lc/dx = ' num2str(M.Lc/dx)]);
    end
    disp(' ');
    disp(['Lc    = ' num2str(M.Lc)]);
    disp(['Tfric = ' num2str(M.tfric)]);
    disp(['Lcrs  = ' num2str(M.Lcrs)]);
    disp(['N     = ' num2str(M.N)]);
    disp(['W/Lc  = ' num2str(M.W/M.Lc)]);
    disp(['W0/Lc = ' num2str(W0/M.Lc)]);
    disp(['Wrs/Lcrs = ' num2str(M.Wrs/M.Lcrs)]);
    disp(['Lup/Lc = ' num2str(M.Lup/M.Lc)]);
    disp(' ');
    disp('Starting Solver...');
    opt = odeset('abstol',atol,'reltol',rtol,'OutputFcn',outfunc);
end


%
% Run the solver
%
% sol =  ode15s(@ode,[0 M.tmax],InitialCond,opt,M);
sol =  ode15s(@ode,[0 M.tmax],InitialCond,opt,M);
disp(['     Finished ODE15s @ ' num2str(toc)]);

% t   =  sol.x; 
% This saves a ton of space.  Round to the nearest second:
t   =  unique( floor(sol.x) ); 
Y   =  deval(t,sol);
U = Y(1:M.N,:); V = Y(M.N+1:2*M.N,:); Q = Y(2*M.N+1:3*M.N,:);
disp(['     Finished DEval @ ' num2str(toc)]);

% Make the file name
if tides>0
    ti='';
    switch tides 
        case 1
            ti='+t';
        case 2
            ti='+tn';
        case 3
            ti='+tx';
        case 4
            ti='+tz';
        case 5
            ti='+tl';
    end
    if M.TideHeight ~= 1e-3
        ti = [ti '_h' num2str(M.TideHeight)];
    end
    if M.TidePeriod ~= 86400
        ti = [ti '_t' num2str(M.TidePeriod)];
    end
    
else
    ti='';
end
if M.W~=50, pa=['w_' num2str(M.W)];else pa=''; end
if M.Wrs~=0, wr=['+wrs' num2str(M.Wrs)];else wr=''; end
if res~=1, rz=['+r' num2str(res)]; else rz='';end
if block, bl=['+block']; else bl='';end
if M.BigLeft, lft = 'bl'; else lft='';end
if M.BigLeft==2, lft = 'bl2'; end
if M.brs ~= 0.005, brs=['+brs' num2str(M.b(round(M.N/2)))]; else brs=''; end
if M.prs ~= 1, prs = ['+p' num2str(M.prs)]; pprs=M.prs; else prs=''; pprs=1; end
if M.drs ~= 0, drs = ['+drs' num2str(M.drs)]; else drs=''; end
if M.arw ~= 0.05, aa = ['+a' num2str(M.arw(1))]; else aa=''; end
if M.Lup~=150,lup=['+lup' num2str(M.Lup)]; else lup=''; end
if M.Ldown~=100,ldown=['+ldown' num2str(M.Ldown)]; else ldown=''; end
if M.Vup~=0.01,vup=['+vup' num2str(M.Vup)]; else vup=''; end
if rtol<1e-9,hires='+hires'; else hires=''; end
if W_domain~=200,lll=['+L' num2str(W_domain)]; else lll=''; end
name = [lft pa 'n' num2str(M.sigma(1)) 'kpa_dc' ...
	    num2str(M.Dc(1)) 'mm' ti wr rz brs bl prs drs aa lup ldown vup hires lll '.mat'];
thetitle=[num2str(M.sigma(1)) ' kPa ' num2str(M.Dc(1)) ' mm Wrs ' ...
            num2str(M.Wrs) ' brs' num2str(M.brs) ' prs' num2str(M.prs)];

% Set the path for saving files

% path = '/n/regal/rice_lab_seas/lipovsky/';
SavePath = '~/Desktop/simple_out/';
        
if fig
    figure; plot(t,V(round([150 200 250]/400*M.N),:));
    title(thetitle);
    drawnow;
end




% Return some of the event properties
% range = round(numel(t)/4):numel(t);
range = round(3*numel(t)/4):numel(t);
MPD = 1e3;
MPH = max(V(round(100/200*M.N),range)) / 100;
[PKS,LOCS,W] = findpeaks(V(round(M.N/2),range),t(range),'MinPeakDistance',MPD,'MinPeakHeight',MPH);
[PKS2,LOCS2,W2] =findpeaks(V(round(150/400*M.N),range),t(range),...
    'MinPeakDistance',MPD,'MinPeakHeight',MPH);

out.w1=mean(W);
out.v1=mean(PKS);
out.std1=std(PKS);
out.t1=mean(diff(LOCS))/3600;

out.w2=mean(W2);
out.v2=mean(PKS2);
out.std2=std(PKS2);
out.t2=mean(diff(LOCS2))/3600;

out.trec = diff(LOCS)/3600;
out.pks=PKS;
out.locs=LOCS;

if histo
    figure; 
    subplot(121); 
    title(thetitle);
    hist ( diff(LOCS2)/3600);
    subplot(122); 
    hist (PKS2);
end

if saveit    
    save([SavePath name],'t','V','M','Q','U','name','opt','-v7.3');
    disp(['Finished save:  ' SavePath name]);
    toc;
else
    if ~quiet, disp('Not saving result...'); end
end
