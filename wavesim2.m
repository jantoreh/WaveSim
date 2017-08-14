%
% WaveSim - Matlab format
% Jan-Tore H. Horn, Feb 2016
%
% Version: 2.00
%
% Update history:
%
% 25.10.16 JT Horn - Skrevet om til 2D spectrum og 3D grid
% 28.10.16 JT Horn - FAST format utskrift i 3D grid (+vpOne selvsagt)
% 29.10.16 JT Horn - Simulering av partikkel-hastigheter
% 
%
%
%*************************************************************************%
function OUT=wavesim2(varargin)

%*************************************************************************%
% Check input

dir = './'; % Initialize the present location as the output directory

if rem(nargin,2)
    error('Number of input arguments must be ZERO or EVEN')
end



%*************************************************************************%
% Global variables
g           = 9.81;
format long                 % Printed precision
file_gwf    = 'gwf.w132';
file_spec   = 'spec.txt';

%*************************************************************************%
% Read input-file

file='wavesim2.inp';
fid=fopen(file,'r');
if fid<0; error('Cannot find input-file at current location..');end
fgetl(fid);                                          % Header
RunID  = fscanf(fid,'%f',1);fgetl(fid);              % Wave option
FORMAT = fscanf(fid,'%f',1);fgetl(fid);              % Output format, vpOne/USFOS=1, FASTv8=2
Sim    = fscanf(fid,'%f',1);fgetl(fid);             % Simulation option
Plott = Sim;
Time   = fscanf(fid,'%f',1);fgetl(fid);              % Simulation time
dt     = fscanf(fid,'%f',1);fgetl(fid);
Time = Time+dt; % One extra time-step
seed   = fscanf(fid,'%f',1);
if exist('seed.txt','file'); seed=load('seed.txt','-ascii');end % Import seed.txt if it exists in this folder
if isempty(seed); seed_file=fscanf(fid,'%s',1);     % Then path to seed file is given: Uniform(0,2pi)
seed=load(seed_file);end                            % Read seed file
fgetl(fid);                                         % Continue
h      = fscanf(fid,'%f',1);fgetl(fid);
Ht     = fscanf(fid,'%f',1);fgetl(fid);fgetl(fid);
D      = fscanf(fid,'%f',1);fgetl(fid);fgetl(fid); % MacCamy Fuchs diameter
if D>0; MacCam=1; else MacCam=0; end

WS    = struct(); % Struct for wind sea/total sea
Swell = struct(); % Struct for swell
WS.Hs       = fscanf(fid,'%f',1);fgetl(fid);
WS.Tp       = fscanf(fid,'%f',1);fgetl(fid);
WS.Gamma    = fscanf(fid,'%f',1);fgetl(fid);
WS.Heading  = eval(fscanf(fid,'%s',1))*pi/180;fgetl(fid);
WS.Spread  = fscanf(fid,'%f',1);fgetl(fid);
WS.Ndir     = fscanf(fid,'%f',1);fgetl(fid);fgetl(fid);
Swell.Hs = fscanf(fid,'%f',1);fgetl(fid);
Swell.Tp = fscanf(fid,'%f',1);fgetl(fid);
Swell.Gamma = fscanf(fid,'%f',1);fgetl(fid);
Swell.Heading = eval(fscanf(fid,'%s',1))*pi/180;fgetl(fid);
Swell.X    = fscanf(fid,'%f',1);fgetl(fid);fgetl(fid);

% Current
Current.Heading = eval(fscanf(fid,'%s',1))*pi/180;fgetl(fid);
Current.Velocity = fscanf(fid,'%f',1);fgetl(fid);
Current.Power = fscanf(fid,'%f',1);fgetl(fid);fgetl(fid);


ORDER  = fscanf(fid,'%f',1);fgetl(fid);
if ORDER == 1.2; wheeler=1; end % Wheeler stretching activated
WAMIT2s = fscanf(fid,'%d',1);fgetl(fid);
WAMIT2s_file = fscanf(fid,'%s',1);fgetl(fid);
WAMIT2s_node = fscanf(fid,'%d',1);fgetl(fid);fgetl(fid);
randAmp = fscanf(fid,'%d',1);fgetl(fid);
if randAmp ~= 1; randAmp=0; end


% X-distribution
Nx = fscanf(fid,'%d',1);
x1 = fscanf(fid,'%d',1);
x2 = fscanf(fid,'%d',1);fgetl(fid);
if Nx==1
    x = mean([x1,x2]);
else
    x  = linspace(x1,x2,Nx);
end

% Y-distribution
Ny = fscanf(fid,'%d',1);
y1 = fscanf(fid,'%d',1);
y2 = fscanf(fid,'%d',1);fgetl(fid);
if Ny==1
    y=mean([y1,y2]);
else
    y  = linspace(y1,y2,Ny);
end

% Z-distribution
Z_OPTION  = fscanf(fid,'%d',1);fgetl(fid);
Z_MAX   = fscanf(fid,'%d',1);fgetl(fid);
Nz     = fscanf(fid,'%f',1);fgetl(fid);
if Z_OPTION==3 && Nz > 0
    for i=1:Nz
        z(i) = fscanf(fid,'%f',1);fgetl(fid);
    end
    z=sort(z);
end

% Adjust Z-distribution for tidal elevation
z=z+Ht; % Caluculate for larger/smaller z-values
z(z<-h)=[]; % Remove evaluation points below surface
Nz = length(z);

% Close inputfile
fclose(fid);



%*************************************************************************%
% Default values only changed with varargin
WS.FreqSpread=0; % Frequency dependent spreading for wind sea
Swell.FreqSpread=0; % Frequency dependent spreading for swell
Spectrum = 0; % No spectrum is given

%*************************************************************************%
% Adjust variables based on input


for i = 1:2:nargin
    
    name = varargin{i};
    value = varargin{i+1};
    
    switch name
        
        case 'directory' % Output directory is given
            
            dir = [value,'/']; % Output directory in first entry
            if ~exist('dir','dir')  % Create output directory if it does not exist
               mkdir(dir);
            end
            
        case 'Switch1'
            Switch1 = value;
        case 'Hs'
            WS.Hs = value;
        case 'Tp'
            WS.Tp = value;
        case 'Heading'
            WS.Heading = value*pi/180;
        case 'Hs_swell'
            Swell.Hs = value;
        case 'Tp_swell'
            Swell.Tp = value;
        case 'Heading_swell'
            Swell.Heading = value*pi/180;
        case 'NSPREAD'
            WS.Spread = value;
        case 'NDIR'
            WS.Ndir = value;
        case 'X'
            Swell.X = value;
        case 'D'
            D = value;
        case 'h'
            h = value;  
        case 'Time'
            Time = value;
        case 'RunID'
            RunID = value;
        case 'seed'
            seed = value;
        case 'FreqSpread' % Frequency-dependent wave spreading
            WS.FreqSpread=value;
            Swell.FreqSpread=value;
        case 'Spectrum' % Full spectrum is given as struct
            Spectrum = value;
    end
    
    
    
end


%*************************************************************************%
% Basic check on inputs compared to program limitations

if RunID < 1 || RunID > 3 % Not implemented option
    disp('Nothing to do...')
    return
end


if ORDER > 1.5 && WS.Spread > 0
    error('This input is expected to be too time-consuming.. Consider switching off second order waves or spreading.')
end

if ORDER > 1.5 && Ny > 1 
    error('Second order kinematics are not implemented for multi-directional waves, only one y-value.')
end



if ~isstruct(Spectrum) % Spectrum is not given and must be calculated
    
%*************************************************************************%
% Create wave spectrum based on current Hs and Hs_swell
wavespectrum;


   
else % Spectrum is given as input
   omega = Spectrum.omega;
   Spec = Spectrum.S;
   dw = Spectrum.dw;
   theta = Spectrum.theta;
   Ntheta = length(theta);
   wcut1 = Spectrum.wcut1;
   M = Spectrum.M;
   N = Spectrum.N;
   
   
end

    % Assign output
    OUT.S      = Spec;
    OUT.omega  = omega;
    OUT.theta  = theta;
    OUT.wcut1  = wcut1;
    OUT.M      = M;
    OUT.N      = N;
    OUT.dw = dw;
    OUT.Tp = WS.Tp;
    
if RunID == 1 % Escape from WaveSim 
    return
end

%*************************************************************************%
% Create z-distribution based on extreme value statistics:
if ~exist('z','var')
    plott=0;
    z = z_distribution(Z_MAX,Nz,h,Z_OPTION,plott);                 %0=log, 1=uniform
end

%*************************************************************************%
% Modify wave amplitudes
zeta_a1=sqrt(2*Spec*dw);                                   % Amplitudes
id1=find(omega>wcut1,1);                                % Lower cut
zeta_a1(1:id1-1,:) = 0;                                 % Lower cut
k = 2*pi./wave_length(omega',h,g)';    % Wave numbers

%*************************************************************************%
% Generate random seeds
phase=zeros(M,Ntheta);
if size(seed,1)>1 && size(seed,2) == 1 && Ntheta==1
    phase(id1:id1+N-1)=seed(:,1);       % Use provided phases for non-zero components
    rand('seed',1)                      % Set random amplitude seed
    A=rand(M,1);                        % Draw for random amplitudes
elseif size(seed,2) == 2 && Ntheta==1   % Two random variables per wave component
    phase(id1:id1+N-1,:) = -atan2(seed(:,2),seed(:,1)); % Phase
    A = sqrt(seed(:,1).^2+seed(:,2).^2)/sqrt(2);        % Amplitude correction
else
    rand('seed',seed)                   % Seed for angular variation
    phase=rand(M,Ntheta)*2*pi;          % Random phase
    rand('seed',round(seed*1000))       % Seed for amplitude variation
    A=rand(M,Ntheta);                   % Draw for random amplitudes
end

%*************************************************************************%
% Randomize amplitudes

if randAmp==1
    zeta_a = wblinv(A,sqrt(zeta_a1.^2/2)*sqrt(2),2);
    %zeta_a = raylinv(A,sqrt(zeta_a1.^2/2));    % Random amplitudes
    zeta_a(~isfinite(zeta_a))=0;               % Remove roundoff error
elseif size(seed,2) == 2 && Ntheta==1         % Amplitudes from 2D random seed
    zeta_a = zeta_a1;
    zeta_a(id1:id1+N-1,:) = A.*zeta_a(id1:id1+N-1,:);
else
    zeta_a=zeta_a1;
end


%*************************************************************************%
% Time-vector for total simulation
t=0:dt:(Time-dt);    % Simulation + startup
Nt=length(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunID == 2                              % Only proceed for GWF

%*************************************************************************%
% Time-vector for FFT realization
Nt_fft=ceil(2*pi/(dw*dt))-1;
t_fft=0:dt:Nt_fft*dt-dt;


if length(t)>length(t_fft) % FFT signal has to be repeated for Tstart>0
    repeat=1;
else
    repeat=0;
end

%*************************************************************************%
% Adjust MacCamy & Fuchs phase for first order kinematics
if MacCam==1
       r=D/2;
       phase_inertia=phase-repmat((-5.5*k.^3*r^3+11100*k.*k.*r*r-1817*k*r+119)...
       ./(k.^3*r^3+6164*k.^2*r^2-1514*k*r+6276),1,size(phase,2));
else
    phase_inertia=phase;
end

%*************************************************************************%
% Perform FFT for first order kinematics
zeta1=zeros(Nt_fft,Ntheta,length(x),length(y));
u=zeros(Nt_fft,Ntheta,length(x),length(y),Nz);
a=zeros(Nt_fft,Ntheta,length(x),length(y),Nz);


for i=1:Ntheta % Loop over all directions
      [zeta1(:,i,:,:),u(:,i,:,:,:),a(:,i,:,:,:),~]=first_order_fft(zeta_a(:,i),...
          omega,k,phase(:,i),phase_inertia(:,i),Nt_fft,g,x,y,z,h,M,ORDER,...
          MacCam,D/2,0,theta(i)); % Transfer from NB to XY coordinate system
end


%*************************************************************************%
% Second order 2D FFT for incident wave


if ORDER > 1.5 % Second order is wanted  
    zeta_sum=zeros(Nt_fft,Ntheta,Nx,1);                               % Pre-allocate
    a_sum=zeros(Nt_fft,Ntheta,Nx,1,Nz);
    u_sum=a_sum;
    for j1 = 1:Ntheta % For all directions
      for jx = 1:Nx % For all x-values
        for jz = 1:Nz % For all z-values
            if jz==1 % Only find elevation for first call
                [zeta_sum(:,j1,jx,1),u_sum(:,j1,jx,1,jz),a_sum(:,j1,jx,1,jz)]=second_order_fft(zeta_a(:,j1),...
            omega,k,phase(:,j1),h,x(jx),z(jz),Nt_fft,1);
            else % Find velocities and accelerations for remaining calls
                [~,u_sum(:,j1,jx,1,jz),a_sum(:,j1,jx,1,jz)]=second_order_fft(zeta_a(:,j1),...
            omega,k,phase(:,j1),h,x(jx),z(jz),Nt_fft,0);
            end
        end
      end
    end

    if ORDER > 2                                   % Assign zeta
        zeta=zeta1+zeta_sum;
    else
        zeta=zeta1;
    end
    
    if ORDER < 3                  % Cut off second order potential at z=0;
        id=find(z>0,1);
        u_sum(:,:,:,:,id:end)=0;
        a_sum(:,:,:,:,id:end)=0;
    end
    
      % Sum first and second order potential contributions
      u=u+u_sum;
      a=a+a_sum;
else
    zeta=zeta1;
end

%*************************************************************************%
% Second order Wamit loads

if (WAMIT2s == 1) && (ORDER <= 1) && (ORDER > 0)
    if Ntheta>1; error('WAMIT forces not implemented for multiple directions');end
    [Fwam2_1]=SOWL(WAMIT2s_file,omega,phase,zeta_a,floor(Time/dt),1);
    [Fwam2_5]=SOWL(WAMIT2s_file,omega,phase,zeta_a,floor(Time/dt),5);
    file_wam='wamit2.fem'; % Output in current folder
    if repeat==1
        n=Nt-Nt_fft;
        Fwam2_1(Nt_fft+1:Nt,:)=Fwam2_1(1:n,:);
        Fwam2_5(Nt_fft+1:Nt,:)=Fwam2_5(1:n,:);
    end
    write_timehist(4,4,WAMIT2s_node,t,Fwam2_1+Fwam2_5/h,'',file_wam,1,'w');
end

%*************************************************************************%
% Repeat FFT-signal if necessary
if repeat==1
    n=Nt-Nt_fft;
    zeta(Nt_fft+1:Nt,:,:,:)=zeta(1:n,:,:,:);
    u(Nt_fft+1:Nt,:,:,:,:)=u(1:n,:,:,:,:);
    a(Nt_fft+1:Nt,:,:,:,:)=a(1:n,:,:,:,:);
end

% No elevation for ORDER 1 or smaller
if ORDER <= 1
    zeta = zeros(size(zeta));
end


%*************************************************************************%
% Current velocity
% Stretching of current profile to free surface
z_tide = (h+sum(zeta,2)).*(1+permute(z,[1,3,4,5,2])/h)-h;
u_tide = Current.Velocity*((h+z_tide)/h).^Current.Power;
u_tidex = u_tide.*cos(Current.Heading);
u_tidey = -u_tide.*sin(Current.Heading);

dummy = zeros(size(cat(2,u_tidex,u_tidey)));
u=cat(2,u,u_tidex,u_tidey); % Add current velocities
a = cat(2,a,dummy);
zeta=cat(2,zeta,dummy(:,:,:,:,1));
theta=[theta,Current.Heading,Current.Heading]; % Add current directions
Ntheta=length(theta);



%*************************************************************************%
% Plot to gridwave file

% Shape data to x,y,z directions
datashape;



if FORMAT == 1 % USFOS/vpOne grid: .w33 file                                                  % Here
    write_gwf_usfos(elev,ux,uy,uz,ax,ay,az,t,x,[x1,x2],y,[y1,y2],z,dir,file_gwf);
elseif FORMAT == 2 % FAST Grid: .Vxi,.Vyi,.Vzi,.Axi,.Ayi,.Azi,.Elev
    % Specify filename without extension
    file = 'gwf';
    write_gwf_fast(elev,ux,uy,uz,ax,ay,az,dir,file);
else
    error('Grid format not compatible')
end
    
    
else  % Only print spectrum and wave elevation in (0,0)

    
%*************************************************************************%
% Write spectrum to USFOS-format (only non-zero amplitudes)
% Modify phase to coincide with GWF calculations: phase_mod= -(phase+pi/2).
% verified with usfos version 874

Ncomp = numel(zeta_a(id1:end,:));
n1 = size(zeta_a(id1:end,:),1);
n2 = size(zeta_a(id1:end,:),2);

write_spectrum(reshape(zeta_a(id1:end,:),Ncomp,1),repmat(omega(id1:end),n2,1),...
    180-reshape(repmat(theta*180/pi,n1,1),Ncomp,1),reshape(-phase(id1:end,:)*180/pi-90,Ncomp,1),file_spec); % Convert theta to degrees

% Replace NFREQ in head.fem (if it exists) with number of elements in S
if exist('./head.fem','file')
    substitute('head.fem','NFREQ',Ncomp)
end



% Print wave elevation at (0,0)
    Nt_fft=ceil(2*pi/(dw*dt));
    t_fft=(0:dt:(Nt_fft-1)*dt)';
    zeta1=zeros(Nt_fft,1);
  for j=1:Ntheta
    [zeta1(1:Nt_fft,j),~,~,~]=first_order_fft(zeta_a(:,j),omega,k,phase(:,j),phase(:,j),...
     Nt_fft,g,0,0,0,h,M,ORDER,MacCam,D/2,0,theta(j));
  end 

    % Write time series to file STIME_HS_TP_SEED.Elev
    file_elev=[num2str(Time,'%.0f') '_' num2str(Hs,'%.1f'),'_',num2str(Tp,'%.1f'),'_',num2str(seed,'%d'),'.Elev'];
    tmp=[t_fft,sum(zeta1,2)]; % Sum over theta
    save(file_elev,'tmp','-ascii');

end


%*************************************************************************%
% Simulation
if Sim == 1 && RunID>1 && min([Ny,Nx]) > 1 && ~isOctave% Perform live wave simulation in figure
    wavesimulation;
elseif isOctave && Sim == 1 && RunID>1
    disp('Simulation is not performed in Octave')
elseif min([Ny,Nx])==1
    disp('Too few gridpoints in X or Y-direction for simulation')
end

return