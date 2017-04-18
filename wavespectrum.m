
% Create wave spectrum using WAFO routines

warning('off','all') % Suppress WAFO warnings

%*************************************************************************%
% Physical frequency cut-off and frequency parameters
wcut2 = max(min(3.0,sqrt(2*g/WS.Hs)),2.2);  % Upper cut
wcut1 = 0.2;                       % Lower cut
dw = 2*pi/Time;                     % delta omega

%*************************************************************************%
% Edit dw and Time after number of provided components
if size(seed,1)>1; dw = (wcut2 - wcut1)/size(seed,1);end
N=ceil((wcut2-wcut1)/dw);          % Non-zero components
omega=(0:dw:wcut2)';               % Only omega below upper cut
M=length(omega);
omega(1)=1e-20;                    % "almost" zero omega for kinematics
domega = omega(3)-omega(2);        % delta omega

%*************************************************************************%
% Check on input
if ~rem(WS.Ndir,2) % Then WS.Ndir is even
    WS.Ndir=WS.Ndir+1;
    disp(['Must have odd number of wave directions, adding one to NDIR=' num2str(WS.Ndir)])
end


%*************************************************************************%
% Wave spreading angles relative to heading angle
theta_tmp = linspace(-pi,pi,WS.Ndir); 


if WS.Hs == 0 || WS.Tp == 0 % Do not create spectrum
    WS.theta  = [];
    WS.dtheta = 1;
    WS.D      = 1;
    WS.Ntheta = 0;
    WS.Spec.S = [];

else
%*************************************************************************%
% Generate Jonswap spectrum for wind sea or total sea
WS.S = JONSWAP(omega,[WS.Hs,WS.Tp,WS.Gamma]);


%*************************************************************************%
% Spreading for wind-sea
if WS.Spread > 1 && WS.Ndir > 1 % Short-crested sea
    
    WS.theta = theta_tmp + WS.Heading; % +/- 180 degrees from Heading
    WS.dtheta = (WS.theta(2)-WS.theta(1)); % Delta theta in radians
    WS.Ntheta = length(WS.theta);
        
    if WS.FreqSpread ~= 1
        % WAFO spreading - frequency independent
        WS.D = spreading(theta_tmp,'cos2s',0,[WS.Spread,WS.Spread],[],0);
    else
        % WAFO spreading - frequency dependent - Tucker (1991) and Krokstad (1998)
        WS.D = spreading(theta_tmp,'cos2s',0,[6.97,9.77,2*pi/WS.Tp,4.06,-1.91,0,1.05,inf],omega,1);
    end

    WS.Spec = mkdspec(WS.S,WS.D,Plott); 
    
else % Long crested
    WS.theta  = WS.Heading;
    WS.dtheta = 1;
    WS.D      = 1;
    WS.Ntheta = 1;
    WS.Spec.S = WS.S.S';
    
end
end


%*************************************************************************%
% Generate Jonswap spectrum for swell
if Swell.Hs > 0 && Swell.Tp > 0
    
    Swell.S = JONSWAP(omega,[Swell.Hs,Swell.Tp,Swell.Gamma]);
    
    % Include swell directions
    
    if Swell.X > 0 && WS.Ndir > 1 % Short-crested swell
        
        Swell.theta = theta_tmp + Swell.Heading;
        
        if Swell.FreqSpread ~= 1
            % Spreading for swell - frequency independent
            Swell.D=spreading(theta_tmp,'poisson',0,[Swell.X,Swell.X],[],0);
        else
            % Spreading for swell - frequency dependent - Bitner-Gregersen (2002)
            Swell.D=spreading(theta_tmp,'poisson',0,[Swell.X,Swell.X,2*pi/Swell.Tp,2.21,-0.35],omega,1);
        end
        % Create swell spectrum
        Swell.Spec = mkdspec(Swell.S,Swell.D,Plott);
        Swell.dtheta = Swell.theta(2)-Swell.theta(1);
        Swell.Ntheta = length(Swell.theta);
    else
        
        Swell.theta = Swell.Heading;
        % Create long-crested swell
        Swell.Spec.S = Swell.S.S';
        Swell.dtheta = 1;
        Swell.Ntheta = 1;
    end
    

else % No swell component
    
    Swell.theta  = [];
    Swell.Ntheta = 0;
    Swell.Spec.S = [];
    Swell.dtheta = 0;
    
end


% All directions
theta = [WS.theta,Swell.theta];
Ntheta = WS.Ntheta+Swell.Ntheta;

% Total spectrum
Spec = [WS.Spec.S'*WS.dtheta,Swell.Spec.S'*Swell.dtheta]; % Omega in dim 1, theta in dim 2


% Total energy in spectrum
Energy = sum(sum(Spec*domega));

