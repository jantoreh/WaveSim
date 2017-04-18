%*************************************************************************%
% Jan-Tore H. Horn
% August 2015


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Time-Domain-simulation in irregular seas  %
%     Calculating second order kinematics     %
% ------------------------------------------- %
%  Calculates second order sum frequency term %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT:
% zeta_a    Wave amplitudes
% omega     Wave component frequencies with constant delta_omega. First
%           value must be zero.
% k         Wave component wave number
% phase     Wave component phases
% h         Water depth
% x         x-component, scalar value
% z         z-value (scalar) for kinematics calculation, taken as 0 if z>0.
% n_fft     Number of FFT components, to obtain the correct length of the time series
% ELEV      Logical. 1=calculate wave elevation, else only kinematics

% OUTPUT:
% F_plus    Sum frequency wave elevation in vector.
% u_plus    Sum frequency velocity components in vector.
% a_plus    Sum frequency acceleration components in vector.

function [F_plus,u_plus,a_plus]=second_order_fft(zeta_a,omega,k,phase,h,x,z,n_fft,ELEV)

% Define variables
g = 9.81;
M = length(omega);

nx=length(x);
nz=length(z);

% Basic check
if nx>1||nz>1
    error('Too many x or z values')
end
if size(phase,2)>1
    phase = phase';
end
if size(omega,2)>1
    omega=omega';
end

%Setting up important values to be used in surface elevation and kinematic calculations, the D's
% For FFT:
R = k.*tanh(k.*h);
ki=k;
kj=k';
Ri=R;
Rj=R';
kij=ki*kj;
Rij=Ri*Rj;
omegai=omega;
omegaj=omega';
omegaij=omegai*omegaj;
zeta_ai=zeta_a;
zeta_aj=zeta_a';
zeta_aij=zeta_ai*zeta_aj;
phasei=phase - ki*x;
phasej=phase' - kj*x;

% k and omega sum
k_plus=(repmat(ki,1,M)+repmat(kj,M,1));
w_plus=repmat(omegai,1,M)+repmat(omegaj,M,1);

% Denominator
denom_plus=(repmat(sqrt(Ri),1,M)+repmat(sqrt(Rj),M,1)).^2-k_plus.*tanh(k_plus*h);

% D for sum frequency
D_plus=((repmat(sqrt(Ri),1,M)+repmat(sqrt(Rj),M,1)).*((ki.*ki-Ri.*Ri)*sqrt(Rj)+sqrt(Ri)*(kj.*kj-Rj.*Rj))+...
    2*(repmat(sqrt(Ri),1,M)+repmat(sqrt(Rj),M,1)).^2.*(ki*kj-Ri*Rj))./denom_plus;

%clear denom_plus

% Removing round-off errors in matlab creating inf and nan numbers
D_plus(~isfinite(D_plus))=0;

% The QTF for wave elevation
QTF_plus=0.25*(D_plus-(kij-Rij))./sqrt(Rij)+0.25*repmat(Rj,M,1,1)+0.25*repmat(Ri,1,M,1);

% Removing round-off errors in matlab creating inf and nan numbers
QTF_plus(~isfinite(QTF_plus))=0;

if ELEV==1
    
% 2-D Fourier coefficients in 3D:
f_plus=zeta_aij.*QTF_plus.*(exp(-1i*phasei)*exp(-1i*phasej));


% FFT over row and column to abtain time series
F_plus=diag(real(fft2(f_plus,n_fft,n_fft)));

else
   F_plus=0; 
end
clear f_plus

% Corrected z, max(z)=0
if z>0
    z_corr=0;
else
    z_corr=z;
end

%Calculate amplitudes for second order kinematics
Z_plus=0.25*g*g*k_plus.*zeta_aij./omegaij./cosh(k_plus*h).*D_plus./w_plus.*cosh(k_plus.*(z_corr+h));
clear D_plus
U_plus=Z_plus.*(exp(-1i*phasei)*exp(-1i*phasej));
A_plus=-U_plus.*w_plus;

% Padding with zeros
U_plus(n_fft,n_fft)=0;
A_plus(n_fft,n_fft)=0;


% FFT for sum frequency kinematics for the complete grid
if isOctave % Seems to be fastest in Octave
    a_plus=diag(-imag(fft2(A_plus,n_fft,n_fft)));
    u_plus=diag(real(fft2(U_plus,n_fft,n_fft)));
else % Seems to be fastest in Matlab
    a_plus=diag(-imag(fft(fft(A_plus,n_fft).',n_fft).'));
    u_plus=diag(real(fft(fft(U_plus,n_fft).',n_fft).'));
end




end