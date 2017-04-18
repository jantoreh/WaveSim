function lambda=wave_length(wfreq,h,g)
% function lambda=wave_length(omega,h,g)
% --- Description: ---------------------------------
%
% --- Parameters: ----------------------------------
% I/O	Name:		Size/Description:
% I	x 
%  
% --- Method: --------------------------------------
%
% --- Uses: ----------------------------------------
%
% --- Programed by: --------------------------------

%finner b?lgelengden ved hjelp av newtons metode
%wfreq=b?lgefrekvensen i rad/sec.
%h=vanndybden.
%
% Torgrim Driveklepp, sommer 2000
% modifisert JRK 10 Nov 2005. Rutinen er skrevet om til ? ta vektorer 
%                             av b?lgelengder. 
% Mod JT. Horn Jun 2015 til baade vektorer og skalarer.
%
if length(wfreq)>1
wfreq(1) = 0.000001;
end
T = 2*pi./wfreq;
lambda=(g/(2*pi)*T.^2);
feilmargin=0.0001;
feil=feilmargin+100;
while feil>=feilmargin
lambda_ny=1/2*lambda*g.*T.^2.*(-2*pi*h+2*pi*h.*tanh(2*pi*h./lambda).^2 ...
         -lambda.*tanh(2*pi*h./lambda))./...
         (pi*(-g*T.^2*h+g*T.^2*h.*tanh(2*pi*h./lambda).^2-lambda.^2));
   feil=abs(lambda_ny-lambda)./lambda;
   feil=max(feil);
   lambda=lambda_ny;
end
