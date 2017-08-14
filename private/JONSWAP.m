%*************************************************************************%
% Jan-Tore H. Horn
% August 2015
%*************************************************************************%


function Spec=JONSWAP(omega,dat)


if size(omega,2)>1
    omega=omega';
end

Hs=dat(1);
Tp=dat(2);
gamma=dat(3);

%Calculate gamma
if gamma == 0
if Tp/sqrt(Hs) <= 3.6
    gamma = 5;
elseif Tp/sqrt(Hs) >= 5
    gamma = 1;
else
    gamma = exp(5.75 - 1.15*Tp/sqrt(Hs));
end
end

wp=(2*pi)/Tp;

id=find(omega>wp,1);
sigma=ones(length(omega),1)*0.07;
sigma(id:end,1)=0.09;


W = wp./omega;

S = ((5/(32*pi))*Hs^2*Tp)*W.^5.*exp(-1.25*W.^4).*((1 - ...
        0.287.*log(gamma))).*gamma.^(exp(-1./(2*sigma.^2).*(1./W-1).^2));

S(isnan(S))=0;

% Assign output to struct
Spec=struct();
Spec.S=S;
Spec.w=omega;
Spec.h=Inf; % Infinite water depth assumed
Spec.note=['JONSWAP, Hm0 = ',num2str(Hs),', Tp = ',num2str(Tp),', gamma = ',num2str(gamma)];
Spec.date=datestr(now);
Spec.type='freq';
Spec.phi=0;
Spec.norm=0;
Spec.tr=[];




return
