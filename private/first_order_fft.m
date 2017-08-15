%*************************************************************************%
% Jan-Tore H. Horn
% August 2015
% Perform first order surface simulation with kinematics using FFT
%*************************************************************************%

function [zeta,u,a,FNV3]=first_order_fft(zeta_a,omega,k,phase,phase2,n_fft,g,x,y,z,h,M,ORDER,MACCAM,r,opt_fnv3,theta)

%*************************************************************************%
% Pre-process
%*************************************************************************%



% x and z in rows
if size(x,1) ~= 1;x = x';end
if size(y,1) ~= 1;y = y';end
if size(z,1) ~= 1;z = z';end
if size(omega,2) ~= 1;omega = omega';end

nx=length(x);
ny=length(y);
nz=length(z);


%*************************************************************************%
% Fourier coefficients, x in columns (dim2)
%*************************************************************************%

phi=repmat(phase,[1,nx,ny])-repmat(k*x*cos(theta),[1,1,ny])+repmat(permute(k*y*sin(theta),[1,3,2]),[1,nx,1]);
phi_MACCAM=repmat(phase2,[1,nx,ny])-repmat(k*x*cos(theta),[1,1,ny])+repmat(permute(k*y*sin(theta),[1,3,2]),[1,nx,1]); % y in dim 3
C=repmat(zeta_a,[1,nx,ny]).*exp(-1i*phi);
C_MACCAM=repmat(zeta_a,[1,nx,ny]).*exp(-1i*phi_MACCAM);

%*************************************************************************%
% Elevation, x values in columns, y-values in dim3
%*************************************************************************%
zeta=real(fft(C,n_fft));

%*************************************************************************%
% Correction of z, constant extrapolation above mean waterline (z=0), z in
% dim 4
%*************************************************************************%

    z_corr = z;
    z_corr(z_corr>0) = 0;
    index = find(z>0,1);

    Decay = repmat(permute(cosh(k*(h+z_corr))./repmat(cosh(k*h),1,nz),[1 4 3 2]),[1,nx,ny]); % omega in dim 1, x in dim 2, y in dim 3, z in dim 4
    
    
    % Velocity
    Cu_tmp=C*g.*repmat(k./omega,[1,nx,ny]); % 3 dims
    Cu=repmat(Cu_tmp,[1,1,1,nz]).*Decay; % z in dim4
    
    
    % Acceleration
    Ca_tmp=C_MACCAM*g.*repmat(k./omega,[1,nx,ny]);
    Ca_tmp=repmat(Ca_tmp,[1,1,1,nz]).*Decay;
    Ca_tmp=Ca_tmp.*repmat(omega,[1,nx,ny,nz]);
    Ca=Ca_tmp.*exp(-1i*pi/2); % Phase of acceleration is pi/2 before elevation and velocity !!!
    
    % Stretching for higher order kinematics
    if ORDER > 2
        Decay_stansberg = repmat(sinh(k*h)./cosh(k*h),[1,nx,ny,nz]);
        Decay_stansberg(:,:,:,1:index-1) = 0;
        
        Cu_dz=repmat(Cu_tmp.*repmat(k,[1,nx,ny]),[1,1,1,nz]).*Decay_stansberg;
        Ca_dz=Cu_dz.*repmat(omega,1,[nx,ny,nz]).*Decay_stansberg;
        Cu=Cu+Cu_dz.*repmat(reshape(z,1,1,1,nz),[M,nx,ny]);
        Ca=Ca+Ca_dz.*repmat(reshape(z,1,1,1,nz),[M,nx,ny]);
    end


    if MACCAM==1 % MacCamy and Fuchs correction to acceleration
        Ca=Ca.*repmat(1./(2).*(0.581.*k.*k.*r^2+0.718*k*r+0.78)...
           ./(k.^3.*r^3-0.256*k.^2*r^2+0.381*k*r+0.389),[1,size(Ca,2),size(Ca,3),size(Ca,4)]);
    end
    
    % Perform FFT
    u=real(fft(Cu,n_fft));    
    a=real(fft(Ca,n_fft));
    
utz_complete=imag(fft(Cu.*repmat(omega.*k,[1,nx,ny,nz]),n_fft));


%*************************************************************************%
% FNV3 if applicable
%*************************************************************************%
if opt_fnv3>0
    u=real(fft(repmat(k*g.*zeta_a./omega,[1,nx,ny]).*exp(-1i*phi),n_fft,1)); % Verified

    ut=imag(fft(repmat(k*g.*zeta_a,[1,nx,ny]).*exp(-1i*phi),n_fft,1));       % Verified
    
    ux=imag(fft(repmat(-k.^2./omega.*g.*zeta_a,[1,nx,ny]).*exp(-1i*phi),n_fft,1)); % Verified
    
    utz=imag(fft(repmat(k.*omega.^2.*zeta_a,[1,nx,ny]).*exp(-1i*phi),n_fft,1)); % Verified
    
    w=imag(fft(repmat(-g*k.*zeta_a./omega,[1,nx,ny]).*exp(-1i*phi),n_fft,1));  % Verified
    
    wx=real(fft(repmat(-omega.*k.*zeta_a,[1,nx,ny]).*exp(-1i*phi),n_fft,1));  % Verified
    
    wt=real(fft(repmat(g.*k.*zeta_a,[1,nx,ny]).*exp(-1i*phi),n_fft,1));   % Verified
    

    FNV3=struct('u',u(:,1,end),'ut',ut,'utz',utz,'w',w,'wx',wx,'wt',wt,'ux',ux,'utz_complete',utz_complete);
else
    FNV3=0;
end

% Now permute x,y,z from dim 2 3 4 to dim 3 4 5;
zeta=permute(zeta,[1,5,2,3,4]);
a = permute(a,[1,5,2,3,4]);
u = permute(u,[1,5,2,3,4]);




end