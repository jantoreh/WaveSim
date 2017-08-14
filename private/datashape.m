
% Data shaping for writing to file efficiently

    if size(theta,2) ~= size(zeta,2)
        theta=theta'; % Theta must be a row vector
    end
    
    
    % Repeat theta to account for whole grid
    Theta = repmat(theta,Nt,1,Nx,Ny,Nz);
    
    ncol = Nz*Ny*Nx; % Total number of columns in files

    % Create zero-dummy
    dummy = zeros(Nt,ncol);
    
    %*********************************************************************%
    % Correct for wave elevation
    elev = sum(zeta(:,:,:,:),2); % Sum over all directions in all points
    Z = repmat(permute(z,[1,3,4,5,2]),Nt,1,Nx,Ny,1);
    id = (Z>repmat(elev,1,1,1,1,Nz)); % Find points larger than elevation
    id = repmat(id,1,Ntheta,1,1,1);
    u(id)=0;
    a(id)=0;
    
    
    %*********************************************************************%
    % Re-arrange data
    
    % Re-arrange zeta
    
    elev=permute(elev,[1,4,3,2]); % Permute to facilitate reshape, now z is in dim 2, y in dim 3, x in dim 4
    elev = reshape(elev,Nt,Nx*Ny);% Reshape so that (X1,Y1..Yn),(X2,Y1..Yn)
    
    
    % Velocities
    ux = sum(u.*cos(Theta),2);
    ux = permute(ux,[1,5,4,3,2]);
    ux = reshape(ux,Nt,ncol);
    
    uy = sum(u.*sin(-Theta),2);
    uy = permute(uy,[1,5,4,3,2]);
    uy = reshape(uy,Nt,ncol);
    
    uz = dummy;
    
    
    % Accelerations
    ax = sum(a.*cos(Theta),2);
    ax = permute(ax,[1,5,4,3,2]);
    ax = reshape(ax,Nt,ncol);
    
    ay = sum(a.*sin(-Theta),2);
    ay = permute(ay,[1,5,4,3,2]);
    ay = reshape(ay,Nt,ncol);
    
    az = dummy;
    
    
    % EOF %