% Perform simulation
f1=figure;
[X,Y,Z]=meshgrid(x,y,z);
Theta = repmat(theta,1,1,Nx,Ny,Nz);

for i=1:Nt
    
    elev1 = permute(sum(zeta(i,:,:,:),2),[4,3,2,1]);
    elev = repmat(elev1,1,1,Nz);
    id = (Z>elev);
    
    Ux = permute(sum(u(i,:,:,:,:).*cos(Theta),2),[4,3,5,1,2]);
    Uy = permute(sum(u(i,:,:,:,:).*sin(-Theta),2),[4,3,5,1,2]);

    Ax = permute(sum(a(i,:,:,:,:).*cos(Theta),2),[4,3,5,1,2]);
    Ay = permute(sum(a(i,:,:,:,:).*sin(-Theta),2),[4,3,5,1,2]);
    
    Ux(id)=0;
    Uy(id)=0;
    Ax(id)=0;
    Ay(id)=0;    
    
    figure(f1);
    quiver3(X,Y,Z,Ux,Uy,zeros(size(Z)))
    hold on;
    if min(size(elev1))>1; % Plot surface and monopile
        s=surf(X(:,:,1),Y(:,:,1),elev1);
        s.LineStyle='none';
        s.FaceColor='blue';
        s.FaceAlpha=0.3;
    
        if D>0; % Plot monopile
           [x2,y2,z2]=cylinder(D/2);
           xx = repmat(x2(1,:),Nz,1);
           yy = repmat(y2(1,:),Nz,1);
           zz = repmat(z',1,size(xx,2));
           s2 = surf(xx,yy,zz);
           s2.FaceColor = [0.8,0.8,0.8];
           s2.FaceAlpha=0.6;

        end
        fac=1.3;
        axis([fac*min(x),fac*max(x),fac*min(y),fac*max(y),min(z),1.2*max(z)])
        view([30,10])

    else
        view([0,0])
        fac=1.3;
        axis([fac*min(x),fac*max(x),-5,5,min(z),1.2*max(z)])
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Particle velocities')
    hold off;
    
   
    pause(0.05);

end


