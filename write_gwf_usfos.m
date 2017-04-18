%*************************************************************************%
% Jan-Tore H. Horn
% October 2016
% Write GWF input file to USFOS

function write_gwf_usfos(elev,ux,uy,uz,ax,ay,az,t,x,X,y,Y,z,directory,file)

nz=length(z);
nt=length(t);
ny=length(y);
nx=length(x);

printformat = ' %10.6f %10.6f %10.6f \n';

filename = [directory,file];
fid = fopen(filename,'W');
fprintf(fid,'''\n');
fprintf(fid,'''\n');
fprintf(fid,' GWF Type 132 \n');
fprintf(fid,'''\n');
fprintf(fid,[' GWF   KeyData Size ',num2str(nx),' ',num2str(ny),' ',num2str(nz),' ', num2str(nt),'    1101 \n']); % 1101 = elevation, acc. and velocity
fprintf(fid,'''\n');
fprintf(fid,'''\n');
fprintf(fid,[' GWF   Grid   X   Constant  ' num2str(X(1)) ' ' num2str(X(end)) ' \n']);
fprintf(fid,[' GWF   Grid   Y   Constant '  num2str(Y(1)) ' '  num2str(Y(end)) ' \n']);
   
    
fprintf(fid,' GWF   Grid   Z   Free    \n');

% Print z values
for i=1:nz
    fprintf(fid,' %.3f  ',z(i));
    if rem(i,10) == 0
    fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,'''\n');
fprintf(fid,'''\n');

% Loop for time domain simulation

for i=1:nt

    % Time
    fprintf(fid,[' GWF    Time    ',num2str(t(i)),'\n']);
    fprintf(fid,'''\n');
    

    fprintf(fid,' GWF    Elevation \n');
    fprintf(fid,'''\n');
    fprintf(fid,' %.4f ',elev(i,:));
    fprintf(fid,'\n');
    fprintf(fid,'''\n');
    
    % Velocity   
    fprintf(fid,' GWF    Velocity \n');

    % Print arrays
    fprintf(fid,printformat,[ux(i,:)',uy(i,:)',uz(i,:)']');

    fprintf(fid,'''\n');
    
    % Acceleration
    fprintf(fid,'''\n');
    fprintf(fid,' GWF    Acceleration\n');
    fprintf(fid,'''\n');
    
    % Print arrays
    fprintf(fid,printformat,[ax(i,:)',ay(i,:)',az(i,:)']');

    fprintf(fid,'''\n');
    fprintf(fid,'''\n');
    
    
end
fprintf(fid,''' \n');
fprintf(fid,' GWF    End \n');
fprintf(fid,''' \n');
fprintf(fid,''' \n');
fprintf(fid,''' -------------- e o f --------------');
fclose(fid);

