
% Jan-Tore Horn, Oct 2016

% Write gid kinematics to FAST format

function write_gwf_fast(elev,ux,uy,uz,ax,ay,az,dir,file)


    
    % Filenames, 8 in total
    fvx = [dir,file,'.Vxi'];
    fvy = [dir,file,'.Vyi'];
    fvz = [dir,file,'.Vzi'];
    fax = [dir,file,'.Axi'];
    fay = [dir,file,'.Ayi'];
    faz = [dir,file,'.Azi'];
    fp  = [dir,file,'.DynP'];
    fel = [dir,file,'.Elev'];
    dummy = zeros(size(ux));
    
    
    %*********************************************************************%  
    % Save data
    
    % Save elevation to file
    save(fel,'elev','-ascii');
    
    % Save dummy to dynamic pressure file
    save(fp,'dummy','-ascii');
    
    % Save velocities
    save(fvx,'ux','-ascii');
    save(fvy,'uy','-ascii');
    save(fvz,'uz','-ascii');
    
    % Save accelerations
    save(fax,'ax','-ascii');
    save(fay,'ay','-ascii');
    save(faz,'az','-ascii');
    
    

return