% Simple script for writing spectrum to USFOS formatted file
% IN: phase in radians
%     

function write_spectrum(zeta,omega,head,phase,file)

fid=fopen(file,'w');

if length(head)<length(zeta)
    head=repmat(head,length(zeta),1);
end

fprintf(fid,''' ----------- Wave Components ------------\n');
fprintf(fid,'''\n''   Comp  Amplitude[m]  Period[s]  Dir[deg]  Phase[deg]\n');

for i = 1:length(omega)
    fprintf(fid,'%7d  %7.5e  %7.5e  %7.2f  %7.2f\n'...
        ,i,zeta(i),2*pi/omega(i),head(i),phase(i));
end
fclose(fid);
return