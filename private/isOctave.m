
% Function to determine wheter Octave is being used
% Jan-Tore Horn

function flag = isOctave()

    persistent cacheval;
    
    if isempty(cacheval)
        cacheval = (exist ('OCTAVE_VERSION','builtin') > 0 );
    end
    
    flag = cacheval;
end