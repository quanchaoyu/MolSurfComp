function read_inputs
% c
% c   Read inputs
% c

global Para; % parameters

disp('-------------------');
disp('! Input Parameters');
disp('-------------------');

fid = fopen('Inputs.txt','r');

line = fgetl(fid);

while ischar(line)
    [token,remain] = strtok(line);
    
    flag = 1;
    switch(token)
        % Filename & Option
        case 'filename', Para.filename = strtok(remain);
        case 'format', Para.format = strtok(remain);
        case 'option', Para.option = strtok(remain);
        
        % Parameters    
        case 'radius_probe', Para.radius_probe = sscanf(remain,'%f');
        case 'arg_meshing', Para.arg_meshing = sscanf(remain,'%d%f')';
        case 'arg_fillholes', Para.arg_fillholes = sscanf(remain,'%d');
        
        % Display
        case 'out_AreaVol', Para.out_AreaVol = sscanf(remain,'%d');
        
        % Outputs    
        case 'out_MolStrc', Para.out_MolStrc = sscanf(remain,'%d');
        case 'out_STL', Para.out_STL = sscanf(remain,'%d'); 
            
        % Visualization
        case 'arg_viz', Para.arg_viz = sscanf(remain,'%d');
        case 'arg_ext', Para.arg_ext = sscanf(remain,'%d');
        case 'arg_int', Para.arg_int = sscanf(remain,'%d');
        
        case 'viz_SAS', Para.viz_SAS = sscanf(remain,'%d%f')';
        case 'viz_SES', Para.viz_SES = sscanf(remain,'%d%f')';
            
        case 'viz_sing', Para.viz_sing = sscanf(remain,'%d');
        case 'viz_CirSeg', Para.viz_CirSeg = sscanf(remain,'%d%d')';

        otherwise flag = 0;
    end
    
    if flag == 1
        fprintf([token,': ',num2str(eval(['Para.',token])),'\n']);
    end
    
    line = fgetl(fid);
end

fclose(fid);


end