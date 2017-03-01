function [lmk_coords, f3d_idxs, lmk_name, found_flag] = Read_dFL_Landmarks( fileName )
    
[fid, msg] = fopen( fileName, 'rt' );
if fid == -1
    error( msg );
end

nLandmarks = 0;
lmk_coords = [];
f3d_idxs = [];
found_flag = [];
lmk_name = cell(0);

while not( feof( fid ))
    newLine = fgetl(fid);
    while newLine(1) == ' '        
        newLine(1) = [];                  
    end
    
    % Skip comment lines
    if newLine(1) == '@'
        continue;
    end
    
    if newLine(1) == 'L'
        nLandmarks = nLandmarks + 1;        
        
        % Lxx - face 3D id
        [token, remain] = strtok( newLine );    
        f3d_idxs = [f3d_idxs, str2num( token(2:end) )];
    
        % x \t y \t z
        new_landmark = zeros(3, 1);
        [token, remain] = strtok(remain);
        new_landmark(1) = str2num( token );
        [token, remain] = strtok(remain);
        new_landmark(2) = str2num( token );
        [token, remain] = strtok(remain);
        new_landmark(3) = str2num( token );
        lmk_coords = [lmk_coords, new_landmark];
        
        % name \t found/inferred
        [token, remain] = strtok(remain);
        lmk_name{ nLandmarks } = token;
        if strfind( remain, 'FOUND' )
            found_flag = [found_flag, 1];
        else
            found_flag = [found_flag, 0];
        end
    end
end

fclose( fid );
