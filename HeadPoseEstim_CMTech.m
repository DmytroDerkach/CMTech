function [pitch, yaw, roll] = HeadPoseEstim_CMTech( npy_fileName, varargin )
% 
% Default:
% [pitch, yaw, roll] = HeadPoseEstim_CMTech( npy_fileName );
% 
% It is also possible to specify the number of cores/processors used (if
% not specified, default is 4).
% [pitch, yaw, roll] = HeadPoseEstim_CMTech( npy_fileName, nCores );
%

    pitch = NaN;
    yaw = NaN;
    roll = NaN;
    
    % Parameters
    maxPROC = 4;
    if not( isempty( varargin ))
        if not( length( varargin ) == 1 )
            error('Too many parameters');
        end
        maxPROC = varargin{1};
        if maxPROC < 1
            error('At least one processor is required');
        end
    end
    maxSearchTIME = max(20, 45 * 4 / maxPROC);

    % Preliminary checks    
    cfg_file = 'srilf3dFL_F1_K2_155_BQ3.cfg';
    if not( exist( cfg_file, 'file' ))
        fprintf(1, '\n%s file missing\n', cfg_file);
        error('Cannot find landmark configuration file');
    end
    if not( exist('vtkClearDisconnecteedPts.exe', 'file'))
        error('Cannot find vtkClearDisconnecteedPts.exe');
    end
    if not( exist('srilf3dFL.exe'))
        fprintf(1, '\n*** SRILF landmarking software missing ***');
        fprintf(1, '\n*** Please download SRILF 3D Face Landmarker from\n');
        fprintf(1, '%s\n', 'http://fsukno.atspace.eu/Data.htm');
        error('srilf3dFL.exe missing');
    end
    
    % Load landmark regressors
    load LmkRegress.mat;       
    
    % Load input data
    fprintf(1, '\nLoading: %s\n', npy_fileName);
    if not( exist( npy_fileName, 'file' ))        
        error('Cannot find input file');
    end
        
    myMesh = [];
    for n_try = 1 : 4
        try 
            myMesh = extract_HeadMesh_from_npy( npy_fileName );
            break;
        catch
            if n_try > 3
                warning('Cannot build mesh from %s\n', npy_fileName);
            end
            pause( 1 );            
        end
    end
    
    if isempty( myMesh )
        warning('Estimation failed');    
        pitch = 0;
        yaw = 0;
        roll = 0;
        return;
    end
    
    % Auxiliary file names
    last_bar = find( npy_fileName == '\' | npy_fileName == '/', 1, 'last');
    if isempty( last_bar )
        last_bar = 0;
    end
    last_pt = find( npy_fileName == '.', 1, 'last');
    if isempty( last_pt )
        last_pt = length( npy_fileName ) + 1;
    end
    short_name = npy_fileName( last_bar + 1 : last_pt - 1 );
    
    ply_name = sprintf('%s.ply', short_name);
    dFL_name = sprintf('%s.dFL', short_name);
    
    % Landmark detection
    r_val = ExportAndLocateLandmarks( myMesh, ...
        ply_name, dFL_name, cfg_file, maxSearchTIME, maxPROC);
    if r_val == 1
        warning('Using pre-saved landmarks. Make sure they are up-to-date\n');
    end

    
    landmark_estim_ok = false;
    if exist( dFL_name, 'file' )
        % Get coordinates
        lmk_coords = Read_dFL_Landmarks( dFL_name );
        
        % Estimate angles
        [pG, yG, rG] = estimateAnglesFromGeometry( lmk_coords );    
        [pD, yD, rD] = estimateAnglesFromDescr( lmk_coords, ...
            myMesh, LmkRegress);
        
        if abs( pG - pD ) + abs( yG - yD ) + abs( rG - rD ) < 50
            pitch = .5 * (pG + pD);
            yaw = .5 * (yG + yD);
            roll = .5 * (rG + rD);
            landmark_estim_ok = true;
        else
            landmark_estim_ok = false;
        end
    end

    if not( landmark_estim_ok )        
        fprintf(1, '\nDictionary-regression');
        load DictionaryReg.mat
        
        theIndices = mesh_umbrellaSample (myMesh, 7);
        SCFV_512 = mesh_3dSC_FixedView( myMesh, theIndices,...
            [1 0 0]', [0 0 1]', 2, 30, 5, [8 8 8]);        
        [pitch, yaw, roll] = ...
            estimAngles_DiccReg( SCFV_512, centers, dicc_regressor );           
    end
         
end



% ------------------------------------------------------------------
function myMesh = extract_HeadMesh_from_npy( npy_fileName )

    data = readNPY( npy_fileName );
    pcloud = depthToCloud(data);
    
    xx = -1e3 * pcloud(:,:,1);
    yy = -1e3 * pcloud(:,:,2);
    zz = -1e3 * pcloud(:,:,3);
    
    valid_idxs = intersect( find( isfinite( zz )),...
        intersect( find( isfinite( xx )), find( isfinite( yy ))));
    
    xx = xx( valid_idxs );
    yy = yy( valid_idxs );
    zz = zz( valid_idxs );
    
    % We only want the "cluster" that is closest to the camera
    [clus_id, clus_cc] = kmeans( zz, 2 );
    if clus_cc(1) > clus_cc(2)
        valid_idxs = find( clus_id == 1 );
    else
        valid_idxs = find( clus_id == 2 );
    end
    
    xx = xx( valid_idxs );
    yy = yy( valid_idxs );
    zz = zz( valid_idxs );
    
    %  only head
    [clus_id, clus_cc] = kmeans( yy, 2 );
    if clus_cc(1) > clus_cc(2)
        valid_idxs = find( clus_id == 1 );
    else
        valid_idxs = find( clus_id == 2 );
    end
      
    theVerts = [xx( valid_idxs ), yy( valid_idxs ), zz( valid_idxs )]';
    theTriangles = delaunay( xx( valid_idxs ),  yy( valid_idxs ));
    
    % Create the mesh
    rawMesh = struct('verts', theVerts, 'faces', theTriangles');
    myMesh = meshRemoveLargeTriangs( rawMesh, 25 );
    
end

% ------------------------------------------------------------------
function newMesh = meshRemoveLargeTriangs( myMesh, maxTriMaxEdge )
% newMesh = meshRemoveLargeTriangs( myMesh, maxTriMaxEdge )
% 
% Eliminates triangles which have any of its edges larger than the
% specified maxTriMaxEdge
%

    triEdgeLengths2 = zeros( size( myMesh.faces ) );
    for jT = 1 : size( myMesh.faces, 2 )
        triEdgeLengths2(1, jT) = sum((...
            myMesh.verts( :, myMesh.faces(1, jT) ) - ...
            myMesh.verts( :, myMesh.faces(2, jT) )...
            ).^2);
        triEdgeLengths2(2, jT) = sum((...
            myMesh.verts( :, myMesh.faces(2, jT) ) - ...
            myMesh.verts( :, myMesh.faces(3, jT) )...
            ).^2);
        triEdgeLengths2(3, jT) = sum((...
            myMesh.verts( :, myMesh.faces(3, jT) ) - ...
            myMesh.verts( :, myMesh.faces(1, jT) )...
            ).^2);   
    end

    %maxTriMaxEdge = outlier_qThresh( sqrt( triEdgeLengths2(:) ), 10 );
    %maxTriMaxEdge = 30; % mm

    if maxTriMaxEdge == 0
        maxTriMaxEdge2 = outlier_qThresh( triEdgeLengths2(:), 20 );
    else
        maxTriMaxEdge2 = maxTriMaxEdge^2;
    end


    elim_idx = zeros(1, size( myMesh.faces, 1));
    n_to_elim = 0;
    for jT = 1 : size( myMesh.faces, 2 )
        if max( triEdgeLengths2(:, jT) ) > maxTriMaxEdge2
            n_to_elim = n_to_elim + 1;
            elim_idx( n_to_elim ) = jT;
        end
    end

    newMesh = myMesh;
    newMesh.faces( :, elim_idx( 1 : n_to_elim )) = [];

end


% ------------------------------------------------------------------
function r_val = ExportAndLocateLandmarks( myMesh, ...
    ply_name, dFL_name, cfg_file, maxSearchTIME, maxPROC)
    
    if exist( dFL_name, 'file' )
        r_val = 1;
        return;
    end

    fprintf(1, '\nExporting to PLY format ... ');
    for n_try = 1 : 2
        try            
            ply_writeMesh( myMesh, ply_name );
            fprintf(1, 'ok\n');
            break;
        catch
            if n_try > 1
                warning('Cannot write mesh to %s\n', ply_name);
            end
            pause( 1 );
        end
    end

    if not( exist( ply_name, 'file' ))
        r_val = 2;
        return;
    end

    fprintf(1, 'Detecting landmarks...');
    clean_mesh = sprintf('vtkClearDisconnecteedPts.exe %s',...
        ply_name);
    dos( clean_mesh );                

    find_landmarks = sprintf(...
        'srilf3dFL.exe %s %s %s -occ2 10 -forceAnswerYes -maxTime %d -omp %d',...
        cfg_file, ply_name, dFL_name, maxSearchTIME, maxPROC);
    r_val = dos( find_landmarks );
    if not( r_val == 0 )
        if not( r_val == -500 )
            warning('\nThe landmarking software has returned an error');
            fprintf(1, '\nThe most likely reason for this is that either some dLLs');
            fprintf(1, '\nor the Microsoft VC++ 2008 Redistributable package');
            fprintf(1,'\nare missing.');
            fprintf(1, 'Please refer to the technical documentaion of');
            fprintf(1, '\n3D Face Landmarker for details, available at;\n%s\n',...
                'http://fsukno.atspace.eu/srilf3dFL/Documentation_srilf3dFL_v1.0.pdf');
            error('Cannot run SRILF 3D Face Landmarker');
        end        
    end
    delete( ply_name );
  
end



% ------------------------------------------------------------------
function [ pitch, yaw, roll] = estimateAnglesFromGeometry( landmarks )
% landmarks = 3 x n matrix of landmarks coordinates

 indecesForLine = [4,9,5,10]; % coord of eyes line [11 12]
 pointsLine = landmarks(:, indecesForLine)';
% ***** calculate ROLL angle *******************    
    a = mean(pointsLine(1:2,:));
    b = mean(pointsLine(3:4,:));
    lineAB = [a;b];
    f = lineAB(:,1:2);
    f(1,:) = f(1,:) - f(2,:);
    f(2,:) = f(2,:) - f(2,:);
    roll = atand(f(1,2)/f(1,1));
    
     
% ***** calculate PITCH & YAW angle *******************    
    indecesForPlane = [1,3,6:8,11,12];
    pointsPlane = landmarks(:, indecesForPlane)';
    
%     for j = indecesForPlane%1:length(landmarks)
%         plot3(landmarks(1,j), landmarks(2,j), landmarks(3,j), 'or',...
%             'LineWidth', 4); hold on;
%        
%     end
    
    [n, ~, ~] = affine_fit(pointsPlane);
    yaw = atand(n(1)/n(3));
    pitch = -atand(n(2)/n(3));
end

% ------------------------------------------------------------------
function [n,V,p] = affine_fit(X)
    %Computes the plane that fits best (lest square of the normal distance
    %to the plane) a set of sample points.
    %INPUTS:
    %
    %X: a N by 3 matrix where each line is a sample point
    %
    %OUTPUTS:
    %
    %n : a unit (column) vector normal to the plane
    %V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
    %plane
    %p : a point belonging to the plane
    
    
    %the mean of the samples belongs to the plane
    p = mean(X,1);
    
    %The samples are reduced:
    R = bsxfun(@minus,X,p);
    %Computation of the principal directions if the samples cloud
    [V,D] = eig(R'*R);
    %Extract the output from the eigenvectors
    n = V(:,1);
    V = V(:,2:end);
end

% ------------------------------------------------------------------
function [ pitch, yaw, roll] = estimateAnglesFromDescr( ...
    lmk_coords, myMesh, LmkRegress)
    
    % Get closest mesh vertices
    LLverts = zeros( 1, 12 );
    LLonSrf = zeros( 1, 12 );
    for jL = 1 : 12
        dv = myMesh.verts - repmat(...
            lmk_coords(:,jL), [1 size( myMesh.verts, 2 )]);
        [min_d2, LLverts(jL)] = min( sum( dv.^2 ));
        if min_d2 <= 25
            LLonSrf( jL ) = 1;
        end
    end

    % Get the descriptor and regress to get the angles
    Lmk_descrip = mesh_3dSC_FixedView(...
        myMesh, LLverts, [1 0 0]', [0 0 1]', 2, 30, 0, [8 8 8]);

    pitch = 0;
    yaw = 0;
    roll = 0;
    n_adds = 0;
    for jL = 1 : 12        
        if LLonSrf( jL )
            n_adds = n_adds + 1;
            myDescPCA = (Lmk_descrip(jL,:) - LmkRegress(jL).mean) *...
                LmkRegress(jL).FI;
            pitch = pitch + [myDescPCA, 1] * LmkRegress(jL).Pitch;
            yaw = yaw + [myDescPCA, 1] * LmkRegress(jL).Yaw;
            roll = roll + [myDescPCA, 1] * LmkRegress(jL).Roll;
        end
    end
    
    if n_adds
        pitch = pitch / n_adds;
        yaw = yaw / n_adds;
        roll = roll / n_adds;
    end
end


% ------------------------------------------------------------------
function theIndices = mesh_umbrellaSample (myMesh, r)
%
% 

    myKdTree = createns(myMesh.verts','nsmethod','kdtree');

    % Compute spin image for each vertex
    disabeled_flag = zeros( 1, size( myMesh.verts, 2 ));
    for j = 1 : size( myMesh.verts, 2 )
        if disabeled_flag(j) 
            continue;
        end

        pj = myMesh.verts(:, j);

        % Find all points within a ball of size r    
        % v_idxs = kdtree_ball_query(myKdTree, pj, r);
        v_idxs = rangesearch(myKdTree, pj', r);
        disabeled_flag( v_idxs{:} ) = 1;
        disabeled_flag( j ) = 0;

    end

  
   theIndices = find( disabeled_flag == 0 );

end


% ------------------------------------------------------------------
function [pitch, yaw, roll] = ...
    estimAngles_DiccReg( mesh_descriptors, centers_train, dicc_reg )
%

   
  % build hist for input
  d = pdist2(centers_train',mesh_descriptors);
  % exp normalization
  d = exp(-d);  
  d_norm = bsxfun(@rdivide, d, sum(d));
  setHist_mean_test = mean(d_norm');  
  setHist_mean_test = [setHist_mean_test,  ones(size(setHist_mean_test, 1),1)];
  
  predic_Angles = setHist_mean_test * dicc_reg;
  
  pitch = predic_Angles(1);
  yaw = predic_Angles(2);
  roll = predic_Angles(3);
end