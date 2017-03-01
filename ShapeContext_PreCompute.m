function common3dSC = ShapeContext_PreCompute(...
        myMesh, r_min, r_MaX, r_density, nBinsEAR, use_kd_tree, varargin)
%
% common3dSC = ShapeContext_PreCompute(...
%      myMesh, r_min, r_MaX, r_density, nBinsEAR, use_kd_tree)
%
% or
%
% common3dSC = ShapeContext_PreCompute(...
%      myMesh, r_min, r_MaX, r_density, nBinsEAR, use_kd_tree, exter_tree)
%       
% Pre-compute the common information needed for all vertices in the
% computation of a 3D-Shape Context [Frome et al, 2004]
%
% myMesh     Is a structure that describes the triangulated mesh with fields
%            'verts' (of size 3 x num_of_points) and 'faces' (of size 3 x
%            num_of_triangles)
% r_min      The minumum radius for the descriptor bins
% r_MaX      The maximum radius for the descriptor bins (controls the
%            neighborhood size - Typically r_MaX / r_min is around 20 or 30)
% r_dens     The radius to compute the sampling density 
% nBinsEAR   The number of bins for elevation, azimuth and radius, in that
%            order. The default from Frome et al. is nBinsEAR = [11,12,15]
% use_kd_tree Flag to enable the use of a kdTree 
%
% exter_tree is an optional parameter to specify an external kdTree %
%
% Copyright (C) 2017 Federico Sukno, Pompeu Fabra University
%

% Optional use of a kd-tree
myVerts = myMesh.verts;
if use_kd_tree
    if isempty( varargin )
        %fprintf(1, '\tBuilding kdtree...');        
        myKdTree = createns( myVerts','nsmethod','kdtree');
        %fprintf(1, '\n');
    else
        myKdTree = varargin{1};
    end        
end


% Initialize
% ---------------------------------------------------------------------
NV = size( myVerts, 2);

% 3dSC parameters
common3dSC = struct();
common3dSC(1).r_density = r_density;
common3dSC(1).r_min = r_min;
common3dSC(1).r_MaX = r_MaX;
common3dSC(1).nBinsEAR = nBinsEAR;

if use_kd_tree    
    common3dSC(1).KdTree = myKdTree;    
end

% Auxiliar variabiles
common3dSC(1).azimuth_interval = 2*pi / nBinsEAR(2);
common3dSC(1).elevation_interval = pi / nBinsEAR(1);
common3dSC(1).radialIntervals = exp(...
    log( r_min ) + (0 : nBinsEAR(3)) * log( r_MaX / r_min ) / nBinsEAR(3));
common3dSC(1).radialInt_invSize = 1.0 ./ (...
    common3dSC(1).radialIntervals(2 : end) -...
    common3dSC(1).radialIntervals(1 : end-1));
common3dSC(1).thetaDiv = (0 : nBinsEAR(1) + 1) *...
    common3dSC(1).elevation_interval;

% Inverse bin volume
common3dSC(1).invBinVol = zeros( nBinsEAR([2,1,3]) );
for jR = 1 : nBinsEAR(3)
    integr_r = (common3dSC(1).radialIntervals(jR+1) ^ 3) / 3.0 -...
        (common3dSC(1).radialIntervals(jR) ^ 3) / 3.0;

	for kE = 1 : nBinsEAR(1)
		integr_theta = cos( common3dSC(1).thetaDiv(kE) ) -...
            cos ( common3dSC(1).thetaDiv(kE+1) );
        theVolume = common3dSC(1).azimuth_interval * integr_theta * integr_r;
	
        % Compute cube root of the computed volume
        common3dSC(1).invBinVol( :, kE, jR ) =...
            ones( nBinsEAR(2), 1 ) / ( theVolume ^ (1/3) );

    end
end



        
% Vertex density
% ---------------------------------------------------------------------
if r_density > 0    
    common3dSC(1).vertDensity = zeros(NV, 1);
    fprintf(1, '\tComputing vertex density ');
    for jV = 1 : NV
        if mod( jV, ceil(NV / 40) ) == 0
            fprintf(1, '.');
        end

        if use_kd_tree
            v_idxs = rangesearch(myKdTree, myVerts(:, jV)', r_density);
            v_idxs = v_idxs{:};
        else        
            dd_v = myVerts - repmat(myVerts(:,jV), [1 NV]);
            v_idxs = find( sum( dd_v.^2, 1 ) <= r_density^2);
        end

        nNeigs = length( v_idxs );
        if length( v_idxs ) < 2
            % Not enough neighbors
            common3dSC(1).vertDensity(jV) = 1;
        else
            common3dSC(1).vertDensity(jV) = 1.0 / nNeigs;
        end        
    end
else
    common3dSC(1).vertDensity = ones(NV, 1);
end

fprintf(1,' \n');
