function [my3dSC, common3dSC] = ...
    ShapeContext3D( myMesh, kV, ref_dir, kV_normal, r_min, r_MaX, r_dens, nBinsEAR,...
    common3dSC, varargin)
%
% Compute 3dSC descriptor as defined in [Frome et al, 2004] addiont
% tri-linear interpolation to compute the histogram bins.
% This function computes the descriptor for a single vertex of the mesh.
%
% Syntax options:
% 1) ShapeContext3D( myMesh, kV, ref_dir, kV_normal, r_min, r_MaX, r_dens,
% nBinsEAR, common3dSC)
% 2) ShapeContext3D( myMesh, kV, ref_dir, kV_normal, r_min, r_MaX, r_dens,
% nBinsEAR, common3dSC, external_kdTree)
%
% myMesh     Is a structure that describes the triangulated mesh with fields
%            'verts' (of size 3 x num_of_points) and 'faces' (of size 3 x
%            num_of_triangles)
% kV         The vertex number for which we want to compute the descriptor
% ref_dir    The reference direction to set the Azimuth origin (in principle
%            it should be set randomly, unless there is extra information
%            that can be used)
% kV_normal  The normal verctor to the vertex kV
% r_min      The minumum radius for the descriptor bins
% r_MaX      The maximum radius for the descriptor bins (controls the
%            neighborhood size - Typically r_MaX / r_min is around 20 or 30)
% r_dens     The radius to compute the sampling density 
% nBinsEAR   The number of bins for elevation, azimuth and radius, in that
%            order. The default from Frome et al. is nBinsEAR = [11,12,15]
% common3dSC The common part that could be pre-computed for a given mesh
%            and descriptor parameters (i.e. it is the same for all
%            vertices of the mesh). If not available, specify [] instead
%            and it will be computed internally (and returned at the output).
%            It is highly recomended to avoid computing it for every new
%            vertex.
% use_kd_tree Flag to enable the use of a kdTree 
%
% Copyright (C) 2017 Federico Sukno, Pompeu Fabra University
%

% Optional: use a kd-tree
use_kd_tree = true;



% The common info for 3dSC computation can be passed externally
if not( isempty( varargin ))    
    myKdTree = varargin{1};
    use_kd_tree = true;
    if isempty( common3dSC )
        common3dSC = ShapeContext_PreCompute(...
            myMesh, r_min, r_MaX, r_dens, nBinsEAR, true, myKdTree);
    end
else
    if isempty( common3dSC )
        common3dSC = ShapeContext_PreCompute(...
            myMesh, r_min, r_MaX, r_dens, nBinsEAR, true);        
    end
    myKdTree = common3dSC.KdTree;
end
    
% Some pre-definitions for speed-up
myVerts = myMesh.verts;
SC_radialIntervals = common3dSC(1).radialIntervals;
SC_elevation_interval = common3dSC(1).elevation_interval;
SC_azimuth_interval = common3dSC(1).azimuth_interval;
SC_radialInt_invSize = common3dSC(1).radialInt_invSize;

% Initialize descriptor
my3dSC = zeros( nBinsEAR([2,1,3]) );
NV = size( myVerts, 2);


% Compute the X-axis orientation by computing the input direction into the
% plane defined by the normal
% ------------------------------------------------------------------------
% - The x_dir must have null dot product with N to belong to the plane
% Thus x1 n1 + x2 n2 + x3 n3 = 0 => There are 3 options
% - But we want to respect the input dir, 
% Hence we want to minimize || x_dir - ref_dir ||
% - A simple soltion is to compute all 3 options and test how much change
% of direction there is w.r.t. the ref_direction

x_options = zeros( 3 );
ref_dir = ref_dir / norm( ref_dir );

if abs( kV_normal(3) ) > 1e-6
    x_options(:, 1) = [ref_dir(1:2);...
        (-kV_normal(1:2)' * ref_dir(1:2)) / kV_normal(3)];
    x_options(:, 1) = x_options(:, 1) / norm( x_options(:, 1));
else
    x_options(:, 1) = [0 0 1]';
end

if abs( kV_normal(2) ) > 1e-6
    x_options(:, 2) = [ref_dir(1);...
        (-kV_normal([1,3])' * ref_dir([1,3])) / kV_normal(2);...
        ref_dir(3)];
    x_options(:, 2) = x_options(:, 2) / norm( x_options(:, 2));
else
    x_options(:, 2) = [0 1 0]';
end

if abs( kV_normal(1) ) > 1e-6
    x_options(:, 3) = [...
        (-kV_normal(2:3)' * ref_dir(2:3)) / kV_normal(1),...
        ref_dir(2:3)'];
    x_options(:, 3) = x_options(:, 3) / norm( x_options(:, 3));
else
    x_options(:, 3) = [1 0 0]';
end


% Since we have normalized the vectors, the dot product of each option with
% the ref_direction will directly tell the cosine of the difference angle.
% Then, the highest dot product os the one we're looking for
[nada, best_opt] = max( dot( x_options, repmat(ref_dir, [1 3]) ));
x_axis = x_options(:, best_opt);


% Compute the descriptor
% ------------------------------------------------------------------------
vertWeights = common3dSC(1).vertDensity;

% Get the neighbors
if use_kd_tree
    %v_idxs = kdtree_ball_query(...
    %    myKdTree, myVerts(:, kV), r_MaX);
    v_idxs = rangesearch(myKdTree, myVerts(:, kV)', r_MaX);
    v_idxs = v_idxs{:};
else        
    dd_v = myVerts - repmat(myVerts(:,kV), [1 NV]);
    v_idxs = find( sum( dd_v.^2, 1 ) <= r_MaX^2);
end

if length( v_idxs ) > 1

    for j_idx = 1 : length( v_idxs )
       jv = v_idxs( j_idx );
       if jv == kV
           continue;
       end
       
       diffV = myVerts(:, jv) - myVerts(:, kV);
       rDiff = norm( diffV );
       
       N_dot_diffV = diffV(1)*kV_normal(1) +...
           diffV(2)*kV_normal(2) + diffV(3)*kV_normal(3);
            
       tpV = myVerts(:, jv) - N_dot_diffV .* kV_normal -...
           myVerts(:, kV);
       
       % The norm, if no repeated points, cannot be zero
       tpV = tpV / norm( tpV );
       
       % Azimuth (angle of tpV with X-axis)
       cp = [x_axis(2).*tpV(3)- x_axis(3).*tpV(2);
           x_axis(3).*tpV(1)-x_axis(1).*tpV(3);
           x_axis(1).*tpV(2)-x_axis(2).*tpV(1)];
            
       dot_xAxis_tpV = x_axis(1) * tpV(1) + ...
           x_axis(2) * tpV(2) + x_axis(3) * tpV(3);       
       
       phi = atan2 (norm( cp ), dot_xAxis_tpV);
       
       N_dot_cp = kV_normal(1) * cp(1) + ...
           kV_normal(2) * cp(2) + kV_normal(3) * cp(3);              
       if N_dot_cp < 0
           phi = 2*pi - phi;
       end
       
       % Elevation (angle to the normal) in [0, PI]
       theta = acos (min (1.0,... 
            max( -1.0, N_dot_diffV / rDiff )));

       % Determine the bin coordinates
       % ------------------------------------------------------------------
       jRad = 1;
       for rad = 2 : nBinsEAR(3) + 1
           if rDiff <= SC_radialIntervals( rad )
               jRad = rad - 1;
               break;
           end
       end
       
       xElev = theta / SC_elevation_interval;
       jElev = 1 + min( floor( xElev ), nBinsEAR(1) - 1 );
       
       xAzim = phi / SC_azimuth_interval;
       jAzim = 1 + min( floor( xAzim ), nBinsEAR(2) - 1);
       
       % Neighbors for interpolation
       jR2 = min( jRad + 1, nBinsEAR(3) );
       jE2 = jElev + 1;
       if jE2 > nBinsEAR(1)
           jE2 = 1;
       end
       
       jA2 = jAzim + 1;
       if jA2 > nBinsEAR(2)
           jA2 = 1;
       end
       
       % Deltas       
       dE = 1 + xElev - jElev;
       dA = 1 + xAzim - jAzim;
       dR = 0;
       if rDiff > SC_radialIntervals(1)
           dR = (rDiff - SC_radialIntervals( jRad )) *...
               SC_radialInt_invSize( jRad );
       end
       
       % Update weights
       my3dSC( jAzim, jElev, jRad ) = my3dSC( jAzim, jElev, jRad ) + ...
           (1-dA) * (1-dE) * (1-dR) * vertWeights(jv);
       my3dSC( jAzim, jElev, jR2 ) = my3dSC( jAzim, jElev, jR2 ) + ...
           (1-dA) * (1-dE) * (dR) * vertWeights(jv);       
       
       my3dSC( jA2, jElev, jRad ) = my3dSC( jA2, jElev, jRad ) + ...
           (dA) * (1-dE) * (1-dR) * vertWeights(jv);
       my3dSC( jA2, jElev, jR2 ) = my3dSC( jA2, jElev, jR2 ) + ...
           (dA) * (1-dE) * (dR) * vertWeights(jv);       
       
       my3dSC( jAzim, jE2, jRad ) = my3dSC( jAzim, jE2, jRad ) + ...
           (1-dA) * (dE) * (1-dR) * vertWeights(jv);
       my3dSC( jAzim, jE2, jR2 ) = my3dSC( jAzim, jE2, jR2 ) + ...
           (1-dA) * (dE) * (dR) * vertWeights(jv);    
       
       my3dSC( jA2, jE2, jRad ) = my3dSC( jA2, jE2, jRad ) + ...
           (dA) * (dE) * (1-dR) * vertWeights(jv);
       my3dSC( jA2, jE2, jR2 ) = my3dSC( jA2, jE2, jR2 ) + ...
           (dA) * (dE) * (dR) * vertWeights(jv);    
        
    end    
end

% Normaliztion
my3dSC = my3dSC .* common3dSC.invBinVol;

