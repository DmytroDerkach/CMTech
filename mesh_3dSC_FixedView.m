function SCFV =...
    mesh_3dSC_FixedView( myMesh, v_idxs, ref_dir, ref_normal, r_min, r_MaX, r_dens, nBinsEAR )    
%
% Compute 3dSC descriptors as defined in [Frome et al, 2004] from a fixed 
% viewpoint
% This function computes the descriptor for all vertices indicated in
% 'v_idxs' vector of the given mesh
%
% Syntax :
% 1) ShapeContext3D( myMesh, v_idxs, ref_dir, ref_normal, r_min, r_MaX, r_dens,
% nBinsEAR)
%
% myMesh     Is a structure that describes the triangulated mesh with fields
%            'verts' (of size 3 x num_of_points) and 'faces' (of size 3 x
%            num_of_triangles)
% v_idxs     The vertex indices for which to compute the descriptors
%            Use [] if you want descriptors for all mesh vertices
% ref_dir    The reference direction to set the Azimuth origin 
% kV_normal  The normal vector (common to all points, i.e. cammera viewpt)
% r_min      The minumum radius for the descriptor bins
% r_MaX      The maximum radius for the descriptor bins (controls the
%            neighborhood size - Typically r_MaX / r_min is around 20 or 30)
% r_dens     The radius to compute the sampling density 
% nBinsEAR   The number of bins for elevation, azimuth and radius, in that
%            order. The default from Frome et al. is nBinsEAR = [11,12,15]
%
% Copyright (C) 2017 Federico Sukno, Pompeu Fabra University
%

[one3dSC, common3dSC] = ShapeContext3D( ...
    myMesh, v_idxs(1), ref_dir, ref_normal, ...
    r_min, r_MaX, r_dens, nBinsEAR, []);

t = one3dSC(:);
SCFV = zeros( length( v_idxs ), length( t ));
SCFV(1,:) = t;

fprintf(1, '\tComputing vertex descriptors ');
for j = 2 : length( v_idxs )
    if mod( j, ceil(length( v_idxs ) / 40) ) == 0
        fprintf(1, '.');
    end
%     try
    t = ShapeContext3D( ...
        myMesh, v_idxs(j), ref_dir, ref_normal, ...
        r_min, r_MaX, r_dens, nBinsEAR, common3dSC );    
%     catch
%         disp('')
%     end
    SCFV(j, :) = t(:);
end
fprintf(1,' \n');


