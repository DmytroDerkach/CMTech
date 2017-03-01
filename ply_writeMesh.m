function ply_writeMesh ( myMesh, fileName, varargin )
%
% PLY_WRITEMESH writes 3D data as a PLY file.
%  This is just a wrapper function that accepts a mesh structure as input
%  (containing fields .verts and .faces as 3*Nv and 3*nF matrices.
%
%  PLY_WRITEMESH (MYMESH, FILENAME)
%  PLY_WRITEMESH (MYMESH, FILENAME, FORMAT)
%       FORMAT can be ascii (default), binary_little_endian or
%       binary_big_endian
%
%  See PLY_WRITE for further details
%

writeFormat = 'ascii';
if not( isempty( varargin ))
    if length( varargin ) > 1
        error('Too many input parameters');
    end
    writeFormat = varargin{1};
end

%   A common PLY data structure has the following fields:
%      DATA.vertex.x = x coordinates, [Nx1] real array
%      DATA.vertex.y = y coordinates, [Nx1] real array
%      DATA.vertex.z = z coordinates, [Nx1] real array
DATA.vertex.x = myMesh.verts(1, :)';
DATA.vertex.y = myMesh.verts(2, :)';
DATA.vertex.z = myMesh.verts(3, :)';

%      DATA.face.vertex_indices = vertex index lists,
%         an {Mx1} cell array where each cell holds a one-
%         dimesional array (of any length) of vertex indices.
% (remember that indexes need to be (-1) shifted for C)
DATA.face.vertex_indices = mat2cell(...
    myMesh.faces - 1, 3, ones(1, size(myMesh.faces, 2)))';

ply_write(DATA, fileName, writeFormat);
