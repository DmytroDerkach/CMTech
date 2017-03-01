function [pcloud, distance] = depthToCloud(depth, topleft)
% depthToCloud.m - Convert depth image into 3D point cloud
% Author: Liefeng Bo and Kevin Lai
%
% Input: 
% depth - the depth image
% topleft - the position of the top-left corner of depth in the original depth image. Assumes depth is uncropped if this is not provided
%
% Output:
% pcloud - the point cloud, where each channel is the x, y, and z euclidean coordinates respectively. Missing values are NaN.
% distance - euclidean distance from the sensor to each point
%

if nargin < 2
    topleft = [1 1];
end

depth= double(depth);
depth(depth == 0) = nan;

% RGB-D camera constants
%-----------------------------------------
[imh, imw] = size(depth);
center = [imw / 2 imh / 2]; %[320 240] [263 222]
constant = 365.7;%570.3;
%-----------------------------------------


MM_PER_M = 1000;

% convert depth image to 3d point clouds
pcloud = zeros(imh,imw,3);
xgrid = ones(imh,1)*(1:imw) + (topleft(1)-1) - center(1);
ygrid = (1:imh)'*ones(1,imw) + (topleft(2)-1) - center(2);
pcloud(:,:,1) = xgrid.*depth/constant/MM_PER_M;
pcloud(:,:,2) = ygrid.*depth/constant/MM_PER_M;
pcloud(:,:,3) = depth/MM_PER_M;
distance = sqrt(sum(pcloud.^2,3));

%% does not work
% 
%  center = [263 230];
% depth= double(depth);
% depth(depth == 0) = nan;
% [imh, imw] = size(depth);
% pcloud = zeros(imh,imw,3);
% 
% MM_PER_M = 1000;
% xgrid = ones(imh,1)*(1:imw) - center(1);%(imw-1)/2;
% ygrid = (1:imh)'*ones(1,imw) - center(2);%(imh-1)/2;
% %ygrid =  - ygrid ;% (imh-1)/2;
% 
% pcloud(:,:,1) = xgrid.* depth/MM_PER_M / (center(1)*2)  * tan(degtorad(70/2)); %
% pcloud(:,:,2) = ygrid.* depth/MM_PER_M / (center(2)*2)  * tan(degtorad(60/2)); %
% pcloud(:,:,3) = depth/MM_PER_M;
% distance = sqrt(sum(pcloud.^2,3));
% 
