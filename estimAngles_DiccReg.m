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

