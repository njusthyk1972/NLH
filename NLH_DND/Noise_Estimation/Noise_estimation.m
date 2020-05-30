
function  sigma_est = Noise_estimation(im)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Original Noise level  estimation

[dis_map,sigma_tb] = NL_distance(8,16,2,39,single(32),double(im));

sigma_est = double(sigma_tb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refined Noise level estimation

if sigma_tb > 70.0 
[dis_map1,sigma_tb1] = NL_distance_refine(8,16,2,39,single(32),double(im),double(dis_map),double(sigma_tb));
sigma_est  = double(sigma_tb1);
elseif sigma_tb <= 70.0
[dis_map1,sigma_tb1] = NL_distance_refine(8,16,2,39,single(32),double(im),double(dis_map),double(sigma_tb));
[dis_map1,sigma_est] = NL_distance_refine(8,16,2,39,single(32),double(im),double(dis_map1),double(sigma_tb1));
end














