
function  sigma_est = Noise_estimation_1(im)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Original Noise level  estimation

[dis_map,sigma_tb] = NL_distance(8,16,2,39,single(32),double(im));

sigma_est = double(sigma_tb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refined Noise level estimation

end














