
clear all
close all


image_name = [

   '5dma_iso3200_1_real.png';

   '5dma_iso3200_2_real.png';

   '5dma_iso3200_3_real.png';

   'd600_iso3200_1_real.png';

   'd600_iso3200_2_real.png';

   'd600_iso3200_3_real.png';

   'd800_iso1600_1_real.png';
 
   'd800_iso1600_2_real.png';

   'd800_iso1600_3_real.png';

   'd800_iso3200_1_real.png';

   'd800_iso3200_2_real.png';
     
   'd800_iso3200_3_real.png';

   'd800_iso6400_1_real.png';

   'd800_iso6400_2_real.png';

   'd800_iso6400_3_real.png'

    ];

im_noise = cell(15,1);
Idenoised = cell(15,1);
im_mean = cell(15,1);

tic,
parfor ii = 1 :15
   im_noise{ii,1} = im2double(imread(image_name(ii,:)));  %% read a noisy image
    
   Idenoised{ii,1} = NLH_CC(im_noise{ii,1});
   
end
toc,

PSNR_SUM = 0.0;
MSSIM_SUM = 0.0;
for ii = 1:15
 im3 = im2double(imread( fullfile('Z:\NLH_CC\real_images\',[image_name(ii,1:end-8),'mean.png'])));  
 imr_final = Idenoised{ii,1};
 PSNR = 10*log10(1/mean((im3(:)-double(imr_final(:))).^2));
 figure,imshow(double(imr_final));
 PSNR_SUM = PSNR_SUM + PSNR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MMSIM value
K = [0.01 0.03];
window = fspecial('gaussian', 11, 1.5);
L = 1;
[mssim1, ssim_map1] = ssim_index(imr_final(:,:,1),im3(:,:,1),K,window,L);
[mssim2, ssim_map2] = ssim_index(imr_final(:,:,2),im3(:,:,2),K,window,L);
[mssim3, ssim_map3] = ssim_index(imr_final(:,:,3),im3(:,:,3),K,window,L);
mssim = (mssim1+mssim2+mssim3)/3.0;
fprintf('Final Estimation, PSNR: %.2f dB, MSSIM: %.4f \n',PSNR, mssim);   
MSSIM_SUM = MSSIM_SUM + mssim;    
end

PSNR_AVE = PSNR_SUM/15
SSIM_AVE = MSSIM_SUM/15

   
   
   
   
   
   
   
  
   
   
   
   
   
   
   
   
   