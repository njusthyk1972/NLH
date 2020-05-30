
clear all
close all

clc


image_name = [
%     'Cameraman256.png'
    'house.png'
%     'peppers256.png'
%     'montage.png'
%     'Lena512.png'
%     'barbara.png'
%     'boat.png'
%     'man.png'
%     'couple.png'
%     'hill.png' 
    ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Specify the std. dev. of the corrupting noise
%%%%
if (exist('sigma') ~= 1),
    sigma               = 50; %% default standard deviation of the AWGN
end


    im = im2double(imread(image_name));  %% read a noise-free image and put in intensity range [0,1]

    im1= im;

    [M,N]=size(im);
    
     
    randn('seed', 0);                          %% generate seed
    im = im + (sigma/255)*randn(size(im)); %% create a noisy image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure,imshow(im);
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation

[sigma_tb sigma_est dis_map] = Noise_estimation(im); %% estimate noisy level

sigma_est = double(sigma_est)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dis_map,sigma_tb1] = NL_distance(8,16,2,39,single(8),double(im));
dis_map = double(dis_map);
%%%%%%%%%%%%%%%%%% 
Ns     = 43;
N3     = 4;
N2     = 16;

%%%%%%%%%%%%%%%%%
tic,


%%%%%%% Main function

imr = im;  

alpha = 0.618; 

lamda = 0.8;

Thr = 1.45;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   iteration times
if sigma_est < 7.5
    k = 2;
    N_step = 5;
    N1     = 9;
elseif sigma_est >= 7.5 && sigma_est < 12.5
    k = 3;
    N_step = 5;
    N1     = 9;
 elseif sigma_est >= 12.5 && sigma_est < 35
    k = 3;
    N_step = 6;
    N1     = 9;
 elseif sigma_est >= 35 && sigma_est < 55
    k = 4;
    N_step = 6;
    N1     = 9;
  elseif sigma_est >= 55 && sigma_est < 75
      k = 5;
     N_step = 7;
     N1     = 10;
 elseif sigma_est >= 75 && sigma_est < 85
     k = 6;
      N_step = 8;
      N1     = 11;
else 
      k = 7;
      N_step = 9;
      N1     = 11;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ll = 1:k 

     
        imr = NLH_AWGN_Gray(N1,N2,N3,Ns,N_step,double(alpha*imr+(1-alpha)*im),single(Thr),...
        sigma_est, double(lamda*imr+(1-lamda)*im));
     
  
          PSNR = 10*log10(1/mean((im1(:)-double(imr(:))).^2));
 
         figure,imshow(double(imr));

         fprintf( 'Iter = %2.0f, PSNR = %2.4f \n\n\n', ll, PSNR );

         N_step = N_step -1;
         if N_step <= 3
             N_step = 3;
         end
 
end 



toc,


imr = double(imr);


PSNR = 10*log10(1/mean((im1(:)-imr(:)).^2));
 figure,imshow(imr);title(sprintf('Denoised %s, PSNR: %.3f dB', ...
        image_name(1:end-4), PSNR));
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation

[sigma_tb1 sigma_est1 dis_map1] = Noise_estimation(imr); %% estimate noisy level

sigma_est1 = double(sigma_est1)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% NLH_wiener
N2 = 64;
N3 = 8;
Ns = 129;


if sigma_est < 25
    N1      = 8;
    N_step1 = 8;
    N_step2 = 5;
elseif sigma_est >= 25 && sigma_est < 75
    N1      = 16;
    N_step1 = 16;
    N_step2 = 13;    
else
    N1      = 20;
    N_step1 = 20;
    N_step2 = 13;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sigma_est < 7.5
gamma = 2.0 + abs(sigma_est1-sigma_est*0.05); 
elseif sigma_est >= 7.5 && sigma_est < 25
gamma = 2.0 + abs(sigma_est1-sigma_est*0.05); 
else
gamma = 2.0 + abs(sigma_est1-sigma_est*0.05);
end

beta = 0.8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stage two: Wiener filtering

tic,
       
         y_est = NLH_AWGN_Wiener_Gray(N1,N2,N3,Ns,N_step1,double(im),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(imr*beta+im*(1-beta)));            
    
                y_est=double(y_est);
         
         
          PSNR = 10*log10(1/mean((im1(:)-y_est(:)).^2));        
          figure,imshow(y_est);title(sprintf('Denoised  PSNR: %.3f dB',PSNR));
 
                       
        y_est = NLH_AWGN_Wiener_Gray(N1,N2,N3,Ns,N_step2,double(y_est),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(y_est*beta+im*(1-beta)));
     
         
         PSNR = 10*log10(1/mean((im1(:)-y_est(:)).^2));        
         fprintf(sprintf('The final result , PSNR: %.2f dB\n',  PSNR));
         
                   
      y_est=double(y_est);
      figure,imshow(y_est);title(sprintf('Denoised  PSNR: %.3f dB',PSNR));

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MMSIM value
K = [0.01 0.03];
window = fspecial('gaussian', 11, 1.5);
L = 1;
[mssim, ssim_map] = ssim_index(im1,y_est,K,window,L);

fprintf('FINAL ESTIMATE, MSSIM: %.4f dB\n', mssim);

% figure,imshow(ssim_map);


 














