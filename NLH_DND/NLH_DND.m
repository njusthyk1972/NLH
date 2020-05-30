
% A blind pixel-level nonlocal method for real-world image denoising on DND
% dataset
%
% denoiser      Function handle
%               It is called as Idenoised = NLH_DND(Inoisy) where Inoisy is the real-world noisy image
%              

% Author: Yingkun Hou, Taishan University (ykhou@tsu.edu.cn)
%         Jun Xu,      Nankai  University
%
% This file is to implement real-world image denoising algorithm as described in the following paper:
% Yingkun Hou, Jun Xu, Mingxia Liu, Guanghai Liu, Li Liu, Fan Zhu, and Ling Shao, "NLH: A Blind Pixel-level
% Non-local Method for Real-world Image Denoising", IEEE Transactions on Image Processing, vol. 29, pp. 5121-5135, 2020.

function imr_final = NLH(Inoisy_crop)

   im_n = Inoisy_crop;
   
   [M1,M2,M3]=size(im_n);
 
%    img=rgb2gray(im_n);
   
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Colorspace in which we perform denoising. BM is applied to the first
%%%%  component and the matching information is re-used for the other two.
%%%%
if (exist('colorspace') ~= 1),
    colorspace              = 'opp'; %%% (valid colorspaces are: 'yCbCr' and 'opp')
end
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change colorspace, compute the l2-norms of the new color channels
%%
[zColSpace l2normLumChrom] = function_rgb2LumChrom(im_n, colorspace);
im=zColSpace; 

% % 
sigma_est = Noise_estimation(im_n(:,:,1)); %% estimate noisy level
sigma_est1 = double(sigma_est);

sigma_est = Noise_estimation(im_n(:,:,2)); %% estimate noisy level
sigma_est2 = double(sigma_est);

sigma_est = Noise_estimation(im_n(:,:,3)); %% estimate noisy level
sigma_est3 = double(sigma_est);



sigma_est = sqrt((sigma_est1^2+sigma_est2^2+sigma_est3^2)/3);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%% 
Ns     = 43;   % Searching neighbour radius
N3     = 4;    % Row number in a similar pixel group
N_step = 6;    % Sliding step
N1     = 6;    % Block size
N2     = 16;   % Image block number in a similar block group
%%%%%%%%%%%%%%%%




beta = 1.0;
alpha = 1.0; 

sigma_est_b = sigma_est*4.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


  
    
    Thr = 2.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
imr1 = imp1;
imr2 = imp2;
imr3 = imp3;

%%%%%%% Main function

% tic,
for ll = 1:2 
        [imr1,imr2,imr3] = NLH_Ruf(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(beta*imr2+(1-beta)*imp2),double(beta*imr3+(1-beta)*imp3),single(Thr),sigma_est_b/255,double(im_n(:,:,2)));
        N_step = 6;
end 

% toc,

imr=zeros(M1,M2,M3);
imr(:,:,1)=imr1(:,:);
imr(:,:,2)=imr2(:,:);
imr(:,:,3)=imr3(:,:);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert back to RGB colorspace
%%%
imr_basic = function_LumChrom2rgb(imr, colorspace);

imr_basic = double(imr_basic);

[dis_map,sigma_tb1] = NL_distance(8,16,2,39,single(8),double(imr_basic(:,:,2)));
dis_map = double(dis_map);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation
sigma_est = Noise_estimation_1(imr_basic(:,:,1)); %% estimate noisy level
sigma_est1 = double(sigma_est);
sigma_est = Noise_estimation_1(imr_basic(:,:,2)); %% estimate noisy level
sigma_est2 = double(sigma_est);
sigma_est = Noise_estimation_1(imr_basic(:,:,3)); %% estimate noisy level
sigma_est3 = double(sigma_est);

sigma_est1 = sqrt((sigma_est1^2+sigma_est2^2+sigma_est3^2)/3)*4.0;


beta = 1.0;
alpha = 0.618; 

sigma_est_b = sigma_est_b - sigma_est1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


   
    Thr = 2.0; % hard-thresholding parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
imr1 = imp1;
imr2 = imp2;
imr3 = imp3;

%%%%%%% Main function

if sigma_est_b < 10
    N1 = 6;
    N_step = 6;
elseif  sigma_est_b >= 10 && sigma_est_b < 20
    N1 = 7;
    N_step = 7;
elseif  sigma_est_b >= 20 && sigma_est_b < 30
    N1 = 8; 
    N_step = 8;
elseif  sigma_est_b >= 30 && sigma_est_b < 40    
    N1 = 8;
    N_step = 8;
elseif  sigma_est_b >= 40 && sigma_est_b < 50
    N1 = 9;
    N_step = 9;
else
    N1 = 9;
    N_step = 9;
end


        [imr1,imr2,imr3] = NLH_Basic(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(beta*imr2+(1-beta)*imp2),double(beta*imr3+(1-beta)*imp3),single(0.75),sigma_est_b/255,double(im_n(:,:,3)));
        
        N_step = N_step -1;
        
        [imr1,imr2,imr3] = NLH_Basic(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(beta*imr2+(1-beta)*imp2),double(beta*imr3+(1-beta)*imp3),single(0.75),sigma_est_b/255,double(im_n(:,:,1)));
        
        N_step = 3;
        
        [imr1,imr2,imr3] = NLH_Basic(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(beta*imr2+(1-beta)*imp2),double(beta*imr3+(1-beta)*imp3),single(1.5),sigma_est_b/255,double(im_n(:,:,2)));
        
      

imr=zeros(M1,M2,M3);
imr(:,:,1)=imr1(:,:);
imr(:,:,2)=imr2(:,:);
imr(:,:,3)=imr3(:,:);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert back to RGB colorspace
%%%
imr_basic = function_LumChrom2rgb(imr, colorspace);

imr_basic = double(imr_basic);


if sigma_est1 > 0.0

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% NLH_wiener
N2 = 64;
N3 = 8;
Ns = 129;



N1     = 16;
N_step = 16;



imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
 

  

beta = 0.8;


if sigma_est_b < 10
    iter = 1;
elseif sigma_est_b >= 10 && sigma_est_b < 20
    iter = 1;
 elseif sigma_est_b >= 20 && sigma_est_b < 30
    iter = 1;   
  elseif sigma_est_b >= 30 && sigma_est_b < 40
    iter = 1; 
    elseif sigma_est_b >= 40 && sigma_est_b < 50
    iter = 2;  
else 
    iter = 2;
end


gamma = 5.0 + iter + 10.0*sigma_est1/(sigma_est_b+0.00000001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic,

    [y_est1,y_est2,y_est3] = NLH_Wiener(8,N2,N3,Ns,8,double(imp1),double(imp2),double(imp3),single(gamma),...
                sigma_est_b/255, double(dis_map*255),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(imr2*beta+imp1*(1-beta)));            
     for ii = 1:iter
    [y_est1,y_est2,y_est3] = NLH_Wiener(8,N2,N3,Ns,8,double(y_est1),double(y_est2),double(y_est3),single(gamma),...
                sigma_est_b/255, double(dis_map*255),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(imr1*beta+imp1*(1-beta)));
     end
           
      [y_est1,y_est2,y_est3] = NLH_Wiener(N1,N2,N3,59,N_step,double(y_est1),double(y_est2),double(y_est3),single(gamma),...
                sigma_est_b/255, double(dis_map*255),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(y_est1));
     
     
         
% toc,

y_est = zeros(M1,M2,M3);
y_est(:,:,1)=y_est1(:,:);
y_est(:,:,2)=y_est2(:,:);
y_est(:,:,3)=y_est3(:,:);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert back to RGB colorspace
%%%%
imr_final = function_LumChrom2rgb(y_est, colorspace);

else
    
imr_final = imr_basic;

end

