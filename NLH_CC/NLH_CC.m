function imr_final = NLH_CC(im_n)

    
   [M1,M2,M3]=size(im_n);

   
   
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


[sigma_tb2 sigma_est dis_map2] = Noise_estimation(im_n(:,:,2)); %% estimate noisy level
sigma_est = double(sigma_est);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dis_map,sigma_tb1] = NL_distance(8,16,2,39,single(8),double(im_n(:,:,2)));
dis_map = double(dis_map);

%%%%%%%%%%%%%%%%%% 
Ns     = 43;
N3     = 4;
N_step = 6;
N1     = 6;
N2     = 16;
%%%%%%%%%%%%%%%%





alpha = 0.618; 

sigma_est_b = sigma_est*4.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   µü´ú´ÎÊýÓëãÐÖµÉèÖÃ


    k   = 2;
    
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
for ll = 1:k 


        [imr1,imr2,imr3] = NLH_Color(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(alpha*imr2+(1-alpha)*imp2),double(alpha*imr3+(1-alpha)*imp3),single(Thr),...
        sigma_est_b/255,double(im_n(:,:,2)));
        N_step = 3;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation

imr = rgb2gray(imr_basic);

[sigma_tb11 sigma_est1 dis_map11] = Noise_estimation(imr); %% estimate noisy level

sigma_est1 = double(sigma_est1);


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
 

gamma = 2.0;  
 
beta = 0.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       [y_est1,y_est2,y_est3] = NLH_Wiener_Color(N1,N2,N3,Ns,N_step,double(imp1),double(imp2),double(imp3),single(gamma),...
                sigma_est_b/255, double(dis_map),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(imr2*beta+imp1*(1-beta)));            

                        
        [y_est1,y_est2,y_est3] = NLH_Wiener_Color(N1,N2,N3,Ns,N_step,double(y_est1),double(y_est2),double(y_est3),single(gamma),...
                sigma_est_b/255, double(dis_map),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(y_est1*beta+imp1*(1-beta)));
         
     
            
y_est = zeros(M1,M2,M3);
y_est(:,:,1)=y_est1(:,:);
y_est(:,:,2)=y_est2(:,:);
y_est(:,:,3)=y_est3(:,:);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert back to RGB colorspace

 imr_final = function_LumChrom2rgb(y_est, colorspace);


 
 


 
 







