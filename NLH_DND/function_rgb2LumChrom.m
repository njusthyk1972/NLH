function [y A l2normLumChrom]=function_rgb2LumChrom(x,colormode)
% Forward color-space transformation   ( inverse transformation is function_LumChrom2rgb.m )
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2016   Public release v1.43 (15 May 2016)
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
% SYNTAX:
%
%   [y A l2normLumChrom] = function_rgb2LumChrom(x,colormode);
%
% INPUTS:
%   x  is RGB image with range [0 1]^3
%
%   colormode = 'dct', 'opp', 'yuv', 'pca', or a custom 3x3 matrix
%
%       'dct' or 'opp'  Opponent color space ('dct' is orthonormal version, 'opp' is equirange version)
%       'yuv'           standard YUV (e.g. for JPEG images)
%       'pca'           Principal components   (note that this transformation is renormalized to be equirange) 
%
% OUTPUTS:
%   y  is color-transformed image (with range typically included in or equal to [0 1]^3, depending on the transformation matrix)
%
%   l2normLumChrom (optional) l2-norm of the transformation (useful for noise std calculation)
%   A  transformation matrix  (used necessarily if colormode='pca')
%
%   NOTES:  -  If only two outputs are used, then the second output is l2normLumChrom, unless colormode='pca';
%           -  'opp' is used by default if no colormode is specified.
%
%
% USAGE EXAMPLE FOR PCA TRANSFORMATION:
%  %%%%  -- forward color transformation --
%    if colormode=='pca'
%       [zLumChrom colormode]=function_rgb2LumChrom(zRGB,colormode);
%    else
%       zLumChrom=function_rgb2LumChrom(zRGB,colormode);
%    end
%
%  %%%% [ ... ]  Processing  [ ... ]
%
%  %%%%  -- inverse color transformation --
%    zRGB=function_LumChrom2rgb(zLumChrom,colormode);
%

if nargin==1
    colormode='opp';
end
change_output=0;
if size(colormode)==[3 3]
    A=colormode;
    l2normLumChrom=sqrt(sum(A.^2,2));
else
    if strcmp(colormode,'dct')
        %  A=dct(eye(3))/sqrt(3);
        A=[1/3 1/3 1/3;1/sqrt(6) 0 -1/sqrt(6);1/sqrt(18) -sqrt(2)/3 1/sqrt(18)];
    end
    if strcmp(colormode,'opp')
        %  A=dct(eye(3));
        %  A=A./repmat(sum(A.*(A>0),2)-sum(A.*(A<0),2),[1 3]);  %% ranges are normalized to unitary length;
        A=[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
    end
    if strcmp(colormode,'yuv')
        A=[0.29900000000000   0.58700000000000   0.11400000000000;   -0.16873660714285  -0.33126339285715   0.50000000000000;   0.50000000000000  -0.41868750000000  -0.08131250000000];
    end
    if strcmp(colormode,'pca')
        A=princomp(reshape(x,[size(x,1)*size(x,2) 3]))';
        A=A./repmat(sum(A.*(A>0),2)-sum(A.*(A<0),2),[1 3]);  %% ranges are normalized to unitary length;
    else
        if nargout==2
            change_output=1;
        end
    end
end
l2normLumChrom=sqrt(sum(A.^2,2)); %% l2-norms

%%% COLOR TRANSFORMATION %%%%%%%
y(:,:,1)=x(:,:,1)*A(1,1) + x(:,:,2)*A(1,2) + x(:,:,3)*A(1,3) + (1-sum(A(1,:).*(A(1,:)>0),2)-sum(A(1,:).*(A(1,:)<0),2))/2;
y(:,:,2)=x(:,:,1)*A(2,1) + x(:,:,2)*A(2,2) + x(:,:,3)*A(2,3) + (1-sum(A(2,:).*(A(2,:)>0),2)-sum(A(2,:).*(A(2,:)<0),2))/2;
y(:,:,3)=x(:,:,1)*A(3,1) + x(:,:,2)*A(3,2) + x(:,:,3)*A(3,3) + (1-sum(A(3,:).*(A(3,:)>0),2)-sum(A(3,:).*(A(3,:)<0),2))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if change_output
    A=l2normLumChrom;
end
