function y=function_LumChrom2rgb(x,colormode)
% Inverse color-space transformation   ( forward transformation is function_rgb2LumChrom.m )
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2016   Public release v1.43 (15 May 2016)
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
% SYNTAX:
%
%   x = function_LumChrom2rgb(y,colormode);
%
% INPUTS:
%  y  is color-transformed image (with range typically included in or equal to [0 1]^3, depending on the transformation matrix)
%
%  colormode = 'dct', 'opp', 'yuv', or a custom 3x3 matrix (e.g. provided by the forward transform when 'pca' is selected)
%
%       'dct' or 'opp'  Opponent color space ('dct' is orthonormal version, 'opp' is equirange version)
%       'yuv'           standard YUV (e.g. for JPEG images)
%
% OUTPUTS:
%   x  is RGB image (with range [0 1]^3)
%
%
% NOTE:    'dct' is used by default if no colormode is specified
%

if nargin==1
    colormode='dct';
end
if size(colormode)==[3 3]
    A=colormode;
    B=inv(A);
else
    if strcmp(colormode,'dct')
        %  A=dct(eye(3))/sqrt(3);
        %  B=inv(A);
        A=[1/3 1/3 1/3;1/sqrt(6) 0 -1/sqrt(6);1/sqrt(18) -sqrt(2)/3 1/sqrt(18)];
        B=[1 sqrt(3/2) 1/sqrt(2);1 0 -sqrt(2);1 -sqrt(3/2) 1/sqrt(2)];
    end
    if strcmp(colormode,'opp')
        %  A=dct(eye(3));
        %  A=A./repmat(sum(A.*(A>0),2)-sum(A.*(A<0),2),[1 3]);   %% ranges are normalized to unitary length;
        %  B=inv(A);
        A=[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
        B=[1 1 2/3;1 0 -4/3;1 -1 2/3];
    end
    if strcmp(colormode,'yuv')
        A=[0.29900000000000   0.58700000000000   0.11400000000000;   -0.16873660714285  -0.33126339285715   0.50000000000000;   0.50000000000000  -0.41868750000000  -0.08131250000000];
        B=inv(A);
    end
end

x1=x(:,:,1)-(1-sum(A(1,:).*(A(1,:)>0),2)-sum(A(1,:).*(A(1,:)<0),2))/2;
x2=x(:,:,2)-(1-sum(A(2,:).*(A(2,:)>0),2)-sum(A(2,:).*(A(2,:)<0),2))/2;
x3=x(:,:,3)-(1-sum(A(3,:).*(A(3,:)>0),2)-sum(A(3,:).*(A(3,:)<0),2))/2;
y(:,:,1)=x1*B(1,1) + x2*B(1,2) + x3*B(1,3);
y(:,:,2)=x1*B(2,1) + x2*B(2,2) + x3*B(2,3);
y(:,:,3)=x1*B(3,1) + x2*B(3,2) + x3*B(3,3);
