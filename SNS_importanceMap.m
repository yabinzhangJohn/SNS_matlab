function [ importance_map ] = SNS_importanceMap(im, SHOW_MAP)
% Summary of this function goes here
%   Detailed explanation goes here
% The Code (Version 1) is created by ZHANG Yabin,
% Nanyang Technological University, 2015-12-30
% which is based on the method described in the following paper 
% [1] Wang, Yu-Shuen, et al. "Optimized scale-and-stretch for image resizing." 
% ACM Transactions on Graphics (TOG) 27.5 (2008): 118. 
% The binary code is provided on the project page:
% http://graphics.csie.ncku.edu.tw/Image_Resizing/
% The Matlab codes are for non-comercial use only.
% Note that the importance maps are slightly different from the original
% ones, and the retargeted images are influenced.


addpath(genpath('.\gbvs\'));
out_itti = ittikochmap(im);
mapbig = out_itti.master_map_resized; 



L = 0.06 * double(im(:,:,1)) + 0.63 * double(im(:,:,2)) + 0.27 * double(im(:,:,3));
dx = [3 0 -3; 10 0 -10;  3  0 -3]/16;
dy = [3 10 3; 0  0   0; -3 -10 -3]/16;

IxL1 = conv2(L, dx, 'same');     
IyL1 = conv2(L, dy, 'same');    
GM = sqrt(IxL1.^2 + IyL1.^2);
GM = GM/max(GM(:));

%     [GM, ~] = imgradient(L,'Sobel');
%     GM = GM/max(GM(:));

    
importance_map = mapbig.*GM;
% importance_map = GM;
% importance_map = mapbig;

if(SHOW_MAP)

    colormap_sp = ...
    [254    11     1
   255    26     0
   254    36     0
   255    47     0
   254    60     0
   254    71     0
   255    83     1
   255    96     0
   254   106     0
   255   118     1
   254   131     0
   255   142     0
   255   153     2
   254   166     0
   254   177     0
   255   188     0
   254   201     0
   255   212     0
   255   223     1
   254   237     0
   254   247     0
   249   255     0
   236   255     0
   225   255     0
   214   255     0
   201   255     0
   190   255     0
   180   255     0
   167   255     0
   155   255     1
   145   255     0
   131   255     0
   120   255     0
   110   254     0
    96   255     0
    85   255     0
    75   255     0
    61   255     0
    50   255     0
    40   254     0
    26   255     0
    16   255     0
     4   255     0
     0   247     7
     0   236    16
     0   224    27
     1   210    45
     0   199    55
     0   188    64
     0   173    79
     0   162    91
     0   151   100
     1   139   115
     0   128   125
     0   116   137
     1   102   154
     0    91   163
     0    81   173
     0    66   189
     1    55   199
     0    43   211
     0    29   225
     0    19   236
     0     6   249];

    colormap_sp = double(colormap_sp)/255;
    colormap_sp = flipud(colormap_sp);
    figure('units' ,'normalized' ,'outerposition' ,[0.01 0.01 0.95 0.95]); 
    imagesc(importance_map); axis equal; axis off; colormap(colormap_sp); % colorbar
end

end

