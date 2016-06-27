% ==== (SNS matlab code)======
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

clear all; clc

im =  imread('tajmahal.png');
im_SNS = imread('tajmahal_0.50_sns.png');

% im =  imread('Brasserie_L_Aficion.png');
% im_SNS = imread('Brasserie_L_Aficion_0.50_sns.png');

% parameters
Ratio = 0.5;
mesh_size = 20; % using mesh_size x mesh_size quad
[h, w, ~] = size(im);
quad_num_h = floor(h/mesh_size);
quad_num_w = floor(w/mesh_size);

% the regular mesh on original image
Vertex_set_org = ImgRegualrMeshGrid(im, mesh_size);

% the importance map generation
importance_map =  SNS_importanceMap(im, true); % generate the importance map
importance_quad = SNS_importanceMap_quad(importance_map, Vertex_set_org);
importance_quad = importance_quad/sum(importance_quad(:)); % the importance weight for the quad

% the naive initialization of the mesh
% retargeting on the width
Vertex_warped_initial = Vertex_set_org;
Vertex_warped_initial(:,:,2) = Vertex_warped_initial(:,:,2)*Ratio;

% the mesh grid optimization
[Vertex_updated] = ...
    SNS_optimization(Vertex_set_org ,Vertex_warped_initial, importance_quad);

% warp the new image
im_warped = MeshBasedImageWarp(im, [1 Ratio], Vertex_set_org, Vertex_updated);
figure; subplot(1,2,1); imshow(im_warped); title(['My warped'], 'FontSize' , 15); 
subplot(1,2,2); imshow(im_SNS); title(['Original SNS warped'], 'FontSize' , 15); 

% show the mesh grid on the original image and retargeted image
% MeshGridImgPlot(im, Vertex_set_org, [0.5 0.0 0.5]);
% title(['Regular mesh grid on original image'], 'FontSize' , 15);
% MeshGridImgPlot(im_warped, Vertex_updated, [0.5 0.0 0.5]);
% title(['Warped image '], 'FontSize' , 15); 




