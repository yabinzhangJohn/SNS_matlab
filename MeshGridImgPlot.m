function [  ] = MeshGridImgPlot(im, Vertex_set, COLOR_V)
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

if(nargin < 3)
    COLOR_V = [0.0 0.5 0.0];
    % COLOR_V = [0.5 0.0 0.5];
end

tic
figure('units' ,'normalized' ,'outerposition' ,[0.01 0.01 0.95 0.95]);
imshow(im); hold on;

Vertex_set( Vertex_set == 0) = 1;

[ver_h, ver_w, ~] = size(Vertex_set);
quad_num_h = ver_h - 2;
quad_num_w = ver_w - 2;

for vect_i = 1:quad_num_h+1
    for vect_j = 1:quad_num_w+1
        pos_ini_x = Vertex_set(vect_i, vect_j, 2);
        pos_ini_y = Vertex_set(vect_i, vect_j, 1);
        pos_end_x = Vertex_set(vect_i, vect_j+1, 2);
        pos_end_y = Vertex_set(vect_i, vect_j+1, 1);
        line([pos_ini_x pos_end_x],...
                [pos_ini_y pos_end_y],...
                'Color', COLOR_V, 'LineWidth',2);
        pos_end_x = Vertex_set(vect_i+1, vect_j, 2);
        pos_end_y = Vertex_set(vect_i+1, vect_j, 1);
        line([pos_ini_x pos_end_x],...
                [pos_ini_y pos_end_y],...
                'Color', COLOR_V, 'LineWidth',2);
            
    end
end
for vect_i = quad_num_h+2    
    for vect_j = 1:quad_num_w+1
        pos_ini_x = Vertex_set(vect_i, vect_j, 2);
        pos_ini_y = Vertex_set(vect_i, vect_j, 1);
        pos_end_x = Vertex_set(vect_i, vect_j+1, 2);
        pos_end_y = Vertex_set(vect_i, vect_j+1, 1);
        line([pos_ini_x pos_end_x],...
                [pos_ini_y pos_end_y],...
                'Color', COLOR_V, 'LineWidth',2);
    end
end
for vect_i = 1:quad_num_h+1
    for vect_j = quad_num_w+2
        pos_ini_x = Vertex_set(vect_i, vect_j, 2);
        pos_ini_y = Vertex_set(vect_i, vect_j, 1);
        pos_end_x = Vertex_set(vect_i+1, vect_j, 2);
        pos_end_y = Vertex_set(vect_i+1, vect_j, 1);
        line([pos_ini_x pos_end_x],...
                [pos_ini_y pos_end_y],...
                'Color', COLOR_V, 'LineWidth',2);
    end
end

mesh_grid_time = toc;
disp(['Mesh grid plot time consuming: ' num2str(mesh_grid_time) '.']);

end

