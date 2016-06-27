function [ Vertex_set ] = ImgRegualrMeshGrid(im, mesh_size)
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

[h, w, ~] = size(im);
quad_num_h = floor(h/mesh_size);
quad_num_w = floor(w/mesh_size);
Remain_h = h - mesh_size*quad_num_h;
Remain_w = w - mesh_size*quad_num_w;

Start_point = 0;

Vertex_set = zeros(quad_num_h+1, quad_num_w+1,2);
for vect_i = 1:quad_num_h+1
    for vect_j = 1:quad_num_w+1
        if(vect_i ~= 1)
            Vertex_set(vect_i, vect_j, 1) = ...
                min(Start_point+mesh_size*(vect_i-1)+round((vect_i-1)/quad_num_h*Remain_h), h);
        else
            Vertex_set(vect_i, vect_j, 1) = Start_point;
        end
        if(vect_j ~= 1)
            Vertex_set(vect_i, vect_j, 2) = ...
                min(Start_point+mesh_size*(vect_j-1)+round((vect_j-1)/quad_num_w*Remain_w), w);             
        else
            Vertex_set(vect_i, vect_j, 2) = Start_point;
        end
    end
end

end

