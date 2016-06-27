function [ importance_quad] = SNS_importanceMap_quad(importance_map, Vertex_set_org)
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

    quad_num_h = size(Vertex_set_org,1) -1;
    quad_num_w = size(Vertex_set_org,2) -1;
    importance_quad = zeros(quad_num_h, quad_num_w);

    % layer 1 --- H (Y) % layer 2 --- W (X)
    for Quad_h =  1:quad_num_h
        for Quad_w = 1:quad_num_w
            V1 = [Vertex_set_org(Quad_h, Quad_w, 2) Vertex_set_org(Quad_h, Quad_w, 1)];
            V2 = [Vertex_set_org(Quad_h+1, Quad_w+1, 2) Vertex_set_org(Quad_h+1, Quad_w+1, 1)];
            foo_map = importance_map(max(V1(2),1):V2(2), max(V1(1),1):V2(1));

            importance_quad(Quad_h, Quad_w) = sum(sum(foo_map));
        end
    end


end

