function [im_warped] = MeshBasedImageWarp(im, ratio ,Vertex_set_org, Vertex_set_warped)
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

% ratio is the target ratio [h_ratio, w_ratio]
tic

[h, w, ~] = size(im);
quad_num_h = size(Vertex_set_org,1) -1;
quad_num_w = size(Vertex_set_org,2) -1;
% GridWarped = NaN(round(h*ratio(1)), round(w*ratio(2)), 2);
GridWarped = zeros(round(h*ratio(1)), round(w*ratio(2)), 2);

for Quad_h =  1:quad_num_h
    for Quad_w = 1:quad_num_w
        V1 = [Vertex_set_org(Quad_h, Quad_w, 2) Vertex_set_org(Quad_h, Quad_w, 1)];
        V2 = [Vertex_set_org(Quad_h, Quad_w+1, 2) Vertex_set_org(Quad_h, Quad_w+1, 1)];
        V3 = [Vertex_set_org(Quad_h+1, Quad_w, 2) Vertex_set_org(Quad_h+1, Quad_w, 1)];
        V4 = [Vertex_set_org(Quad_h+1, Quad_w+1, 2) Vertex_set_org(Quad_h+1, Quad_w+1, 1)];
        
        V1w = [Vertex_set_warped(Quad_h, Quad_w, 2) Vertex_set_warped(Quad_h, Quad_w, 1)];
        V2w = [Vertex_set_warped(Quad_h, Quad_w+1, 2) Vertex_set_warped(Quad_h, Quad_w+1, 1)];
        V3w = [Vertex_set_warped(Quad_h+1, Quad_w, 2) Vertex_set_warped(Quad_h+1, Quad_w, 1)];
        V4w = [Vertex_set_warped(Quad_h+1, Quad_w+1, 2) Vertex_set_warped(Quad_h+1, Quad_w+1, 1)]; 
        
        Vmatrix = [V1 1; V2 1; V3 1; V4 1]';
        Vwmatrix = [V1w 1; V2w 1; V3w 1; V4w 1]';
        
        % compute the projective matrix for the entire quad face
        % svd based solution
        Amatrix = zeros(8,9);
        for i = 1:4
            x = Vmatrix(1,i); y = Vmatrix(2,i);
            xw = Vwmatrix(1,i); yw = Vwmatrix(2,i);
            Amatrix(2*i-1,:) = [-xw -yw -1 0 0 0 x*xw x*yw x];
            Amatrix(2*i,:) = [0 0 0 -xw -yw -1 y*xw y*yw y];
        end
        [~, ~, v] = svd(Amatrix); ColV = v(:,9);
        H = [ColV(1:3)'; ColV(4:6)'; ColV(7:9)']/ColV(9);
        
%         % matrix factorization based solution
%         Amatrix = zeros(8,9);
%         for i = 1:4
%             x = Vmatrix(1,i); y = Vmatrix(2,i);
%             xw = Vwmatrix(1,i); yw = Vwmatrix(2,i);
%             Amatrix(2*i-1,:) = [-xw -yw -1 0 0 0 x*xw x*yw x];
%             Amatrix(2*i,:) = [0 0 0 -xw -yw -1 y*xw y*yw y];
%         end
%         A_m = Amatrix(:,1:8); B_m = -1*Amatrix(:,9);
%         [L, U] = lu(A_m);
%         H_foo = U\(L\B_m); 
%         H = [H_foo(1:3)'; H_foo(4:6)'; H_foo(7:8)' 1];
        

        W_min = min([Vwmatrix(1,:)]); W_max = max([Vwmatrix(1,:)]);
        H_min = min([Vwmatrix(2,:)]); H_max = max([Vwmatrix(2,:)]);
        
        for grid_h = max(floor(H_min),1): ceil(H_max)
            for grid_w = max(floor(W_min),1): ceil(W_max)
                            
                if( (grid_w - V1w(1))*(V3w(2) - V1w(2)) - ...
                        (grid_h - V1w(2))*(V3w(1) - V1w(1)) >= 0 && ...
                    (grid_w - V2w(1))*(V4w(2) - V2w(2)) - ...
                        (grid_h - V2w(2))*(V4w(1) - V2w(1)) <= 0 &&...
                    (grid_w - V1w(1))*(V2w(2) - V1w(2)) - ...
                        (grid_h - V1w(2))*(V2w(1) - V1w(1)) <= 0 && ...
                     (grid_w - V3w(1))*(V4w(2) - V3w(2)) - ...
                        (grid_h - V3w(2))*(V4w(1) - V3w(1)) >= 0);
                    
                    grid_point = [grid_w; grid_h; 1];
                    grid_backward = H*grid_point;
                    grid_backward = grid_backward/grid_backward(3);
                    GridWarped(grid_h, grid_w, 1) = grid_backward(1);
                    GridWarped(grid_h, grid_w, 2) = grid_backward(2);
                end   
            end
        end 
        
    end
end

% deal with the out-of-border situations
GridWarped( GridWarped<1) = 1;
GridWarped_w = GridWarped(:,:,1); GridWarped_w(GridWarped_w > w) = w;
GridWarped(:,:,1) = GridWarped_w;
GridWarped_h = GridWarped(:,:,2); GridWarped_h(GridWarped_h > h) = h;
GridWarped(:,:,2) = GridWarped_h;

[X,Y] = meshgrid(1:w,1:h);
im_warped(:,:,1) = interp2(X,Y,double(im(:,:,1)),GridWarped(:,:,1),GridWarped(:,:,2));
im_warped(:,:,2) = interp2(X,Y,double(im(:,:,2)),GridWarped(:,:,1),GridWarped(:,:,2));
im_warped(:,:,3) = interp2(X,Y,double(im(:,:,3)),GridWarped(:,:,1),GridWarped(:,:,2));
im_warped = uint8(im_warped);
time_consumed = toc;
disp(['Warping time : ' num2str(time_consumed) ' s.']);

end

