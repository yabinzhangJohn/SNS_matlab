function [ Vertex_updated] = SNS_Opt_Iter_M(Vertex_set_org ,Vertex_warped, importance_quad)
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

[h_V, w_V, ~] = size(Vertex_set_org);

% compute the sysmetric matrix
% sf update....
Quad_edge =[0 0 1 0; 0 0 0 1; 0 1 1 1; 1 0 1 1];
Sf_matrix = zeros(h_V-1, w_V-1); % quad based, each quad one sf;
for Quad_h =  1:h_V-1
    for Quad_w = 1:w_V-1 
        for i = 1:4
            v_start = [Vertex_set_org(Quad_h+Quad_edge(1,i), Quad_w+Quad_edge(2,i), 1) ...
                Vertex_set_org(Quad_h+Quad_edge(1,i), Quad_w+Quad_edge(2,i), 2)];
            v_end = [Vertex_set_org(Quad_h+Quad_edge(3,i), Quad_w+Quad_edge(4,i), 1) ...
                Vertex_set_org(Quad_h+Quad_edge(3,i), Quad_w+Quad_edge(4,i), 2)];
            vw_start = [Vertex_warped(Quad_h+Quad_edge(1,i), Quad_w+Quad_edge(2,i), 1) ...
                Vertex_warped(Quad_h+Quad_edge(1,i), Quad_w+Quad_edge(2,i), 2)];
            vw_end = [Vertex_warped(Quad_h+Quad_edge(3,i), Quad_w+Quad_edge(4,i), 1) ...
                Vertex_warped(Quad_h+Quad_edge(3,i), Quad_w+Quad_edge(4,i), 2)];
            upper_term(i) = (vw_start- vw_end)*(v_start- v_end)';
            bottom_term(i) = norm(v_start- v_end)^2;
        end
        Sf_matrix(Quad_h, Quad_w) = sum(upper_term)/sum(bottom_term);
    end
end

% grid line bending term
L_edge_hor = zeros(h_V, w_V-1);
L_edge_ver = zeros(h_V-1, w_V);
% the horizontal edges
for i = 1:h_V
    for j = 1:w_V-1
        v_start = [Vertex_set_org(i, j, 1) Vertex_set_org(i, j, 2)];
        v_end = [Vertex_set_org(i, j+1, 1) Vertex_set_org(i, j+1, 2)];
        vw_start = [Vertex_warped(i, j, 1) Vertex_warped(i, j, 2)];
        vw_end = [Vertex_warped(i, j+1, 1) Vertex_warped(i, j+1, 2)];
        foo_l = norm(vw_start - vw_end)/norm(v_start - v_end);
        L_edge_hor(i, j) = foo_l;
    end
end
% the vertical edges
for i = 1:h_V-1
    for j = 1:w_V
        v_start = [Vertex_set_org(i, j, 1) Vertex_set_org(i, j, 2)];
        v_end = [Vertex_set_org(i+1, j, 1) Vertex_set_org(i+1, j, 2)];
        vw_start = [Vertex_warped(i, j, 1) Vertex_warped(i, j, 2)];
        vw_end = [Vertex_warped(i+1, j, 1) Vertex_warped(i+1, j, 2)];
        foo_l = norm(vw_start - vw_end)/norm(v_start - v_end);
        L_edge_ver(i, j) = foo_l;
    end
end


h_Boundary = Vertex_warped(h_V, w_V, 1);
w_Boundary = Vertex_warped(h_V, w_V, 2);

BENDING_TERM = 1;
Lambda_term = 0.01;

Vertex_updated = Vertex_warped;
% construct the system matrix
% layer 1 --- H (Y) % layer 2 --- W (X)
for Layer_Vertex = [2 1] % for H and W layer
    A_matrix = zeros(h_V*w_V, h_V*w_V);
    B_vector = zeros(h_V*w_V, 1);
    Vect_vertex_warped_old = reshape(Vertex_warped(:,:,Layer_Vertex), [h_V*w_V, 1]);
    for Q_h =  1:h_V
        for Q_w = 1:w_V

            Vector_loc = Q_h + (Q_w-1)*h_V;
            
            % ##################################
            % ##### add the quad deformation part coefficients
            % 1 === the top-left quad
            if( (Q_h - 1) > 0 && (Q_w - 1) > 0)
                wf_quad = importance_quad(Q_h - 1, Q_w - 1);
                sf_quad = Sf_matrix(Q_h - 1, Q_w - 1);

                A_matrix(Vector_loc, Vector_loc) = ...
                    A_matrix(Vector_loc, Vector_loc) + 2*wf_quad;
                % T
                A_matrix(Vector_loc, Vector_loc-1) = ...
                    A_matrix(Vector_loc, Vector_loc-1) - wf_quad;
                % L
                A_matrix(Vector_loc, Vector_loc-h_V) = ...
                    A_matrix(Vector_loc, Vector_loc-h_V) - wf_quad;
                % B_vector
                B_vector(Vector_loc) = ...
                    B_vector(Vector_loc) + wf_quad*sf_quad*(2*Vertex_set_org(Q_h,Q_w,Layer_Vertex) ...
                    - Vertex_set_org(Q_h-1,Q_w,Layer_Vertex) - Vertex_set_org(Q_h,Q_w-1,Layer_Vertex));
            end

            % 2 === the top-right quad
            if( (Q_h - 1) > 0 && (Q_w + 1) <= w_V)
                wf_quad = importance_quad(Q_h - 1, Q_w);
                sf_quad = Sf_matrix(Q_h - 1, Q_w);

                A_matrix(Vector_loc, Vector_loc) = ...
                    A_matrix(Vector_loc, Vector_loc) + 2*wf_quad;
                % T
                A_matrix(Vector_loc, Vector_loc-1) = ...
                    A_matrix(Vector_loc, Vector_loc-1) - wf_quad;
                % R
                A_matrix(Vector_loc, Vector_loc+h_V) = ...
                    A_matrix(Vector_loc, Vector_loc+h_V) - wf_quad;
                % B_vector
                B_vector(Vector_loc) = ...
                    B_vector(Vector_loc) + wf_quad*sf_quad*(2*Vertex_set_org(Q_h,Q_w,Layer_Vertex) ...
                    - Vertex_set_org(Q_h-1,Q_w,Layer_Vertex) - Vertex_set_org(Q_h,Q_w+1,Layer_Vertex));

            end

            % 3 === the bottom-left quad
            if( (Q_h + 1) <= h_V && (Q_w - 1) > 0)
                wf_quad = importance_quad(Q_h, Q_w - 1);
                sf_quad = Sf_matrix(Q_h, Q_w - 1);

                A_matrix(Vector_loc, Vector_loc) = ...
                    A_matrix(Vector_loc, Vector_loc) + 2*wf_quad;
                % D
                A_matrix(Vector_loc, Vector_loc+1) = ...
                    A_matrix(Vector_loc, Vector_loc+1) - wf_quad;
                % L
                A_matrix(Vector_loc, Vector_loc-h_V) = ...
                    A_matrix(Vector_loc, Vector_loc-h_V) - wf_quad;
                % B_vector
                B_vector(Vector_loc) = ...
                    B_vector(Vector_loc) + wf_quad*sf_quad*(2*Vertex_set_org(Q_h,Q_w,Layer_Vertex) ...
                    - Vertex_set_org(Q_h+1,Q_w,Layer_Vertex) - Vertex_set_org(Q_h,Q_w-1,Layer_Vertex));    
            end

            % 4 === the bottom-right quad
            if( (Q_h + 1) <= h_V && (Q_w + 1) <= w_V)
                wf_quad = importance_quad(Q_h, Q_w);
                sf_quad = Sf_matrix(Q_h, Q_w);

                A_matrix(Vector_loc, Vector_loc) = ...
                    A_matrix(Vector_loc, Vector_loc) + 2*wf_quad;
                % D
                A_matrix(Vector_loc, Vector_loc+1) = ...
                    A_matrix(Vector_loc, Vector_loc+1) - wf_quad;
                % R
                A_matrix(Vector_loc, Vector_loc+h_V) = ...
                    A_matrix(Vector_loc, Vector_loc+h_V) - wf_quad;
                % B_vector
                B_vector(Vector_loc) = ...
                    B_vector(Vector_loc) + wf_quad*sf_quad*(2*Vertex_set_org(Q_h,Q_w,Layer_Vertex) ...
                    - Vertex_set_org(Q_h+1,Q_w,Layer_Vertex) - Vertex_set_org(Q_h,Q_w+1,Layer_Vertex));
            end
            
            % ##################################
            if(BENDING_TERM)
                % ##### add the grid line bending part coefficients
                % 1 === the left edge part (horizontal) 
                if( Q_w > 1)
                    L_edge = L_edge_hor(Q_h, Q_w-1);
                    Nb_H = 0; Nb_W = -1;

                    A_matrix(Vector_loc, Vector_loc) = ...
                        A_matrix(Vector_loc, Vector_loc) + 1*Lambda_term;
                    % the left vertex
                    A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) = ...
                        A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) - 1*Lambda_term;
                    % B_vector
                    B_vector(Vector_loc) = B_vector(Vector_loc) + ...
                        L_edge*(Vertex_set_org(Q_h,Q_w,Layer_Vertex) - ...
                        Vertex_set_org(Q_h+Nb_H, Q_w+Nb_W,Layer_Vertex))*Lambda_term;
                end
                % 2 === the top edge part (vertical) 
                if( Q_h > 1)
                    L_edge = L_edge_ver(Q_h-1, Q_w);
                    Nb_H = -1; Nb_W = 0;

                    A_matrix(Vector_loc, Vector_loc) = ...
                        A_matrix(Vector_loc, Vector_loc) + 1*Lambda_term;
                    % the top vertex
                    A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) = ...
                        A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) - 1*Lambda_term;
                    % B_vector
                    B_vector(Vector_loc) = B_vector(Vector_loc) + ...
                        L_edge*(Vertex_set_org(Q_h,Q_w,Layer_Vertex) - ...
                        Vertex_set_org(Q_h+Nb_H, Q_w+Nb_W,Layer_Vertex) )*Lambda_term;
                end
                % 3 === the right edge part (horizontal) 
                if( Q_w < w_V)
                    L_edge = L_edge_hor(Q_h, Q_w);
                    Nb_H = 0; Nb_W = 1;

                    A_matrix(Vector_loc, Vector_loc) = ...
                        A_matrix(Vector_loc, Vector_loc) + 1*Lambda_term;
                    % the right vertex
                    A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) = ...
                        A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) - 1*Lambda_term;
                    % B_vector
                    B_vector(Vector_loc) = B_vector(Vector_loc) + ...
                        L_edge*(Vertex_set_org(Q_h,Q_w,Layer_Vertex) - ...
                        Vertex_set_org(Q_h+Nb_H, Q_w+Nb_W,Layer_Vertex) )*Lambda_term;
                end
                % 4 === the bottom edge part (vertical) 
                if( Q_h < h_V)
                    L_edge = L_edge_ver(Q_h, Q_w);
                    Nb_H = 1; Nb_W = 0;

                    A_matrix(Vector_loc, Vector_loc) = ...
                        A_matrix(Vector_loc, Vector_loc) + 1*Lambda_term;
                    % the top vertex
                    A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) = ...
                        A_matrix(Vector_loc+Nb_H, Vector_loc+Nb_W*h_V) - 1*Lambda_term;
                    % B_vector
                    B_vector(Vector_loc) = B_vector(Vector_loc) + ...
                        L_edge*(Vertex_set_org(Q_h,Q_w,Layer_Vertex) - ...
                        Vertex_set_org(Q_h+Nb_H, Q_w+Nb_W,Layer_Vertex) )*Lambda_term;
                end
                
            end
        end

    end
    
    N = 1;
    %These constraints are simply substituted into the linear system during the optimization
    START_POINT = 0;
    if(Layer_Vertex == 1)
        for Q_h =  1
            for Q_w = 1:w_V
                Vector_loc = Q_h + (Q_w-1)*h_V;
                A_matrix(Vector_loc, :) = 0;
                A_matrix(Vector_loc, Vector_loc) = 1*N;
                B_vector(Vector_loc) = START_POINT*N;
            end
        end
        for Q_h =  h_V
            for Q_w = 1:w_V
                Vector_loc = Q_h + (Q_w-1)*h_V;
                A_matrix(Vector_loc, :) = 0;
                A_matrix(Vector_loc, Vector_loc) = 1;
                B_vector(Vector_loc) = h_Boundary;
            end
        end
    else
        for Q_h =  1:h_V
            for Q_w = 1
                Vector_loc = Q_h + (Q_w-1)*h_V;
                A_matrix(Vector_loc, :) = 0;
                A_matrix(Vector_loc, Vector_loc) = 1*N;
                B_vector(Vector_loc) = START_POINT*N;
            end
        end
        for Q_h = 1:h_V
            for Q_w = w_V
                Vector_loc = Q_h + (Q_w-1)*h_V;
                A_matrix(Vector_loc, :) = 0;
                A_matrix(Vector_loc, Vector_loc) = 1;
                B_vector(Vector_loc) = w_Boundary;
            end
        end
    end
    
    
    
    [L, U] = lu(A_matrix);
    Vect_vertex_warped_factorization = U\(L\B_vector);

    Vect_vertex_warped_factorization = 0.7*Vect_vertex_warped_factorization + ...
        0.3*Vect_vertex_warped_old;
    foo = reshape(Vect_vertex_warped_factorization, [h_V w_V]);
    Vertex_updated(:, :, Layer_Vertex) = foo;
end






end

