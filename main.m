close all; clear all; clc

%% === Load in (1) coordinates of the local nodes (2) element connectivity ===
globalnode_loc_all = dlmread('globalnode_all.dat');
localnode_order_all = dlmread('localnode_order_all.dat');
for ielement = 1:5
    node_order = localnode_order_all(ielement, :);
    coord_element = globalnode_loc_all(node_order, :);
    coord_all{ielement} = coord_element;
end

%% === Compute local stiffness matrix [k] for all elements ===
E = 1e9*1e6; % in N/(mm^2)  Note: Pa = N/(m^2),
mu = 0.3; % Poisson's ratio
Wi = 1; % Weight
thickness = 1; % in mm

for ielement = 1:5
    % --- Compute local [k] for each element ---
    coord = coord_all{ielement};
    int_points_all = [-sqrt(1/3), -sqrt(1/3);
        -sqrt(1/3), sqrt(1/3);
        sqrt(1/3), sqrt(1/3);
        sqrt(1/3), -sqrt(1/3)];
    for i = 1:4
        int_points = int_points_all(i,:);
        [ki, ~, ~] = get_ki( coord, int_points, E, mu, Wi);
        ki_all(:,:,i) = ki;
    end
    k_loc = thickness*sum(ki_all, 3);
    k_loc_all(:,:,ielement) = k_loc;
end

%% === Assembly to global stiffness matrix [K] ======
Num_global_node = size(globalnode_loc_all,1);

for ielement = 1:5
    k_loc = k_loc_all(:,:,ielement);
    node_global = localnode_order_all(ielement, :);
    
    k_loc2global = zeros(Num_global_node*2, Num_global_node*2);
    for i = 1:4
        k_loc2global(node_global(i)*2 - 1: node_global(i)*2, node_global(i)*2 - 1: node_global(i)*2 ) =  k_loc(i*2-1:i*2, i*2-1:i*2);
    end
    k_loc2global_all(:,:,ielement) = k_loc2global;
    
end
K_global = sum(k_loc2global_all, 3);

%% === Get local force vectors {r} ======
% qy imposed on element B & D, node 2, 3
% both elements have f_2y = f_3y = (1/4)*(1+eta)*qy*dxhatdzeta
qy = 3000; % in N/mm
qx = 1000; % in N/mm
eta = 1;
zeta = 1;

element_q1 = [2; 4]; %element B & D 
for ielement = 1:2
    coord_element = coord_all{element_q1(ielement)}  ;
    x1 = coord_element(1, 1); y1 = coord_element(1, 2);
    x2 = coord_element(2, 1); y2 = coord_element(2, 2);
    x3 = coord_element(3, 1); y3 = coord_element(3, 2);
    x4 = coord_element(4, 1); y4 = coord_element(4, 2);
    
    dxhatdzeta = 0.25*(-x1*(1-eta) - x2*(1+eta) + x3*(1+eta) + x4*(1-eta));
    r_q1_y(ielement, 1) = (1/4)*(1+eta)*qy*dxhatdzeta;
end
r_qy2_eleB = r_q1_y(1, 1);
r_qy3_eleB = r_q1_y(1, 1);
r_qy2_eleD = r_q1_y(2, 1);
r_qy3_eleD = r_q1_y(2, 1);

% qx imposed on element E, node 3, 4
% f_3x = (3/2)*C
% f_4x = (1/2)*C
% C = (1/4)*(1+zeta)*qx*dyhatdeta
coord_element = coord_all{5};
x1 = coord_element(1, 1); y1 = coord_element(1, 2);
x2 = coord_element(2, 1); y2 = coord_element(2, 2);
x3 = coord_element(3, 1); y3 = coord_element(3, 2);
x4 = coord_element(4, 1); y4 = coord_element(4, 2);

dyhatdeta = 0.25*(-y1*(1-zeta) + y2*(1-zeta) + y3*(1+zeta) - y4*(1+zeta));
C = (1/4)*(1+zeta)*qx*dyhatdeta;
r_qx3_eleE =  (3/2)*C;
r_qx4_eleE =  (1/2)*C;

% --- assign to local elements ---
r_qy_eleB = zeros(8, 1); r_qy_eleD = zeros(8, 1); r_qx_eleE = zeros(8, 1);
r_qy_eleB(4, 1) = r_qy2_eleB; r_qy_eleB(6, 1) = r_qy3_eleB;
r_qy_eleD(4, 1) = r_qy2_eleD; r_qy_eleD(6, 1) = r_qy3_eleD;
r_qx_eleE(5, 1) = r_qx3_eleE; r_qx_eleE(7, 1) = r_qx4_eleE;

%% === Assembly to global {r} ====
node_global_eleB = localnode_order_all(2, :);
node_global_eleD = localnode_order_all(4, :);
node_global_eleE = localnode_order_all(5, :);

r_eleB_loc2global = zeros(Num_global_node*2, 1);
r_eleD_loc2global = zeros(Num_global_node*2, 1);
r_eleE_loc2global = zeros(Num_global_node*2, 1);
for i = 1:4
    r_eleB_loc2global(node_global_eleB(i)*2, 1 ) =  r_qy_eleB(i*2, 1);
    r_eleD_loc2global(node_global_eleD(i)*2, 1 ) =  r_qy_eleD(i*2, 1);
    r_eleE_loc2global(node_global_eleE(i)*2 - 1, 1 ) =  r_qx_eleE(i*2 - 1, 1);
end
R_global = r_eleB_loc2global + r_eleD_loc2global + r_eleE_loc2global;

%% === Solve for nodal displacement ===
% global [K]{U} = {R}
U_global = K_global\R_global;

% ---- apply boundary condition ---
% u1, v1, u2, u3, v4, v7, v9 = 0
U_global(1:3, 1) = 0;
U_global(5, 1) = 0;
U_global(8, 1) = 0;
U_global(14, 1) = 0;
U_global(18, 1) = 0;

% -- nodal displacement ---
U_global_normalized = U_global./max(abs(U_global));

%% === Plot + display ===
U_global_normalized_reshape = zeros(10, 2);
for i = 1:10
    U_global_normalized_reshape(i, 1) = U_global_normalized(i*2-1,1);  %x
    U_global_normalized_reshape(i, 2) = U_global_normalized(i*2,1);  %y
end

for ielement = 1:5
    node_order = [localnode_order_all(ielement, :), localnode_order_all(ielement, 1)];
    coord_element = globalnode_loc_all(node_order, :);
    u_element = U_global_normalized_reshape(node_order, :);
    coord_original{ielement} = coord_element;
    coord_deformed{ielement} = coord_element + u_element;
end

figure(1); clf
for ielement = 1:5
    coord_element = coord_original{ielement};
    plot(coord_element(:,1), coord_element(:,2), 'k', 'Linewidth', 3); hold on;
end
hold on
for ielement = 1:5
    coord_deformed_ele = coord_deformed{ielement};
    plot(coord_deformed_ele(:,1), coord_deformed_ele(:,2), 'r--', 'Linewidth', 2.5); hold on;
end
hold off
xlim([-0.5, 12]); ylim([-0.5, 6.5])
xlabel('x_1 (mm)'); ylabel('x_2 (mm)')
title('Original (black) vs Deformed (red) structure')
box off

%% === Compute sigma_xx, sigma_yy at element E (enclosed by node 7, 8, 10, 9) ===
node_order_E = localnode_order_all(5, :);
coord_E = coord_all{5};
for i = 1:4
    U_E(i*2-1:i*2, 1) = U_global(node_order_E(i)*2-1:node_order_E(i)*2);
end

for i = 1:4
    int_points = int_points_all(i,:);
    [~, Bi, D] = get_ki( coord_E, int_points, E, mu, Wi);
    strain = Bi*U_E;
    stress = D*strain;
    stress_all(:,i) = stress;
end

stress_xx_all = stress_all(1, :)
stress_yy_all = stress_all(2, :)
