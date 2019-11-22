function [ ki, B, D ] = get_ki( coord, int_points, E, mu, Wi )
% This function calculates location [k], [B], and elastic matrix [D]

% -- coordinate --
%coord = [x1, y1; x2, y2; x3, y3; x4, y4]
x1 = coord(1,1); y1 = coord(1,2);
x2 = coord(2,1); y2 = coord(2,2);
x3 = coord(3,1); y3 = coord(3,2);
x4 = coord(4,1); y4 = coord(4,2);

% -- [zeta, eta] --
zeta = int_points(1, 1); eta = int_points(1, 2);

% -- Jacobian matrix ---
J = NaN(2, 2);
J(1,1) = 0.25*(-x1*(1-eta) - x2*(1+eta) + x3*(1+eta) + x4*(1-eta));  %  dxhat/dzeta
J(2,1) = 0.25*(-x1*(1-zeta) + x2*(1-zeta) + x3*(1+zeta) - x4*(1+zeta)); % dxhat/deta
J(1,2) = 0.25*(-y1*(1-eta) - y2*(1+eta) + y3*(1+eta) + y4*(1-eta));  % dyhat/dzeta
J(2,2) = 0.25*(-y1*(1-zeta) + y2*(1-zeta) + y3*(1+zeta) - y4*(1+zeta));  % dyhat/dzeta
J_inv = inv(J);
J_det = det(J);

% -- local deriv. {dNdzeta, dNdeta} ---
loc_derv = 0.25*[-(1-eta), 0, -(1+eta), 0, (1+eta), 0, (1-eta), 0;
    0, -(1-zeta), 0, (1-zeta), 0, (1+zeta), 0, -(1+zeta)];

% -- [B], [D], [B]^T * [D] * B ---
B = zeros(3, 8);
B_first2rows = J_inv*loc_derv;
B(1:2, :) = B_first2rows;
B(3, :) = circshift(B(1,:), [0 1]) + circshift(B(2,:), [0 -1]);

D = (E*(1-mu)/((1+mu)*(1-2*mu)))*[1, mu/(1-mu), 0; mu/(1-mu), 1, 0; 0, 0, (1-2*mu)/(2*(1-mu))]; % for plane strain
BtDB = transpose(B)*D*B;
ki = Wi*J_det*BtDB;

end

