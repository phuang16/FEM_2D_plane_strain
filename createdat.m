% This script create dat.file that records:
% (1) the coordinate for all nodes
% (2) the order of the nodes in each element

% === coordinate for all nodes ===
% (x, y)
node_loc_all =  [0, 0;
    0, 5/2;
    0, 5;
    10/3, 0;
    10/3, 5/2;
    10/3, 5;
    20/3, 0;
    20/3, 5/2;
    10, 0;
    10, 5]

% == order of the nodes in each element ===
% element A, B, C, D, E
node_order_all = [1, 2, 5, 4;   % corresponding to local nodes 1, 2, 3, 4 of an element
    2, 3, 6, 5;
    4, 5, 8, 7;
    5, 6, 10, 8;
    7, 8, 10, 9];

dlmwrite('globalnode_all.dat', node_loc_all)
dlmwrite('localnode_order_all.dat', node_order_all)