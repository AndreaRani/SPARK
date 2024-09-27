function [lambda, fiedv, label]= algconn_ud(mat, alpha, flag)

% Given a binary graph G(V, E), its algebraic connectivity is the second 
% smallest eigenvalue of its symmetric normalized Laplacian.
% 
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%        - flag, 1 or 0: if set to 1 enables the plot of the eigenvector
%          associated to the algebrai connectivity (i.e. Fiedler vector)
%
% Output: - lambda, algebrai connectivity associated to the input graph
%         - fiedv, eigenvector associated to the algebraic connectivity
%         - label, community labels associated to the Fiedler vector 
%           clustering
%
% 
% Copyright (C) 2024 @SPARK toolbox for (di)graphs
% author: Andrea Ranieri
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


% Ensure mat is binary
mat = double(mat~=0);

[~, symlap, ~]= normlapl_ud(mat, alpha);

[V D W]= eig(symlap);
d= diag(D);
[dmin inddmin]= mink(d, 2);
dmin= sort(dmin, 'ascend');
lambda= dmin(2);

labA= W(:, inddmin(2))>0;
labB= ~labA;
label= cat(2, labA, labB);
fiedv= W(:, inddmin(2));

if flag
    figure; stem(W(:, inddmin(2)));
end
end