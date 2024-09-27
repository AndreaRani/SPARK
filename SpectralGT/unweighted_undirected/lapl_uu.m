function lapl= lapl_uu(mat, alpha)

% Given an undirected graph G(V, E), its Laplacian is a symmetric matrix 
% describing the information flow within the graph.
% 
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%
% Output: - lapl, Laplacian matrix associated to the input graph
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

P= probt_uu(mat, alpha);
statprob= stprob_uu(mat, alpha);

lapl= diag(statprob)*(eye(size(mat))-P);
end