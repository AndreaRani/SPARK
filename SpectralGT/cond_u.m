function [condA, condB]= cond_u(mat, labA, labB)

% Given a directed binary graph G(V, E) and a partition partition of the vertex
% set {A, B} such that V(A) U V(B) = V, the conductance of each subset of
% nodes is the ratio between the cut and the association of that subset.
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%
% Output: - condA, conductance of the nodes in A
%         - condB, conductance of the nodes in B
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

[~, ~, totcut]= cut_ud(mat, labA, labB);
[assocA, assocB]= assoc_u(mat, labA, labB);

% Conductance
condA= totcut/assocA;
condB= totcut/assocB;
end