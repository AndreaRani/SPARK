function normcut= ncut_w(mat, labA, labB)

% Given a weighted directed graph G(V, E) and a partition partition of the vertex
% set {A, B} such that V(A) U V(B) = V, its normalized cut is the sum of
% two terms: each one is given by the ratio of the total cut induced by the
% partition and the association of each subset of nodes.
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%
% Output: - normcut, normalized cut for the vertex set partition V= {A, B}
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


[assocA, assocB]= assoc_w(mat, labA, labB);
[~, ~, totcut]= cut_wd(mat, labA, labB);

% Normalized cut
normcut= (totcut/assocA)+ (totcut/assocB);
end