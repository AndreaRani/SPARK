function [assocA, assocB]= assoc_u(mat, labA, labB)

% Given an unweighted graph G(V, E), the association of a subset of nodes 
% is the sum of the connections involving that subset of nodes.
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%
% Output: - assocA, association of the nodes in A
%         - assocB, association of the nodes in B
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

% Rows and columns permutation according to the affiliation vectors labA
% and labB
[hA, hB, newmat]= matchLab(labA, labB, mat);

% Partition matrix
H= zeros(size(newmat, 1), 2);
H(hA, 1)= 1;
H(hB, 2)=1;

% Association and cuts
assocs= sum(H'*newmat*H, 2);

% [volA, volB]= vol(mat, labA, labB);
% [cutA2B, cutB2A, ~]= cut_ud(mat, labA, labB);

% Association
assocA= assocs(1);
assocB= assocs(2);
end