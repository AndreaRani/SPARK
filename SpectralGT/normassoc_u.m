function normas= normassoc_u(mat, labA, labB)

% Given a binary graph G(V, E), the normalized association of a subset of nodes
% is the ratio between the volume and the association of that subset. When
% a partition {A, B} of the vertex set is given such that V(A) U V(B) = V,
% the normalized association is the sum of the normalized association for A
% and B.
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%
% Output: - normas, normalized association for the vertex set partition 
%                   V= {A, B}
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
assocs= H'*newmat*H;

% Association
[assocA, assocB]= assoc_u(newmat, hA, hB); 

% Normalized association
normas= (assocs(1, 1)/assocA)+ (assocs(2, 2)/assocB);
end