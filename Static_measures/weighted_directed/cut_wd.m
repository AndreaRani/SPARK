function [cutA2B cutB2A totcut]= cut_wd(mat, labA, labB)

% Given a directed graph G(V, E) and a partition partition of the vertex
% set {A, B} such that V(A) U V(B) = V, its cut is defined as the number
% of crossing connections between A and B.
%
% Input: - mat, input weighted adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%
% Output: - cutA2B, number of crossing connections from A to B
%         - cutB2A, number of crossing connections from B to A
%         - totcut, total number of crossing connections between A and B
% 
% 
% Copyright (C) 2024 @SPARK toolbox for (di)graphs
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


% Rows and columns permutation according to the affiliation vectors labA
% and labB
[hA, hB, newmat]= matchLab(labA, labB, mat);

% Partition matrix
H= zeros(size(newmat, 1), 2);
H(hA, 1)= 1;
H(hB, 2)=1;

% Cut
volCutadj= H'*newmat*H;
cutA2B= volCutadj(2, 1);
cutB2A= volCutadj(1, 2);
totcut= cutA2B+ cutB2A;
end