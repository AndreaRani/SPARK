function [newlabN1, newlabN2, newmat]= matchLab(labN1, labN2, mat)

% This function rearranges the adjacency matrix mat according to the
% clusters A and B: labN1 and labN2 are the affiliation vectors for 
% clusters A and B respectively.
% Nodes are supposed to be rearranged in order to have nodes in A first
% and then nodes in B (i.e. nodes= [A1, ..., A_nA, B1, ..., B_nB]). The new
% matrix newmat is a block matrix where connections are rearranged as:
% | A -> A   |   B -> A |
% | A- > B   |   B -> B |
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labN1, binary array labelling the nodes in cluster A
%        - labN2, binary array labelling the nodes in cluster B
%
% Output: - newmat, permuted version of mat
%         - newlabN1, labels for cluster A in the permuted matrix newmat
%         - newlabN2, labels for cluster A in the permuted matrix newmat                 
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

% Dimension check
n1= size(mat, 1); n2= size(mat, 2);
if n1~=n2
    error('The input matrix is not squared');
else
    if ( length(labN1)+length(labN2)~= n1)
        error('The given partition does not cover the entire vertex set');
    else
        labN1= reshape(labN1, 1, []);
        labN2= reshape(labN2, 1, []);
        newmat= mat(cat(2, labN1, labN2), cat(2, labN1, labN2));
        newlabN1= find(cat(2, ones(1, length(labN1)), zeros(1, length(labN2))));
        newlabN2= find(cat(2, zeros(1, length(labN1)), ones(1, length(labN2))));
    end
end