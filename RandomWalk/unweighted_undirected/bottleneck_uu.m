function [bott_A, bott_B]= bottleneck_uu(mat, alpha, labA, labB)

% Given an ergodic chain with stationary probability distribution Wp and
% vertex set partition {A, B}, the bottleneck ratio is defined as the ratio 
% between two terms: the probability to leave a certain set of nodes and the  
% probability to move everywhere starting from the same subset, both 
% calculated starting from Wp.
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%
% Output: - bott_A, bottleneck ratio for the nodes in A
%         - bott_B, bottleneck ratio for the nodes in B
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

[QA2B, QB2A, Qleave]= Qleave_uu(newmat, hA, hB, alpha); % One-step probability to leave the partition @stationary state
statprob= stprob_uu(newmat, hA, hB, alpha);

stats= H'*diag(statprob)*H;

bott_A= QA2B/stats(1, 1);
bott_B= QB2A/stats(2, 2);

end