function [QA2A, QB2B, Qrest]= Qrest_ud(mat, labA, labB, alpha)

% Given an ergodic chain with probability transition matrix P, partition 
% of the vertex set {A, B} and stationary distribution Wp, the probability 
% of a random walker to rest into its starting vertex subset in one step once 
% reached the stationarity, is the sum of all resting probabilities of the
% updated probability matrix Q= diag(Wp)*P.
%
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - labA, binary array labelling for the nodes belonging to the subset A
%        - labB, binary array labelling for the nodes belonging to the subset B
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%
% Output: - onesfs_leaveA, probability to leave cluster A in one step from
%           stationarity
%         - onesfs_leaveB, probability to leave cluster B in one step from
%           stationarity
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


% Ensure mat is binary before compute the degree of each node
mat = double(mat~=0);

% Rows and columns permutation according to the affiliation vectors labA
% and labB
[hA, hB, newmat]= matchLab(labA, labB, mat);
% Partition matrix
H= zeros(size(newmat, 1), 2);
H(hA, 1)= 1;
H(hB, 2)=1;

P= probt_ud(mat, alpha);
statprob= stprob_ud(mat, alpha);

% Update the transition matrix
Q= diag(statprob)*P;
% Probability of the imposed partition
sumPart= H'*Q*H;

QA2A= sumPart(1, 1); 
QB2B= sumPart(2, 2);
Qrest= QB2B+ QA2A;
end