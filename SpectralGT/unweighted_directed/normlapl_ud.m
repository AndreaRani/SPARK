function [normdiplacian, normsymlap, normsksymlap]= normlapl_ud(mat, alpha)

% Given a directed graph G(V, E), its Diplacian matrix is a generalized
% version of the Laplacian matrix for directed graphs.
% The normalized Diplacian is a normalized version of the Diplacian matrix
% with all ones on the main diagonal.
% 
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%
% Output: - normdiplacian, normalized Diplacian associated to the input graph
%         - normsymlap, symmetric part of the normalized Diplacian
%         - normsksymlap, skew-symmetric part of the normalized Diplacian
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

P= probt_ud(mat, alpha);
statprob= stprob_ud(mat, alpha);

aux= diag(sqrt(statprob));
aux2= 1./aux;
aux2(isinf(aux2))=0;

normdiplacian= aux*(eye(size(P))-P)*aux2;
normsymlap= 0.5.*(normdiplacian+ normdiplacian'); % Symmetric Diplacian
normsksymlap= 0.5.*(normdiplacian- normdiplacian'); % Skew-Symmetric Diplacian
end