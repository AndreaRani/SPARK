function [diplacian, symlap, sksymlap]= lapl_wd(mat, alpha)

% Given a directed graph G(V, E), its Diplacian matrix is a generalized
% version of the Laplacian matrix for directed graphs.
% The Diplacian can be decomposed into its symmetric and skew-symmetric
% part, namely symlap and sksymlap.
% 
% Input: - mat, input adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%
% Output: - diplacian, Diplacian matrix associated to the input graph
%         - symlap, symmetric part of the Diplacian matrix
%         - sksymlap, skew-symmetric part of the Diplacian matrix
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


P= probt_wd(mat, alpha);
statprob= stprob_wd(mat, alpha);

diplacian= diag(statprob)*(eye(size(mat))-P);
symlap= 0.5.*(diplacian+ diplacian'); % Symmetric Diplacian
sksymlap= 0.5.*(diplacian- diplacian'); % Skew-Symmetric Diplacian
end