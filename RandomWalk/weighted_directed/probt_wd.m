function P= probt_wd(mat, alpha)

% Given a weighted directed graph G(V, E), its probability transition matrix
% P contains the probability of a random walk to jump from a node to one of 
% its neighbours.
% As to guarantee the existance of a stationary distribution on the
% associated Markov chain and to overcome the problem of dangling nodes
% (i.e. to deal with an irreducible Markov Chain), we preferred to adopt
% the PageRank random walk algorithm with respect to the classic one.
% 
% Input: - mat, input weighted adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. The higher the value of alpha, the more accurately the 
%          topology will be preserved
%
% Output: - P, probability transition matrix of the PageRank random walk
%           on the input graph
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


e= ones(size(mat, 1), 1);
so = sum(mat,2)'; % out-strength
a= ~(so>0)';
n= length(so);

pinvso= pinv(diag(so)); % out-strength pseudo inverse matrix

P= alpha.*((pinvso*mat)+ (1/n).*(a*e'))+ ((1-alpha)*(1/n)).*(e*e');
end