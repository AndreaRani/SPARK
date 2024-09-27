function [relt, specgap]= spgap_wd(mat, alpha)

% Given an irreducible Markov Chain with probability transition matrix P,
% its spectral gap is the difference between 1 and the largest eigenvalue 
% of P without considering the one related to the stationary distribution
% of the chain (i.e. lambda= 1).
% The relaxation time of a Markov Chain gives an estimate of the time
% needed by the chain to reach its stationary distribution: it is the 
% reversal of 1 minus the spectral gap.
% 
% Input: - mat, input binary adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. It is used to implement the probability transition
%          matrix of a PageRank random walker on the input graph
%
% Output: - relt, relaxation time of the Markov Chain associated to the
%           input graph
%         - specgap, spectral gap of the Markov Chain associated to the
%           input graph
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

[V D W]= eig(P);
d= diag(round(abs(D), 4));
[dmax inddmax]= maxk(d, 2);

if (dmax(1)~=double(1))
    relt= 1/(d(inddmax(1)));
    specgap= 1- d(inddmax(1));
else
    relt= 1/(d(inddmax(2)));
    specgap= 1- d(inddmax(2));
end

end