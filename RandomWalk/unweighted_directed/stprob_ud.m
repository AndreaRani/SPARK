function statprob= stprob_ud(mat, alpha)

% Given an irreducible Markov Chain with probability transition matrix P,
% its stationary probability distribuion corresponds to a distribution that
% does not change over time. The stationary distribution is the left
% eigenvector that satisfies the equation Wp'*P= Wp'*lambda for lambda=1.
% 
% Input: - mat, input weighted adjacency matrix of G(V, E)
%        - alpha, parameter in [0,1]: (1-alpha) is a correction term
%          representing the probability to jump to any vertex of the
%          network. It is used to implement the probability transition
%          matrix of a PageRank random walker on the input graph
%
% Output: - statprob, stationary probability distribution of the Markov Chain
%           associated to the input graph
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


P= probt_ud(mat, alpha);
TOLERANCE= 1e-2;

% Looking for left-eigenvectors
[Vp lambda Wp]= eig(P);

% Find the eigenvector with real part equal to 1
lambda= diag(round(real(lambda), 4));
ind= find(lambda== double(1));

% If lambda doesn't exist, take the eigenvalue whose real part is closer
% to 1
if length(ind)==0
    ind= find(lambda>= (1-TOLERANCE));
    warning('The real part of the chosen eigenvalue is not exactlyl 1');
    if length(ind)==0
        warning('No stationary distribution was found');
        ind=0;
    end
end
% If there are more lambdas with real part equal to 1, take one uniformly
% at random
if length(ind)>1
    s = RandStream('mlfg6331_64'); 
    ind= randsample(s,[ind],1);
    warning('More than one eigenvalue shares the same real part');
end

statprob= Wp(:, ind);
% Normalization allows to make statprob look like a probability density
statprob= statprob./sum(statprob);
end