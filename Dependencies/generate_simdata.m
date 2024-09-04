function adj= generate_simdata(nCh, D, ratio, alpha_w, alpha_b)

% This function generates synthetic adjacency matrices mimicking the
% interaction of two clusters in the underlying network. Cluster-to-cluster
% interactions are tuned by the generating parameters Ratio, alpha_w 
% and alpha_b.
% 
% Input: - nCh, number of nodes in the underlying network
%        - D, network's density
%        - Ratio, fraction of within-cluster links
%        - alpha_w, within-cluster connection imbalance
%        - alpha_b, between-cluster connection imbalance
%
% Output: - adj, synthetic adiacency matrix for the underlying network
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


if rem(nCh,2)==0 % control on the number of channels 
    maxconn= nCh*(nCh-1); % maximum number of connections (avoiding self-loops)
    nc= maxconn*(D/100); % number of existing connections
    nc_wit= round(ratio*nc); % number of intra-cluster connections
    nc_bet=round((1-ratio)*nc); % number of inter-cluster connections
    
    %% Within-cluster blocks generation
    A11 = zeros(nCh/2); A22 = zeros(nCh/2);
    
    index_i=1:(nCh/2)^2;                   
    index_i(1:(nCh/2+1):end)=[]; % avoid self-loops
    shuffled_A11 = index_i(randperm(size(index_i,2)));
    shuffled_A22 = index_i(randperm(size(index_i,2)));
    dir_rand_INTRA= [round(alpha_w*nc_wit),round((1-alpha_w)*nc_wit)]; % number of within-cluster connections in each block

    index1=[randi(numel(dir_rand_INTRA)) randi(numel(dir_rand_INTRA))];
    while index1(1)==index1(2)
        index1=[randi(numel(dir_rand_INTRA)) randi(numel(dir_rand_INTRA))];
    end
%     index1= [1 2]; % To maintain the same direction in cluster imbalance
    
    A11(shuffled_A11(1:dir_rand_INTRA(index1(1))))=1;
    A22(shuffled_A22(1:dir_rand_INTRA(index1(2))))=1;
    
    %% Between-cluster blocks generation
    A12=zeros(nCh/2); A21=zeros(nCh/2);
    
    index1_2 = randperm((nCh/2)^2);
    index2_1 = randperm((nCh/2)^2);
    dir_rand_INTER= [round(alpha_b*nc_bet),round((1-alpha_b)*nc_bet)]; % number of between-cluster connections in each block

    index2=[randi(numel(dir_rand_INTER)) randi(numel(dir_rand_INTER))];
    while index2(1)==index2(2)
        index2=[randi(numel(dir_rand_INTER)) randi(numel(dir_rand_INTER))];
    end
%     index2= [2 1]; % To maintain the same direction in cluster imbalance
    
    A12(index1_2(1:dir_rand_INTER(index2(1))))=1;
    A21(index2_1(1:dir_rand_INTER(index2(2))))=1;
    
    %% Adjacency matrix construction
    adj= [A11, A12; A21, A22];    
else
    disp('odd number of channels');
end
