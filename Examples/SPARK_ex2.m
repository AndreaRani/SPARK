clc; clear all; close all;

% Set directory
savedir= uigetdir;

%% Generating parameters
nCh= 50; %number of channels
nITER= 100; %number of iterations
Dtot=[20 30 40]; %Total Density
Ratio=[0.45, 0.6, 0.75]; %Ratio intra/inter
alpha_bet= [0.75, 0.5]; %Between-cluster connections imbalance
alpha_wit= [0.75, 0.5]; %Within-cluster connections imbalance

% Size of each partition
nA= nCh/2;
nB= nCh/2;

% PageRank random walk parameter
alpha= .95;

for ds=1:length(Dtot)
    for rt=1:length(Ratio)
        for w= 1:length(alpha_wit)
            for b= 1:length(alpha_bet)
                for n=1:nITER
                    simdata= generate_simdata(nCh, Dtot(ds), Ratio(rt), alpha_wit(w), alpha_bet(b)); % Surrogate network generation

                    %% A priori cluster analysis
                    labA= [ones((0.5*nCh), 1); zeros((0.5*nCh), 1)]; % cluster A affiliation vector
                    labB= [zeros((0.5*nCh), 1); ones((0.5*nCh), 1)]; % cluster B affiliation vector
                    hA= find(labA); hB= find(labB);

                    normas= normassoc_u(simdata, hA, hB); % Normalized association
                    normcut= ncut_u(simdata, hA, hB); % Normalized cut

                    %% Minimum cut partition
                    [lambda, fiedv, label]= algconn_ud(simdata, alpha, 0); % Algebraic connectivity and Fiedler vector
                    newhA= find(label(:, 1)); % min cut affiliation vector cluster A
                    newhB= find(label(:, 2)); % min cut affiliation vector cluster B

                    normas_fv= normassoc_u(simdata, newhA, newhB); % Normalized association
                    normcut_fv= ncut_u(simdata, newhA, newhB); % Normalized cut


                    %% Storage
                    % A priori partition            
                    normass(n,ds,rt,w,b)= normas;
                    nocut(n,ds,rt,w,b)= normcut;

                    % Fiedler partition
                    normass_fv(n,ds,rt,w,b)= normas_fv;
                    nocut_fv(n,ds,rt,w,b)= normcut_fv;

                    clear simdata
                    clear normas normcut
                end
            end
        end
    end
end

tool.note='nIter x Density x Ratio_IntraInter x dir_intra x dir_inter x partition type (a priori/Fiedler)';
tool.normassoc= cat(6, normass, normass_fv);
tool.nocut= cat(6, nocut, nocut_fv);

savename= sprintf('SPARK_simAdj.mat');
save(fullfile(savedir,savename),'tool');
