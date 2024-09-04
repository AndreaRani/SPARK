clc; clear all; close all;

% Set directory
savedir= uigetdir;

%% Generating parameters
nCh= 50; %number of channels
nITER= 100; %number of iterations
Dtot=[20 30 40]; %Total Density  
Ratio=[0.45, 0.6, 0.75]; %Ratio intra/inter
alpha_bet= [0.5]; %Between-cluster connections imbalance
alpha_wit= [0.5]; %Within-cluster connections imbalance

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

                    % List of measures
                    normas= normassoc_u(simdata, hA, hB); % Normalized association
                    normcut= ncut_u(simdata, hA, hB); % Normalized cut
                    [relt, specgap]= spgap_ud(simdata, alpha); % Relaxation time
                    [~, ~, Qleave]= Qleave_ud(simdata, hA, hB, alpha);% Q leave
                    [lambda, fiedv, label]= algconn_ud(simdata, alpha, 0); % Algebraic connectivity


                    %% Random Graph
                    randnet= makerandCIJ_dir(nCh, sum(sum(simdata))); %Random network generation
                    
                    % List of measures
                    rand_normas= normassoc_u(randnet, hA, hB); % Normalized association
                    rand_normcut= ncut_u(randnet, hA, hB); % Normalized cut
                    [rand_relt, ~]= spgap_ud(randnet, alpha); % Relaxation time
                    [~, ~, rand_Qleave]= Qleave_ud(randnet, hA, hB, alpha);% Q leave
                    [rand_lambda, ~, ~]= algconn_ud(randnet, alpha, 0); % Algebraic connectivity

                    
                    %% Storage
                    normass(n,rt, ds)= normas;
                    algcon(n,rt, ds)= lambda;
                    nocut(n,rt, ds)= normcut;
                    relaxt(n,rt, ds)= relt;
                    edgemeas(n,rt, ds)= Qleave;

                    rand_normass(n,rt, ds)= rand_normas;
                    rand_algcon(n,rt, ds)= rand_lambda;
                    rand_nocut(n,rt, ds)= rand_normcut;
                    rand_relaxt(n,rt, ds)= rand_relt;
                    rand_edgemeas(n,rt, ds)= rand_Qleave;

                    clear simdata
                    clear normas lambda normcut relt Qleave
                    clear rand_normas rand_lambda rand_normcut rand_relt rand_Qleave
                end
            end
        end
    end
end

tool.note='nIter x Ratio x Density x type (modular/rand)';
tool.normassoc= cat(4, normass, rand_normass);
tool.algconn= cat(4, algcon, rand_algcon);
tool.nocut= cat(4, nocut, rand_nocut);
tool.relaxt= cat(4, relaxt, rand_relaxt);
tool.edgemeas= cat(4, edgemeas, rand_edgemeas);

savename= 'SPARK_simVSrand.mat';
save(fullfile(savedir,savename),'tool');
