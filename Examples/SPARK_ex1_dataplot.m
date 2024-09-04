clear all; close all;

% Set directory
datadir= uigetdir; % datadir contains the SPARK_simVSrand.mat file created using SPARK_ex1

filename= 'SPARK_simVSrand.mat';
fields= {'normassoc', 'algconn', 'nocut', 'relaxt', 'edgemeas'};
names= {'Normalized association', 'Algebraic connectivity', 'Normalized cut', 'Relaxation time', 'Edge measure'};

%% Generating parameters
nCh= 50; %number of channels
nITER= 100; %number of iterations
Dtot=[20 30 40]; %Total Density  
Ratio=[0.45, 0.6, 0.75]; %Ratio intra/inter
alpha_bet= [0.5]; %Between-cluster connections imbalance
alpha_wit= [0.5]; %Within-cluster connections imbalance

tb= loadname(fullfile(datadir, filename));

for ds= 1: length(Dtot)
    savedir= fullfile(datadir, strcat('den_', string(Dtot(ds))));
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
    data_mat= []; % empty tensor
    fi= figure('Visible','on'); fi.Position= [50 50 1750 950];
    for rt=1: length(Ratio)
        for f=1: length(fields)
            currfield= eval(strcat('tb.', fields{f}));
            mean_currfield= mean(currfield(:, rt, ds, 1), 1); mean_randcurrfield= mean(currfield(:, rt, ds, 2), 1);
            curr_data(:, f)= [mean_randcurrfield; mean_currfield];
        end
        data_mat= cat(1, data_mat, curr_data);
    end
    data_mat(5, :)=[];  data_mat(3, :)= [];  % delete random repetitions
    spider_plot(data_mat, 'AxesLabels', names, 'FillOption', {'on', 'on', 'on', 'on'}, 'FillTransparency', [0.2, 0.2, 0.2, 0.2], ...
        'Marker', 'none', 'LineWidth', 3, 'Color', [193, 18, 31; 251, 133, 0; 42, 157, 143; 38, 70, 83]./256, 'LabelFontSize', 14);
    legend('Random', 'rho= 0.45', 'rho= 0.65', 'rho= 0.75' , 'Location', 'bestoutside');
    savename= sprintf('starplot_%s.png', string(Dtot(ds)));
    saveas(gcf, fullfile(savedir, savename));
    close all;
end