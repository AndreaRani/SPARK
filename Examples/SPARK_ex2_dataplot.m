clear all; close all; clc;

% Set directory
datadir= uigetdir; % datadir contains the SPARK_simAdj.mat file created using SPARK_ex2

filename= 'SPARK_simAdj.mat';
fields= {'normassoc', 'nocut'};
names= {'Normalized association', 'Normalized cut'};

%% Generating parameters
nCh= 50; %number of channels
nITER= 100; %number of iterations
Dtot=[20 30 40]; %Total Density
Ratio=[0.45, 0.6, 0.75]; %Ratio intra/inter
alpha_bet= [0.75, 0.5]; %Between-cluster connections imbalance
alpha_wit= [0.75, 0.5]; %Within-cluster connections imbalance

tb= loadname(fullfile(datadir, filename));

for f=1: length(fields)
    currsavedir= fullfile(datadir, fields{f});
    if ~exist(currsavedir, 'dir')
        mkdir(currsavedir);
    end
    field= eval(strcat('tb.', fields{f}));
    for d=1: length(Dtot)
        figure('Position',[50 50 1500 900]);
        for w=1: length(alpha_wit)
            for b=1: length(alpha_bet)
                currfield_pr= squeeze(field(:, d, :, w, b, 1));
                currfield_fv= squeeze(field(:, d, :, w, b, 2));
                % Comparative measure
                x= (currfield_fv-currfield_pr)./(currfield_pr+currfield_fv);

                %% Boxplot
                x= x(:);
                ratio= [repmat({'0.45'}, [1, size(currfield_pr, 1)]), repmat({'0.6'}, [1, size(currfield_pr, 1)]), repmat({'0.75'}, [1, size(currfield_pr, 1)])];
                d_intra= repmat([repmat( alpha_wit(w), 1, length(ratio))], 1, length(Dtot));
                d_inter= repmat([repmat( alpha_bet(b), 1, length(ratio))], 1, length(Dtot));

                g(b, w)= gramm('x', ratio, 'y', x(:));
                g(b, w).stat_boxplot();
                g(b, w).set_color_options('map', ([0 0.4470 0.7410])); 
                g(b, w).set_names('x','within/between','y', '', 'color', 'Legend:');
                g(b, w).set_text_options('base_size', 14);
                g(b, w).set_title(append('alpha_bet= ', string(alpha_bet(b)), ', alpha_wit= ', string(alpha_wit(w))), 'FontSize', 18);
            end
        end
        g.set_title(append(names{f}), 'FontSize', 22);
        g.draw();
        exportfilename= strcat(fields{f}, '_', string(Dtot(d)), '_compMeas');
        g.export('file_name', exportfilename,'file_type', 'png', 'export_path', currsavedir);
        clear g; close all;
    end
end
