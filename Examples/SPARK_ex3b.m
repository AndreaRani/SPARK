close all; clear all; clc;

datadir= uigetdir;
filename= uigetfile(datadir);
filepath= fullfile(datadir, filename);
alpha= 0.99;

rawdata= readtable(filepath);
data= rawdata(:, 2:(end-1));
data= table2array(data);

gt= table2array(rawdata(:, end));
cam= strcmp('Cammeo', gt);
gt= cam;

% Discard missings
row_nani= sum(data, 2);
nani= isnan(row_nani);
todel= find(nani);
data(todel, :)=[];
gt(todel, :)= [];

s1= sum(gt); s0= sum(~gt);

%% k-fold cross validation

% Find the largest class to balance the dataset
if (s1-s0)>0
    sG_elem= find(gt); sg_elem= find(~gt);
    sG= sum(gt(sG_elem)); sg= sum(~gt(sg_elem));
else
    sG_elem= find(~gt); sg_elem= find(gt);
    sG= sum(~gt(sG_elem)); sg= sum(gt(sg_elem));
end

ratio= sg/sG;
k= round(1/ratio);

for s=1: 50
    for l=1:k
        clear map_wei map_adj C fied
        start_ind= (l-1)*length(sg_elem)+1;
        end_ind= l*length(sg_elem);
        if end_ind>length(sG_elem)
            end_ind= length(sG_elem);
        end
        sG_lf_ind= sG_elem(start_ind:end_ind);
        data_ind= cat(1, sG_lf_ind, sg_elem);
        data_randind= data_ind(randperm(length(data_ind)));
        coord= data(data_randind, :);
        sigma2= var(coord(:));

        % Data centering
        m1= mean(coord, 1);
        DCmat= repmat(m1, size(coord, 1), 1);
        sigma= std(coord, 1);
        sigma= diag(sigma);
        % Standardized Euclidean distance
        coord2= (coord-DCmat)*inv(sigma);

        % Distance matrix
        for j=1: size(coord2, 1)
            for f=1: size(coord2, 1)
                if f>j
                    map_wei(j, f)= sum((coord2(j, :)- coord2(f, :)).^2);
                else
                    map_wei(j, f)= 0;
                end
            end
        end
        map_wei= map_wei+map_wei';
        % Diffusion kernel
        map_wei= exp(-map_wei);
        map_wei= map_wei-diag(diag(map_wei)); %avoid self-loops
        % Thresholding
        map_wei(map_wei< prctile(map_wei, 50,'all'))= 0;

        % Adjacency matrix
        map_adj= zeros(size(map_wei));
        map_adj(map_wei~=0)= 1;
        map_adj= map_adj-diag(diag(map_adj)); %avoid self-loops

        tic
        % Fiedler vector
        [~, ~, fied]= algconn_uu(map_adj, alpha, 0);
        t(s+ (l-1))= toc;

        % fied actually stores the labels for the two classes separately:in
        % practice, just one of the two columns is needed
        
        if (fied(:, 1)'*gt(data_randind))<(fied(:, 2)'*gt(data_randind)) % gt labels match with the 0s of the Fiedler vector
            fied= fied(:, 2);
        else
            fied= fied(:, 1);
        end
        
        C = confusionmat(gt(data_randind), fied(:, 1));
        confusionchart(C)
        acc(s)= trace(C)/(sum(C, 'all'));
        prec(s)= C(2, 2)/(C(1, 2)+C(2, 2));
        rec(s)= C(2, 2)/(C(2, 1)+C(2, 2));
        F1(s)= (2*prec(s)*rec(s))/(prec(s)+rec(s));

        figure; imagesc(map_adj); title('Raw');
        label_1= find(gt(data_ind)); label_0= find(~gt(data_ind));
        lab_gt= cat(1, label_1, label_0);
        figure; imagesc(map_adj(lab_gt, lab_gt)); title('Label');
        fiedlab1= find(fied>0); fiedlab2= find(~fied>0);
        fiedlab= cat(1, fiedlab1, fiedlab2);
        figure; imagesc(map_adj(fiedlab, fiedlab));
    end
    acc_m= mean(acc);
    prec_m= mean(prec);
    rec_m= mean(rec);
    F1_m= mean(F1);
    compt= mean(t);
end