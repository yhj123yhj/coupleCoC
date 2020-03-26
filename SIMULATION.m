clear
clc
close all

%input data (take setting 1 as an example)
load('data/rna_simu1.mat');
load('data/atac_simu1.mat');
load('data/rna_clu_simu1.mat');
load('data/atac_clu_simu1.mat');

i = 1;ind = ((i-1)*100+1):(i*100); 
X = atac_simu(:,ind); clu_X_truth = atac_clu_simu(:,i); % auxiliary data
Y = rna_simu(:,ind); clu_Y_truth = rna_clu_simu(:,i);  % target data

%%initialize the parameters
nclu_X = 2;nclu_Y = 2;nclu_Z = 3;alpha=3;lambda = [1,2.5,1];
niter_outer = 3;niter_inter = 6;

%coupleCoC
[~,~,~,~,~, clu_X_bes,clu_Y_bes] = coupleCoC(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer, niter_inter, lambda, alpha, i);

%Evaluation
[TAB_X, TAB_Y, Eval_tab] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes, clu_Y_bes);



       