clear
clc
close all

%%input data simu1-->setting1; simu3-->setting2; simu2-->setting3
load('rna_simu3.mat');load('atac_simu3.mat');load('rna_clu_simu3.mat');load('atac_clu_simu3.mat');

%%initialize the parameters
maxiter=50;p=100;nclu_Z = 3;alpha=3;lambda1 = [1,2.5,1];lambda2 = [1,2.5,1];
niter_outer1 = 3;niter_inter1 = 6;nclu_X = 2;nclu_Y = 2;
niter_outer2 = 3;niter_inter2 = 6;
Eval = zeros(40,maxiter);

for i = 44:50
    ind = ((i-1)*p+1):(i*p);
    X = atac_simu(:,ind); clu_X_truth = atac_clu_simu(:,i); % source data
    Y = rna_simu(:,ind); clu_Y_truth = rna_clu_simu(:,i);  % target data
    
    %%stc method
    [~,~,~,~,clu_X_bes1,clu_Y_bes1] = stc(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer1, niter_inter1, lambda1, i);
    [TAB_X1, TAB_Y1, Eval_tab1] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes1, clu_Y_bes1);

    %%coupleCoC method
    [~,~,~,~,~, clu_X_bes2,clu_Y_bes2] = coupleCoC(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer2, niter_inter2, lambda2, alpha, i);
    [TAB_X2, TAB_Y2, Eval_tab2] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes2, clu_Y_bes2);

     %%kmeans
    [idx,~] = kmeans(X,nclu_X,'MaxIter',10000,'Replicates',15);
    [idy,~] = kmeans(Y,nclu_Y,'MaxIter',10000,'Replicates',15);
    [TAB_X3, TAB_Y3, Eval_tab3] = clu_eval(clu_X_truth, clu_Y_truth, idx, idy);
    
    %%Spectral clustering
    [C_X, ~, ~] = SpectralClustering(X,nclu_X,5,'kmean',[2 2]);
    [C_Y, ~, ~] = SpectralClustering(Y,nclu_Y,5,'kmean',[2 2]);
    [TAB_X4, TAB_Y4, Eval_tab4] = clu_eval(clu_X_truth, clu_Y_truth, C_X, C_Y);
    
    %%Hierarchial clustering
    Metric   = @(X,Y)norm(X-Y);Colormap = 'spring';
    [~, Clu_X2, ~] = Hierarchical_clustering(X,'UPGMA',Metric,'Number',nclu_X,Colormap);
    [~, Clu_Y2, ~] = Hierarchical_clustering(Y,'UPGMA',Metric,'Number',nclu_Y,Colormap);
    C_X2 = zeros(size(X,1),1); C_X2(Clu_X2{1,:})=1;C_X2(Clu_X2{2,:})=2;
    C_Y2 = zeros(size(Y,1),1); C_Y2(Clu_Y2{1,:})=1;C_Y2(Clu_Y2{2,:})=2;
    [TAB_X5, TAB_Y5, Eval_tab5] = clu_eval(clu_X_truth, clu_Y_truth, C_X2, C_Y2);
    
    %%save the results
    Eval(1:4,i) = Eval_tab1{:,:}(:,1);Eval(21:24,i) = Eval_tab1{:,:}(:,2);
    Eval(5:8,i) = Eval_tab2{:,:}(:,1);Eval(25:28,i) = Eval_tab2{:,:}(:,2);
    Eval(9:12,i) = Eval_tab3{:,:}(:,1);Eval(29:32,i) = Eval_tab3{:,:}(:,2);
    Eval(13:16,i) = Eval_tab4{:,:}(:,1);Eval(33:36,i) = Eval_tab4{:,:}(:,2);
    Eval(17:20,i) = Eval_tab5{:,:}(:,1);Eval(37:40,i) = Eval_tab5{:,:}(:,2);
end

% case 1: Eval_simu1 50 iters done;
% case 2: Eval_simu2 50 iters done;
% case 3: Eval_simu3 50 iters done;
save('Eval_simu3_9Feb.mat', 'Eval');
%load('Eval_simu1_9Feb.mat');
%mean(Eval(21:40,1:50),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameter tuning
i=1;p=100;
ind = ((i-1)*p+1):(i*p);
X = atac_simu(:,ind); clu_X_truth = atac_clu_simu(:,i); % source data
Y = rna_simu(:,ind); clu_Y_truth = rna_clu_simu(:,i);  
RangeZ = 3; Lamb = 2.5;
obj_re_mat = zeros(size(RangeZ,2),size(Lamb,2));
for j = 1:size(RangeZ,2)
    for k = 1:size(Lamb,2)
        nclu_Z = RangeZ(j);
        lambda1 = [1,Lamb(k),1];
        [~,~,~,obj_re,clu_X_bes1,clu_Y_bes1] = stc_main1(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer1, niter_inter1, lambda1, i);
        [TAB_X1, TAB_Y1, Eval_tab1] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes1, clu_Y_bes1);
        obj_re_mat(j,k) = Eval_tab1{:,:}(1,2);
    end
end

Alpha1 = 3; nclu_Z=3;lambda2=[1,2.5,1];
obj_re_mat2 = zeros(1,size(Alpha1,2));
for s = 1:size(Alpha1,2)
    alpha = Alpha1(s);
 [~,~,~,~,~, clu_X_bes2,clu_Y_bes2] = stc_main2(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer2, niter_inter2, lambda2, alpha, i);
 [TAB_X2, TAB_Y2, Eval_tab2] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes2, clu_Y_bes2);
 obj_re_mat2(s) = Eval_tab2{:,:}(1,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Eval_simu1.mat');
%%input data
load('rna_simu1.mat');load('atac_simu1.mat');load('rna_clu_simu1.mat');load('atac_clu_simu1.mat');

%%initialize the parameters
maxiter=50;p=100;
for i = 1:maxiter
    ind = ((i-1)*p+1):(i*p);
    Y = atac_simu3(:,ind); clu_Y_truth = atac_clu_simu3(:,i); % source data
    X = rna_simu3(:,ind); clu_X_truth = rna_clu_simu3(:,i);  % target data
    
    %%Spectral clustering
    [C_X, ~, ~] = SpectralClustering(X,nclu_X,5,'kmean',[2 2]);
    [C_Y, ~, ~] = SpectralClustering(Y,nclu_Y,5,'kmean',[2 2]);
    [TAB_X4, TAB_Y4, Eval_tab4] = clu_eval(clu_X_truth, clu_Y_truth, C_X, C_Y);
    
    %%Hierarchial clustering
    Metric   = @(X,Y)norm(X-Y);Colormap = 'spring';
    [~, Clu_X2, ~] = Hierarchical_clustering(X,'UPGMA',Metric,'Number',nclu_X,Colormap);
    [~, Clu_Y2, ~] = Hierarchical_clustering(Y,'UPGMA',Metric,'Number',nclu_Y,Colormap);
    C_X2 = zeros(size(X,1),1); C_X2(Clu_X2{1,:})=1;C_X2(Clu_X2{2,:})=2;
    C_Y2 = zeros(size(Y,1),1); C_Y2(Clu_Y2{1,:})=1;C_Y2(Clu_Y2{2,:})=2;
    [TAB_X5, TAB_Y5, Eval_tab5] = clu_eval(clu_X_truth, clu_Y_truth, C_X2, C_Y2);
end















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate data
clu_X_mat = zeros(niter_outer,size(X,1));
clu_Y_mat = zeros(niter_outer,size(Y,1));
clu_Z_mat = zeros(niter_outer,size(X,2));

for l = 1:niter_outer
clu_X_mat(l,:) = randsample(nclu_X, size(X,1),true);
clu_Y_mat(l,:) = randsample(nclu_Y, size(Y,1),true);
clu_Z_mat(l,:) = randsample(nclu_Z, size(X,2),true);
end

%save('clu_X_mat_1to400.mat', 'clu_X_mat');
%save('clu_Y_mat_1to400.mat', 'clu_Y_mat');
%save('clu_Z_mat_1to400.mat', 'clu_Z_mat');

%%%
% Non-JSD method and JSD_v1 method
T1 = zeros(nclu_X,nclu_X,niter_outer);
T2 = zeros(nclu_X,nclu_X,niter_outer);
T3 = zeros(nclu_Y,nclu_Y,niter_outer);

clu_X_mat_re = zeros(niter_inter+2,size(X,1),niter_outer);
clu_Y_mat_re = zeros(niter_inter+1,size(Y,1),niter_outer);
clu_Z_mat_re = zeros(niter_inter+1,size(X,2),niter_outer);
obj_mat_re = zeros(niter_inter+1,niter_outer);

T12 = zeros(nclu_X,nclu_X,niter_outer);
T22 = zeros(nclu_X,nclu_X,niter_outer);
T32 = zeros(nclu_Y,nclu_Y,niter_outer);

clu_X_mat_re2 = zeros(niter_inter+1,size(X,1),niter_outer);
clu_Y_mat_re2 = zeros(niter_inter+1,size(Y,1),niter_outer);
clu_Z_mat_re2 = zeros(niter_inter+1,size(X,2),niter_outer);
obj_mat_re2 = zeros(niter_inter,niter_outer);
JSD_mat_re2 = zeros(niter_inter,niter_outer);

for i = 1:niter_outer
 tic
    [clu_X_Trace, clu_Y_Trace, clu_Z_Trace, obj_val_Trace] = stc_fun1(X, Y, clu_X_mat(i,:), clu_Y_mat(i,:), clu_Z_mat(i,:), nclu_X, nclu_Y, nclu_Z, niter_inter, lambda);
 toc
    clu_X_mat_re(:,:,i) = clu_X_Trace;
    clu_Y_mat_re(:,:,i) = clu_Y_Trace;
    clu_Z_mat_re(:,:,i) = clu_Z_Trace;
    obj_mat_re(:,i) = obj_val_Trace;
    T1(:,:,i) = crosstab(clu_X_truth,clu_X_Trace(niter+1,:));
    T2(:,:,i) = crosstab(clu_X_truth,clu_X_Trace(niter+2,:));
    T3(:,:,i) = crosstab(clu_Y_truth,clu_Y_Trace(niter+1,:));
tic
    [clu_X_Trace, clu_Y_Trace, clu_Z_Trace, obj_val_Trace, JSD_Trace] = stc_fun2(X, Y, clu_X_mat(i,:), clu_Y_mat(i,:), clu_Z_mat(i,:), nclu_X, nclu_Y, nclu_Z, niter, lambda, alpha1);
toc
    clu_X_mat_re2(:,:,i) = clu_X_Trace;
    clu_Y_mat_re2(:,:,i) = clu_Y_Trace;
    clu_Z_mat_re2(:,:,i) = clu_Z_Trace;
    obj_mat_re2(:,i) = obj_val_Trace;
    JSD_mat_re2(:,i) = JSD_Trace;
    T12(:,:,i) = crosstab(clu_X_truth(:,1),clu_X_Trace(niter+1,:));
    T32(:,:,i) = crosstab(clu_Y_truth(:,1),clu_Y_Trace(niter+1,:));
end

%[AR1,RI1,MI1,HI1]=RandIndex(exp_true,clu_X_Trace(niter+1,:));
%[AR2,RI2,MI2,HI2]=RandIndex(exp_true,clu_X_Trace(niter+2,:));
%[AR3,RI3,MI3,HI3]=RandIndex(acc_true,clu_Y_Trace(niter+1,:));

[idx,Cx] = kmeans(X,2);[idy,Cy] = kmeans(Y,2);
crosstab(clu_X_truth,idx);crosstab(clu_Y_truth,idy)
KL_Z = obj_mat_re(niter+1,1:I);
find(KL_Z == min(KL_Z), 1 )
obj_mat_re(:,1:7);
j =18;[T1(:,:,j),T2(:,:,j),T3(:,:,j),T12(:,:,j),T32(:,:,j)];
cMat = zeros(4,7);
[cMat(1,1),cMat(2,1),cMat(3,1),cMat(4,1)]=RandIndex(clu_X_truth,clu_X_mat_re(niter+1,:,j));
[cMat(1,2),cMat(2,2),cMat(3,2),cMat(4,2)]=RandIndex(clu_X_truth,clu_X_mat_re(niter+2,:,j));
[cMat(1,3),cMat(2,3),cMat(3,3),cMat(4,3)]=RandIndex(clu_Y_truth,clu_Y_mat_re(niter+1,:,j));
[cMat(1,4),cMat(2,4),cMat(3,4),cMat(4,4)]=RandIndex(clu_X_truth,clu_X_mat_re2(niter+1,:,j));
[cMat(1,5),cMat(2,5),cMat(3,5),cMat(4,5)]=RandIndex(clu_Y_truth,clu_Y_mat_re2(niter+1,:,j));
[cMat(1,6),cMat(2,6),cMat(3,6),cMat(4,6)]=RandIndex(clu_X_truth,idx);
[cMat(1,7),cMat(2,7),cMat(3,7),cMat(4,7)]=RandIndex(clu_Y_truth,idy);
%%%
% JSD_v2 method
for i = 1:1
    tic
    %[clu_X_Trace, clu_Y_Trace, clu_Z_Trace, obj_val_Trace, JSD_Trace] = stc_fun2(X, Y, clu_X_mat(i,:), clu_Y_mat(i,:), clu_Z_mat(i,:), nclu_X, nclu_Y, nclu_Z, niter, lambda, alpha1);
    toc
    %clu_X_mat_re2(:,:,i) = clu_X_Trace;
    %clu_Y_mat_re2(:,:,i) = clu_Y_Trace;
    %clu_Z_mat_re2(:,:,i) = clu_Z_Trace;
    %obj_mat_re2(:,i) = obj_val_Trace;
    %JSD_mat_re2(:,i) = JSD_Trace;
    %T12(:,:,i) = crosstab(exp_true,clu_X_Trace(niter+1,:));
    %T32(:,:,i) = crosstab(acc_true,clu_Y_Trace(niter+1,:));
    %T12(:,:,i) = crosstab(clu_X_truth,clu_X_Trace(niter+1,:));
    %T32(:,:,i) = crosstab(clu_Y_truth,clu_Y_Trace(niter+1,:));
end
