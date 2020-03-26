function [clu_X_mat_re,clu_Y_mat_re,clu_Z_mat_re,obj_mat_re,JSD_mat_re, clu_X_bes,clu_Y_bes, clu_Z_bes] = coupleCoC(X, Y, nclu_X, nclu_Y, nclu_Z, epochs, niter, lambda, alpha,iter)
%%initialize the input for niter_outer times
clu_X_mat = zeros(epochs,size(X,1));
clu_Y_mat = zeros(epochs,size(Y,1));
clu_Z_mat = zeros(epochs,size(X,2));
for l = 1:epochs
clu_X_mat(l,:) = randsample(nclu_X, size(X,1),true);
clu_Y_mat(l,:) = randsample(nclu_Y, size(Y,1),true);
clu_Z_mat(l,:) = randsample(nclu_Z, size(X,2),true);
end

%%generate empty lists in order to save results obtained by each outer
% For each matrix inside the lists, the first row is the initialization of
% clustering;
clu_X_mat_re = zeros(niter+1,size(X,1),epochs);
clu_Y_mat_re = zeros(niter+1,size(Y,1),epochs);
clu_Z_mat_re = zeros(niter+1,size(X,2),epochs);
obj_mat_re = zeros(niter,epochs);
JSD_mat_re = zeros(niter,epochs);

%%outer iteration
title1 = sprintf('coupleCoC, the NO.%d iteration in progress',iter);
title2 = sprintf('coupleCoC, remaining time for the NO.%d iteration =',iter);
h = waitbar(0,title1);
s = clock;
for i = 1:epochs
    %begin process
    [clu_X_Trace, clu_Y_Trace, clu_Z_Trace, obj_val_Trace, JSD_Trace] = coupleCoC_fun(X, Y, clu_X_mat(i,:), clu_Y_mat(i,:), clu_Z_mat(i,:), nclu_X, nclu_Y, nclu_Z, niter, lambda, alpha);
    clu_X_mat_re(:,:,i) = clu_X_Trace;
    clu_Y_mat_re(:,:,i) = clu_Y_Trace;
    clu_Z_mat_re(:,:,i) = clu_Z_Trace;
    obj_mat_re(:,i) = obj_val_Trace;
    JSD_mat_re(:,i) = JSD_Trace;
    %end process
    %begin estimate remaining time
    if i == 1
        is = etime(clock,s);
        esttime = is*epochs;
    end
    h = waitbar(i/epochs,h,...
    [title2,num2str(esttime-etime(clock,s),'%4.1f'),'sec' ]);
    %end estimate remaining time
end
close(h);

% choose the best result and return
obj = obj_mat_re(niter,:);
bes = find(obj == min(obj), 1 );
clu_X_bes = clu_X_mat_re(niter+1,:,bes);
clu_Y_bes = clu_Y_mat_re(niter+1,:,bes);
clu_Z_bes = clu_Z_mat_re(niter+1,:,bes);


