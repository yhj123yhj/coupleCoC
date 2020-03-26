function [clu_X_Trace, clu_Y_Trace, clu_Z_Trace, obj_val_Trace, JSD_Trace] = coupleCoC_fun(X, Y, clu_X, clu_Y, clu_Z, nclu_X, nclu_Y, nclu_Z, niter, lambda, alpha)
% Corresponding to the objective function: J =
% I(X,Z)-I(X_tilde,Z)+lambda(2)*(I(Y,Z)-I(Y_tilde,Z))+lambda(3)*JSD(p(X_tilde)||p(Y_tilde))

p_XZ = X/sum(sum(X));
q_YZ = Y/sum(sum(Y));
p_XZ_tilde = cal_clu_prob(p_XZ,clu_X,clu_Z,nclu_X, nclu_Z);
q_YZ_tilde = cal_clu_prob(q_YZ,clu_Y, clu_Z, nclu_Y, nclu_Z);

clu_X_Trace = zeros(niter+1, size(clu_X,2));
clu_Y_Trace = zeros(niter+1, size(clu_Y,2));
clu_Z_Trace = zeros(niter+1, size(clu_Z,2));

obj_val_Trace = zeros(niter,1);
JSD_Trace = zeros(niter,1);
p_XZ_tilde_Trace = zeros(size(p_XZ_tilde,1),size(p_XZ_tilde,2), niter+1);
q_YZ_tilde_Trace = zeros(size(q_YZ_tilde,1),size(q_YZ_tilde,2), niter+1);

clu_X_Trace(1,:) = clu_X;
clu_Y_Trace(1,:) = clu_Y;
clu_Z_Trace(1,:) = clu_Z;

p_XZ_tilde_Trace(:,:,1) = p_XZ_tilde;
q_YZ_tilde_Trace(:,:,1) = q_YZ_tilde;

for iter = 1:niter
%Update X
[clu_X_Trace(iter+1,:), ~] = update_clu_jsd(p_XZ, clu_X_Trace(iter,:), clu_Z_Trace(iter,:), q_YZ_tilde_Trace(:,:,iter),nclu_X, nclu_Z, alpha);

%Update Y
[clu_Y_Trace(iter+1,:),obj_Y] = update_clu_jsd(q_YZ, clu_Y_Trace(iter,:), clu_Z_Trace(iter,:), p_XZ_tilde_Trace(:,:,iter), nclu_Y, nclu_Z, alpha);

%Update Z
clu_Z_Trace(iter+1,:) = update_clu1(p_XZ, q_YZ, lambda, clu_X_Trace(iter+1,:), clu_Y_Trace(iter+1,:), clu_Z_Trace(iter,:), nclu_X, nclu_Y, nclu_Z);

%Update p_XZ_tilde, q_YZ_tilde
p_XZ_tilde_Trace(:,:,iter+1) = cal_clu_prob(p_XZ, clu_X_Trace(iter+1,:), clu_Z_Trace(iter+1,:), nclu_X, nclu_Z);
q_YZ_tilde_Trace(:,:,iter+1) = cal_clu_prob(q_YZ, clu_Y_Trace(iter+1,:), clu_Z_Trace(iter+1,:), nclu_Y, nclu_Z);

%Label switch
[CX_temp_best, obj_X2, ~] = label_switch(p_XZ, q_YZ_tilde_Trace(:,:,iter+1), clu_X_Trace(iter+1,:), clu_Z_Trace(iter+1,:), nclu_X, nclu_Z);
clu_X_Trace(iter+1,:) = CX_temp_best;
obj_val_Trace(iter) = obj_X2 + lambda(2)*obj_Y;
end

