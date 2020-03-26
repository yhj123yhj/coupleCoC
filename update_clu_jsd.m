function [clu_X_temp0, obj] = update_clu_jsd(p_XZ0, clu_X0, clu_Z0, q_YZ_tilde_temp0,nclu_X0, nclu_Z0, alpha)
clu_X_temp0 = clu_X0;
for i = 1:size(p_XZ0,1)
    KL_X = zeros(nclu_X0,1);
    JSD_X = zeros(nclu_X0,1);
    for k = 1:nclu_X0
        clu_X_temp0(i) = k;
        p_XZ_tilde_temp0 = cal_clu_prob(p_XZ0, clu_X_temp0, clu_Z0, nclu_X0, nclu_Z0);
        p_tilde_XZ_temp0 = cal_coclu_prob(p_XZ0, p_XZ_tilde_temp0, clu_X_temp0, clu_Z0);
        p_post_zx = p_XZ0(i,:)/sum(p_XZ0(i,:));
        p_tilde_zx_tilde = p_tilde_XZ_temp0(i,:)/sum(p_XZ0(i,:));
        ind = p_XZ0(i,:)~=0;
        % add the Jesen-Shannon Divergence
        pX_labswi = sum(p_XZ_tilde_temp0,2);
        pY_labswi = sum(q_YZ_tilde_temp0,2);
        pXY_mean = (pX_labswi + pY_labswi)/2;
        JSD_X(k) = 0.5*sum(pX_labswi.*log(pX_labswi./pXY_mean),1) + 0.5*sum(pY_labswi.*log(pY_labswi./pXY_mean),1);
        KL_X(k) = sum(p_post_zx(ind).*log(p_post_zx(ind)./(p_tilde_zx_tilde(ind)))) + alpha*JSD_X(k);
    end
    clu_X_temp0(i) = find(KL_X == min(KL_X), 1 );
end

p_XZ_tilde_temp0 = cal_clu_prob(p_XZ0, clu_X_temp0, clu_Z0, nclu_X0, nclu_Z0);
p_tilde_XZ_temp0 = cal_coclu_prob(p_XZ0, p_XZ_tilde_temp0, clu_X_temp0, clu_Z0);
obj = 0;
for j = 1:size(p_XZ0,1)
    p_post_zx = p_XZ0(j,:)/sum(p_XZ0(j,:));
    p_tilde_zx_tilde = p_tilde_XZ_temp0(j,:)/sum(p_XZ0(j,:));
    ind = p_XZ0(j,:)~=0;
    obj = obj + sum(p_post_zx(ind).*log(p_post_zx(ind)./(p_tilde_zx_tilde(ind))));
end
% add the Jesen-Shannon Divergence
pX_labswi = sum(p_XZ_tilde_temp0,2);
pY_labswi = sum(q_YZ_tilde_temp0,2);
pXY_mean = (pX_labswi + pY_labswi)/2;
JSD = 0.5*sum(pX_labswi.*log(pX_labswi./pXY_mean),1) + 0.5*sum(pY_labswi.*log(pY_labswi./pXY_mean),1);
obj = obj + alpha*JSD;


