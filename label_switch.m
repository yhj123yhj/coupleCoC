function [CX_temp_best,obj,JSD_min] = label_switch(p_XZ, q_YZ_tilde_temp0, clu_X, clu_Z, nclu_X, nclu_Z)
per_X = perms(1:nclu_X);
CX_temp = clu_X;
CX_temp_mat = zeros(size(per_X,1), size(CX_temp,2));
pX_tilde_mat = zeros(nclu_X, size(per_X,1));
for j = 1:size(per_X,1)
    CX_temp0 = CX_temp;
    for i = 1:nclu_X
        CX_temp0(CX_temp==i) = per_X(j,i);
    end
    CX_temp_mat(j,:) = CX_temp0;
    p_XZ_tilde_temp00 = cal_clu_prob(p_XZ, CX_temp0, clu_Z, nclu_X, nclu_Z);
    pX_tilde_mat(:,j) = sum(p_XZ_tilde_temp00,2);
end
    
pX_labswi = pX_tilde_mat;
pY_labswi = sum(q_YZ_tilde_temp0,2);
pXY_mean = (pX_labswi + pY_labswi)/2;
JSD = 0.5*sum(pX_labswi.*log(pX_labswi./pXY_mean),1) + 0.5*sum(pY_labswi.*log(pY_labswi./pXY_mean),1);
label_best = find(JSD==min(JSD),1);
JSD_min = min(JSD);
CX_temp_best = CX_temp_mat(label_best,:);

p_XZ_tilde_temp0 = cal_clu_prob(p_XZ, CX_temp_best, clu_Z, nclu_X, nclu_Z);
p_tilde_XZ_temp0 = cal_coclu_prob(p_XZ, p_XZ_tilde_temp0, CX_temp_best, clu_Z);
obj = 0;
for j = 1:size(p_XZ,1)
    p_post_zx = p_XZ(j,:)/sum(p_XZ(j,:));
    p_tilde_zx_tilde = p_tilde_XZ_temp0(j,:)/sum(p_XZ(j,:));
    ind = p_XZ(j,:)~=0;
    obj = obj + sum(p_post_zx(ind).*log(p_post_zx(ind)./(p_tilde_zx_tilde(ind))));
end
    