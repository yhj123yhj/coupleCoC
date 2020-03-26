function clu_Z_temp = update_clu1(p_XZ, q_YZ, lambda, clu_X_temp, clu_Y_temp, clu_Z, nclu_X, nclu_Y,nclu_Z)
clu_Z_temp = clu_Z; % to store the updated value of the target cluster for each x in X;
for i = 1:size(p_XZ,2)
    KL_Z = zeros(nclu_Z,1);
    for k = 1:nclu_Z
        clu_Z_temp(i) = k;
        p_XZ_tilde_temp = cal_clu_prob(p_XZ, clu_X_temp, clu_Z_temp, nclu_X, nclu_Z);
        p_tilde_XZ_temp = cal_coclu_prob(p_XZ, p_XZ_tilde_temp, clu_X_temp, clu_Z_temp);
        p_post_XZ = p_XZ(:,i)/sum(p_XZ(:,i));
        p_tilde_XZ_tilde = p_tilde_XZ_temp(:,i)/sum(p_XZ(:,i));
        
        p_YZ_tilde_temp = cal_clu_prob(q_YZ, clu_Y_temp, clu_Z_temp, nclu_Y, nclu_Z);
        p_tilde_YZ_temp = cal_coclu_prob(q_YZ, p_YZ_tilde_temp, clu_Y_temp, clu_Z_temp);
        p_post_YZ = q_YZ(:,i)/sum(q_YZ(:,i));
        p_tilde_YZ_tilde = p_tilde_YZ_temp(:,i)/sum(q_YZ(:,i));
        
        ind1 = p_XZ(:,i)~=0;
        ind2 = q_YZ(:,i)~=0;
        temp1 = sum(p_post_XZ(ind1).*log(p_post_XZ(ind1)./(p_tilde_XZ_tilde(ind1))));   
        temp2 = sum(p_post_YZ(ind2).*log(p_post_YZ(ind2)./(p_tilde_YZ_tilde(ind2))));
        KL_Z(k) = lambda(1)*temp1*sum(p_XZ(:,i)) + lambda(2)*temp2*sum(q_YZ(:,i));
    end
    clu_Z_temp(i) = find(KL_Z == min(KL_Z), 1 );
end