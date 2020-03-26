function p_XZ_tilde = cal_clu_prob(p_XZ,clu_X,clu_Z,nclu_X, nclu_Z)
%Initialize joint probability distribution of co-clustering, p_tilde_XZ and p_tilde_YZ_tilde
p_XZ_tilde = zeros(nclu_X, nclu_Z);
for i = 1:nclu_X
    for j = 1:nclu_Z
        p_XZ_tilde(i,j) = sum(sum(p_XZ(clu_X==i,clu_Z==j)));
    end
end
