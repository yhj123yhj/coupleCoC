function p_tilde_XZ = cal_coclu_prob(p_XZ, p_XZ_tilde, clu_X, clu_Z)
p_tilde_XZ = zeros(size(p_XZ,1),size(p_XZ,2));
for i = 1:size(p_XZ,1)
    for j = 1:size(p_XZ,2)
        p_tilde_XZ(i,j) = (p_XZ_tilde(clu_X(i),clu_Z(j))*sum(p_XZ(i,:))*sum(p_XZ(:,j)))/(sum(p_XZ_tilde(clu_X(i),:))*sum(p_XZ_tilde(:,clu_Z(j))));
    end
end
