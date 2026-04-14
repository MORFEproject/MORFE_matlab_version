function Indp = change_Mu_Ind(Mu)
% 将多指标转化为指标索引

[n1,n2] = size(Mu);
Indp = zeros(sum(Mu(:,1)),n2);
for ni = 1:n2
    for nj = 1:n1
        Indp(1+sum(Mu(1:nj-1,ni)):sum(Mu(1:nj-1,ni))+Mu(nj,ni),ni) = nj;
    end
end
end
