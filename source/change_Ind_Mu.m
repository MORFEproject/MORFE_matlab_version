function Mu = change_Ind_Mu(Ind,Nmax)
Mu = zeros(Nmax,size(Ind,2));
for ni = 1:Nmax
    for nj = 1:size(Ind,2)
        Mu(ni,nj) = sum(Ind(:,nj) == ni);
    end
end
end