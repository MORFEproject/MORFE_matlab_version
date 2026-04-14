function Mu = change_tensor_Mu(tensor)
ind_tensor = find(tensor ~= 0);

l_tensor = size(ind_tensor, 1);

N1 = size(tensor,1);
Nmax = size(tensor,2);

nm = 1;
Mu.I = [];
Mu.vector = [];
for ni = 1:l_tensor
    Multi_index = change_Ind_Mu(ind_tensor(ni, 2:end)', Nmax);
    if ni == 1
        jni = 0;
    else
        jni = find_index_VinA(Multi_index, Mu.I);
    end
    
    if jni > 0
        Mu.vector(:, jni) = Mu.vector(:, jni)...
            + [zeros(ind_tensor(ni,1)-1,1);1;zeros(N1-ind_tensor(ni,1),1)]*tensor(ind_tensor(ni,:));
    else
        Mu.I(:,nm) = Multi_index;
        Mu.vector(:,nm) = [zeros(ind_tensor(ni,1)-1,1);1;zeros(N1-ind_tensor(ni,1),1)]*tensor(ind_tensor(ni,:));
        nm = nm + 1;
    end
end
end

function Mu = change_Ind_Mu(Ind,Nmax)
Mu = zeros(Nmax,size(Ind,2));
for ni = 1:Nmax
    for nj = 1:size(Ind,2)
        Mu(ni,nj) = sum(Ind(:,nj) == ni);
    end
end
end

function index = find_index_VinA(V, A)
AminusV_norm = sum(abs(A-V),1);
index = find(AminusV_norm == 0);
end