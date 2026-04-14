function CRM = CompactModel(RM)

RMf = RM.f;
RMPsi = RM.Psi;
RMUpsilon = RM.Upsilon;
RMI = RM.I;

n = size(RMf{1},1);
N = size(RMPsi{1},1);
np = length(RMI);

Mu = cell(np,1);
Cf = cell(np,1);
CPsi = cell(np,1);
CUpsilon = cell(np,1);

for p_i = 1:np
    Mu{p_i} = flip(findAllVectors(n,p_i),2);
    l_I = size(RMI{p_i},2);
    l_mu = size(Mu{p_i},2);
    Cf{p_i} = zeros(n, l_mu);
    CPsi{p_i} = zeros(N, l_mu);
    CUpsilon{p_i} = zeros(N, l_mu);

    for l_i = 1:l_I
        MuI = change_Ind_Mu(RMI{p_i}(:,l_i), n);
        Index = find_index_VinA(MuI,Mu{p_i});
        Cf{p_i}(:,Index) = Cf{p_i}(:,Index) + RMf{p_i}(:,l_i);
        CPsi{p_i}(:,Index) = CPsi{p_i}(:,Index) + RMPsi{p_i}(:,l_i);
        CUpsilon{p_i}(:,Index) = CUpsilon{p_i}(:,Index) + RMUpsilon{p_i}(:,l_i);
    end
end
CRM.f = Cf;
CRM.Psi = CPsi;
CRM.Upsilon = CUpsilon;
CRM.I = Mu;
end