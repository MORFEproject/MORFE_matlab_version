function RealModel = getRealModel(Compact_model,T)
f = Compact_model.f;
Psi = Compact_model.Psi;
Upsilon = Compact_model.Upsilon;
I = Compact_model.I;

f_new = cell(length(I),1);
Psi_new = f_new;
Upsilon_new = f_new;

for ni = 1:length(I)
    n_z = size(T,2); 

    Muni = flip(findAllVectors(n_z,ni),2); 

    n_mu = size(Muni,2); 
    n_f = size(f{ni},1); n_Psi = size(Psi{ni},1); n_Upsilon = size(Upsilon{ni},1);
    f_ni = zeros(n_f, n_mu);
    Psi_ni = zeros(n_Psi, n_mu);
    Upsilon_ni = zeros(n_Upsilon, n_mu);

    for nj = 1:n_mu
        f_ninj = zeros(n_f, 1);
        Psi_ninj = zeros(n_Psi, 1);
        Upsilon_ninj = zeros(n_Upsilon, 1);
        
        Mu_new = generateCombinations(ones(1,ni), Muni(:,nj));

        Ind_pre = change_Mu_Ind(I{ni});
        n_pre = size(Ind_pre,2);
        for nk = 1:size(Mu_new, 3)
            Ind_new = change_Mu_Ind(Mu_new(:,:,nk));
            Ind_linear = sub2ind(size(T), reshape(Ind_pre,[],1),...
                repmat(Ind_new',size(Ind_pre,2),1));
            f_ninj = f_ninj + sum(f{ni}.*prod(T(reshape(Ind_linear,[],n_pre)),1),2);
            Psi_ninj = Psi_ninj + sum(Psi{ni}.*prod(T(reshape(Ind_linear,[],n_pre)),1),2);
            Upsilon_ninj = Upsilon_ninj + sum(Upsilon{ni}.*prod(T(reshape(Ind_linear,[],n_pre)),1),2);
        end

        f_ni(:, nj) = f_ninj;
        Psi_ni(:, nj) = Psi_ninj;
        Upsilon_ni(:, nj) = Upsilon_ninj;
    end
    
    f_new{ni}.I = Muni; f_new{ni}.vector = real(inv(T)*f_ni);
    Psi_new{ni}.I = Muni; Psi_new{ni}.vector = real(Psi_ni);
    Upsilon_new{ni}.I = Muni; Upsilon_new{ni}.vector = real(Upsilon_ni);
end

RealModel.f = f_new;
RealModel.Psi = Psi_new;
RealModel.Upsilon = Upsilon_new;
end