function ReducedModel = InvariantManifold(Model,style)
% Use the DPIM to calculate the invariant manifold and corresponding
% reduced-order model

CS = Model.Indp.type;
n_p = length(Model.Indp.I);
Master = Model.Master;

M = Model.M;
C = Model.C;
K = Model.K;
Lambda = Model.spectrum.Lambda;

nonlinear_force_type = Model.nonlinear_force_type;
if strcmpi(nonlinear_force_type, 'multi_index')
    Fn = cell(length(Model.nonlinear_force),1);
    for ni = 1:length(Model.nonlinear_force)
        if isempty(Model.nonlinear_force{ni})
            Fn{ni} = [];
        else
            Fn{ni} = change_tensor_Mu(Model.nonlinear_force{ni});
            Fn{ni}.I = change_Mu_Ind(Fn{ni}.I);
        end
    end
else
    Fn = Model.nonlinear_force;
end

if nargin <= 1
    MS = 'zero';
    CM = 'Step';
else
    MS = style.MS;
    CM = style.CM;
end

f = cell(n_p,1);
Psi = cell(n_p,1);
Upsilon = cell(n_p,1);
Indp = Model.Indp;
N = Model.N;
n = length(Master);
X = Model.spectrum.X;
Y = Model.spectrum.Y;

% Reduced-order model and invariant manifold of order p = 1
f{1} = Model.spectrum.L(Master,Master);
Psi{1} = Model.spectrum.Y(1:N,Master);
Upsilon{1} = Model.spectrum.Y(N+1:end,Master);

for p_i = 2:n_p
% Reduced-order model and invariant manifold of order p
I_pi = Indp.I{p_i};
ResonanceFlag_pi = Indp.ResonanceFlag{p_i};
ResonanceI_pi = Indp.ResonanceI{p_i};

n_Ipi = size(I_pi,2);
f{p_i} = zeros(n, n_Ipi);
Psi{p_i} = zeros(N, n_Ipi);
Upsilon{p_i} = zeros(N, n_Ipi);

switch CM
    case 'Step'
        for l_i = 1:n_Ipi
            I = I_pi(:, l_i);
            switch CS
                case 'index'
                    theta_I = change_Ind_Mu(I, n)'*Lambda(Master);
                case 'multi_index'
                    theta_I = I'*Lambda(Master);
            end
            mu_I = getRInd(Psi, f, Indp, I, CS);
            nu_I = getRInd(Upsilon, f, Indp, I, CS);
            G_I = getFnI(Fn, Psi, Indp, I, 2, nonlinear_force_type, CS);
            H_I = getFnI(Fn, Psi, Indp, I, 3, nonlinear_force_type, CS);
            Xi_I = - G_I - H_I - M*nu_I - (theta_I*M + C)*mu_I;
            switch MS
                case 'zero'
                    if ResonanceFlag_pi(l_i) == 1
                        % Exist resonance
                        IR = ResonanceI_pi{l_i};
                        Ur = Y(1:N,Master(IR)); Vr = Y(N+1:end,Master(IR));
                        Xr1 = X(1:N,Master(IR)); Xr2 = X(N+1:end,Master(IR));
                        A11 = theta_I^2*M + theta_I*C + K;
                        A12 = (theta_I*M + C)*Ur+M*Vr;
                        A21 = Xr1'*(C+theta_I*M) + Xr2'*M;
                        A22 = Xr1'*M*Ur;
                        A = [A11 A12;A21 A22];
                        B = [Xi_I; - Xr1'*M*mu_I];
                        Phf = A\B;
                        Psi{p_i}(:,l_i) = Phf(1:N);
                        f{p_i}(IR,l_i) = Phf(N+1:end);
                        Upsilon{p_i}(:,l_i) = theta_I*Phf(1:N) + Ur*Phf(N+1:end) + mu_I;
                    else
                        % Non resonance
                        Psi{p_i}(:,l_i) = (theta_I^2*M + theta_I*C + K)\Xi_I;
                        Upsilon{p_i}(:,l_i) = theta_I*Psi{p_i}(:,l_i) + mu_I;
                    end
                case 'lsqminnorm'
                    if ResonanceFlag_pi(l_i) == 1
                        % Exist resonance
                        IR = ResonanceI_pi{l_i};
                        Ur = Y(1:N,Master(IR)); Vr = Y(N+1:end,Master(IR));
                        Xr = X(:,Master(IR));
                        f{p_i}(IR,l_i) = Xr'*[-M*nu_I-C*mu_I-G_I-H_I;-M*mu_I];
                        A11 = theta_I^2*M + theta_I*C + K;
                        A12 = (theta_I*M + C)*Ur+M*Vr;
                        A = A11;
                        B = Xi_I - A12*f{p_i}(IR,l_i);
                        Phf = lsqminnorm(A,B);
                        Psi{p_i}(:,l_i) = Phf;
                        Upsilon{p_i}(:,l_i) = theta_I*Phf + Ur*f{p_i}(IR,l_i) + mu_I;
                    else
                        % Non resonance
                        Psi{p_i}(:,l_i) = (theta_I^2*M + theta_I*C + K)\Xi_I;
                        Upsilon{p_i}(:,l_i) = theta_I*Psi{p_i}(:,l_i) + mu_I;
                    end
            end
        end
    case 'Total'
        mu_I = zeros(N, n_Ipi);
        nu_I = mu_I;
        G_I = mu_I;
        H_I = mu_I;
        theta_I = zeros(n_Ipi,1);
        for l_i = 1:n_Ipi
            I = I_pi(:, l_i); 
            switch CS
                case 'index'
                    theta_I(l_i) = change_Ind_Mu(I, n)'*Lambda(Master);
                case 'multi_index'
                    theta_I(l_i) = I'*Lambda(Master);
            end
            mu_I(:,l_i) = getRInd(Psi, f, Indp, I, CS);
            nu_I(:,l_i) = getRInd(Upsilon, f, Indp, I, CS);
            G_I(:,l_i) = getFnI(Fn, Psi, Indp, I, 2, nonlinear_force_type, CS);
            H_I(:,l_i) = getFnI(Fn, Psi, Indp, I, 3, nonlinear_force_type, CS);

            if ResonanceFlag_pi(l_i) == 1
                % Exist resonance
                IR = ResonanceI_pi{l_i};
                Xr = X(:,Master(IR));
                f{p_i}(IR,l_i) = Xr'*[-M*nu_I(:,l_i)-C*mu_I(:,l_i)-G_I(:,l_i)-H_I(:,l_i);-M*mu_I(:,l_i)];
            end
        end
        A = kron(diag(theta_I),[C M;M zeros(N)]) + kron(eye(n_Ipi),blkdiag(K,-M));
        B = -kron(eye(n_Ipi),[C M;M zeros(N)]*Y(:,Master))*reshape(f{p_i},[],1)...
            -reshape([M*nu_I+C*mu_I+G_I+H_I;M*mu_I],[],1);
        Mf = reshape(lsqminnorm(A,B),2*N,[]);
        Psi{p_i} = Mf(1:N,:);
        Upsilon{p_i} = Mf(N+1:end,:);
end

end
ReducedModel.f = f;
ReducedModel.Psi = Psi;
ReducedModel.Upsilon = Upsilon;
ReducedModel.I = Indp.I;

if strcmpi(CS, 'index')
    ReducedModel = CompactModel(ReducedModel);
end

end