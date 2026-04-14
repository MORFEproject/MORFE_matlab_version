function Model = spectrum_analysis(Model,n,type)
% 对Model进行特征谱分析
M = Model.M;
C = Model.C;
K = Model.K;
N = size(M,2);

A = [C M;M zeros(size(M))];
B = blkdiag(-K,M);

if nargin <= 2
    type = 'RealMode';
else
    type = type.MAS;
end

switch type
    case 'RealMode'
        [PhTol,Lamb] = eigs(K,M,n,'smallestabs');
        omega = sqrt(diag(Lamb));
        [~, so] = sort(real(omega),1,"ascend");
        Ph = PhTol(:,so);
        Ph = Ph/sqrt(diag(diag(Ph'*M*Ph)));
        Lambdaj = (-diag(Ph'*C*Ph) + sqrt(diag(Ph'*C*Ph).^2 - 4*diag(Ph'*K*Ph)))/2;
        LambdajN = (-diag(Ph'*C*Ph) - sqrt(diag(Ph'*C*Ph).^2 - 4*diag(Ph'*K*Ph)))/2;
        
        L = diag([Lambdaj; LambdajN]);
        Yj = [Ph;Ph*diag(Lambdaj)];
        YjN = [Ph;Ph*diag(LambdajN)];
        Y = [Yj YjN];
        X = conj(Y);
        X = X*diag(1./diag(X'*A*Y)');

    case 'ComplexMode'
        [PhTol,Lamb] = eigs(B,A,2*n,'smallestabs');
        
        omega = diag(Lamb);
        [~,so] = sort(imag(omega),'ascend');
        Omega = omega([so(n+1:end);flip(so(1:n))]);
        Ph = PhTol(:,[so(n+1:end);flip(so(1:n))]);
        for ni = 1:size(Ph,2)
            Ph(:,ni) = Ph(:,ni)/sqrt(Ph(1:N,ni)'*M*Ph(1:N,ni));
        end
        
        Y = Ph;
        [PhToll,Lambl] = eigs(B',A',2*n,'smallestabs');
        omegal = diag(Lambl)';
        so = zeros(length(omegal),1);
        for ni = 1:length(Omega)
            so(ni) = find(abs(omegal - Omega(ni))<1e-3*abs(Omega(ni)));
        end
        Phl = PhToll(:,so);
        X = Phl/diag(diag(Phl'*A*Y)');
        L = diag(diag(X'*B*Y));
end

spectrum = struct('A',A,'B',B,'L',L,'X',X,'Y',Y,'Lambda',diag(L));
Model.spectrum = spectrum;
end