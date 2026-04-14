function [R,dR,Q] = ...
    HB_residual_mix(X,system,H,N,analysis_type,varargin)
%% Handle input variables depending on the modus

type = system.type;

% System matrices
if strcmpi(type,'first order')
    A = system.A; B = system.B; 
    M = zeros(size(A)); D = A; K = -B; 
    V = system.Ve; n = size(A,1); % Dimension of second-order model
else
    M = system.M; D = system.D; K = system.K;
    V = system.Ve; n = system.ndof; % Dimension of first-order model
end
    

% Conversion from sine-cosine to complex-exponential
% representation, where X consists of real numbers and Q consists of complex numbers
I0 = 1:n; ID = n+(1:H*n);
IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); IS = IC+n;
dX = eye(length(X));
Q = zeros(n*(H+1),1);   dQ = zeros(size(Q,1),size(dX,2));
Q(I0) = X(I0);          dQ(I0,:) = dX(I0,:);
Q(ID) = X(IC)-1i*X(IS); dQ(ID,:) = dX(IC,:)-1i*dX(IS,:);

% Handle analysis type
if nargin<=4 || isempty(analysis_type)
    % Default analysis: frequency response
    analysis_type = 'frf';
end
switch lower(analysis_type)
    case {'frf','frequency response'}
        % Frequency response analysis: X = [Q; Om], where Q consists of harmonic coefficients 
        % of the first H orders, arranged in the order of cosine terms followed by sine terms.
        if isfield(system, 'fextm')
            fextm = system.fextm;
        else
            fextm = 1; % Excitation
        end
        Fex = zeros(n*(H+1), 1);
        h = 1; % Force acts on the h-th harmonic component.
        Fex(h*n + (1:n)) = fextm*V;
        dFex = zeros(size(Fex,1), length(X));

        % Derivative of damping (does not depend on unknowns in the case of
        % frequency response)
        dD_dalpha = 0*M; dalpha = zeros(1,length(X));

        % Excitation frequency
        Om = X(end);
        dOm = dX(end,:);

        % Scaling of dynamic force equilibrium
        if length(varargin)<2 || isempty(varargin{2})
            fscl = 1;
        else
            fscl = varargin{2};
        end
    case {'nma','nonlinear modal analysis'}
        % Nonlinear modal analysis:  X = [Psi;Om;del;log10a]
        
        % Interpret additional input
        inorm = varargin{1};
        
        % Modal mass
        a = exp(log(10)*X(end));
        da = log(10)*exp(log(10)*X(end))*dX(end,:);
        
        % In this case, the harmonics are amplitude-normalized. We thus
        % have to scale them by the amplitude a.
        Psi = Q; dPsi = dQ;
        Q = Psi*a;
        dQ = dPsi*a + Psi*da;
        
        % Modal frequency
        Om = X(end-2);
        dOm = dX(end-2,:);
        
        % Modal damping ratio
        del = X(end-1);
        ddel = dX(end-1,:);
        
        % Extended Periodic Motion concept: artifical negative viscous 
        % mass-proportional damping
        alpha = 2*Om*del;
        dalpha = 2*dOm*del+2*Om*ddel;
        D = D-alpha*M;
        dD_dalpha = -M;
        
        % No external forcing
        Fex = zeros(n*(H+1),1);
        dFex = zeros(size(Fex,1),length(X));
        
        % Scaling of dynamic force equilibrium
        if length(varargin)<2 || isempty(varargin{2})
            fscl = 1;
        else
            fscl = varargin{2};
        end
    otherwise
        error(['Unknown analysis type ' analysis.type '.']);
end
%% Computation of the Fourier coefficients of the nonlinear forces and the 
% Jacobian using AFT
[Fnl,dFnl] = HB_nonlinear_forces_AFT(Q,dQ,Om,dOm,H,N,...
    system.nonlinear_elements,type,n);

%% Assembly of the residual and the Jacobian

% Dynamic force equilibrium

Rc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*Q + Fnl - Fex;
dRc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*dQ + dFnl - dFex  + ...
    ( -2*Om*kron(diag((0:H).^2),M) + 1i*kron(diag(0:H),D) )*Q*dOm + ...
    1i*Om*kron(diag(0:H),dD_dalpha)*Q*dalpha;

% Scale dynamic force equilibrium (useful for numerical reasons)
Rc = 1/fscl*(Rc);
dRc = 1/fscl*(dRc);

% Conversion from complex-exponential to sine-cosine representation
R = zeros(size(X,1)-1,1); dR = zeros(size(X,1)-1,size(X,1));
R(I0) = real(Rc(I0)); dR(I0,:) = real(dRc(I0,:));
R(IC) = real(Rc(ID)); dR(IC,:) = real(dRc(ID,:));
R(IS) = -imag(Rc(ID)); dR(IS,:) = -imag(dRc(ID,:));

if strcmpi(analysis_type,'nma') || ...
        strcmpi(analysis_type,'nonlinear modal analysis')
    % Scale dynamic force equilibrium by modal amplitude
    % NOTE: We first evaluate the derivative, as we then overwrite R!
    dR(1:end-2,:) = dR(1:end-2,:)/a-R(1:end-2)/a^2*da;
    R(1:end-2) = R(1:end-2)/a;
    
    % Amplitude normalization: The mass of the nonlinear mode shape (all
    % harmonics) is enforced to be one.
    R(end-1) = real(Psi'*kron(eye(H+1),M)*Psi-1);
    dR(end-1,:) = real(2*(Psi'*kron(eye(H+1),M))*dPsi);
    
    % Phase normalization: Velocity of coordinate 'inorm' is enforced to be 
    % zero at t=0.
    R(end) = (1:H)*imag(Psi(inorm+(n:n:H*n)));
    dR(end,:) = (1:H)*imag(dPsi(inorm+(n:n:H*n),:));
end
end
%% Computation of the Fourier coefficients of the nonlinear forces and the 
% Jacobian using AFT
function [F,dF] = ...
    HB_nonlinear_forces_AFT(Q,dQ,Om,dOm,H,N,nonlinear_elements,type,n)
%% Initialize output
if strcmpi(type,'first order')
    F = zeros(n*(H+1),1);
    dF = zeros(n*(H+1),size(dQ,2));
    for nl = 1:length(nonlinear_elements)
        % Sampling frequency
        tau = (0:2*pi/N:2*pi-2*pi/N)';
    
        if nonlinear_elements{nl}.islocal
            %% Calculate the nonlinear restoring force in the time domain.
            switch lower(nonlinear_elements{nl}.type)
                case 'custom'
                    % add local nonlinearity
            end
    
        else % Global nonlinearity
            switch lower(nonlinear_elements{nl}.type)
                case 'polynomialstiffness'
                    % Inverse FFT
                    Qc = transpose(reshape(Q,[],H+1)); n0 = size(Qc,2);

                    if isfield(nonlinear_elements{nl}, 'Ind')
                        Indx = nonlinear_elements{nl}.Ind;
                        nx = length(Indx);
                    else
                        Indx = 1:n0;
                        nx = n0;
                    end

                    Qc = [Qc(1,Indx); Qc(2:end,Indx)/2; zeros(N-H-1,nx)];
                    Qc(end-H+1:end,:) = flipud(conj(Qc(2:H+1,:)));
                    q = ifft(Qc)*N;

                    % Exponents and coefficients
                    pp = nonlinear_elements{nl}.exponents;
                    Et = -nonlinear_elements{nl}.coefficients;
                    
                    % Evaluate polynimials
                    z = prod(kron(q,ones(size(pp,1),1)) .^ repmat(pp,N,1),2);
                    
                    % Evaluate forces
                    nz = size(Et,2);
                    fnl = (Et*reshape(z,nz,N))';
                    
                    % Forward FFT
                    Fnl_ind = fft(fnl)/N;
                    Fnl_ind = [Fnl_ind(1,:);Fnl_ind(2:H+1,:)*2];
                    
                    Fnl = zeros(H+1,n0);
                    Fnl(:,Indx) = Fnl_ind;

                    % Add forces to global force vector
                    F = F + reshape(transpose(Fnl),[],1);
                    %% Analytical gradients 
                    % (If the double loop bothers you, vectorize!)
                     for l=1:nx
                        % Apply IDFT to Jacobian of Q
                        ndx = size(dQ,2);
                        dQlc = [dQ(Indx(l),:); dQ(n0+Indx(l):n0:end,:)/2; zeros(N-H-1,ndx)];
                        dQlc(end-H+1:end,:) = flipud(conj(dQlc(2:H+1,:)));
                        dql = ifft(dQlc)*N;
                        
                        % Derive
                        notl = setdiff(1:nx,l);
                        dzql_dql = repmat(pp(:,l),N,1).*...
                            (kron(q(:,Indx(l)),ones(size(pp,1),1)).^repmat(pp(:,l)-1,N,1));
                        dzql_dql(isnan(dzql_dql)) = 0;
                        dz_dql = prod([ dzql_dql ...
                            kron(q(:,Indx(notl)),ones(size(pp,1),1)).^...
                            repmat(pp(:,notl),N,1)], 2);
                        df_dql = (Et*reshape(dz_dql,nz,N))';
                        
                        % Derive
                        for i=1:nx
                            dfi = repmat(df_dql(:,i),1,size(dql,2)).*dql;
                            dFi = fft(dfi(end-N+1:end,:))/N;
                            dFi = [dFi(1,:);dFi(2:H+1,:)*2];
                            dF(Indx(i):n0:end,:) = dF(Indx(i):n0:end,:) + dFi;
                        end
                        
                    end
                otherwise
                    error(['Unknown global nonlinearity type ' ...
                        nonlinear_elements{nl}.type '.']);
            end
        end
    end
else

    F = zeros(n*(H+1),1);
    dF = zeros(n*(H+1),size(dQ,2));
    for nl = 1:length(nonlinear_elements)
        % Sampling frequency
        tau = (0:2*pi/N:2*pi-2*pi/N)';
    
        if nonlinear_elements{nl}.islocal
            switch lower(nonlinear_elements{nl}.type)
                case 'custom'
                    
            end
    
        else % Global nonlinearity
            switch lower(nonlinear_elements{nl}.type)
                case 'polynomialstiffness'
                    % Inverse FFT
                    Qc = transpose(reshape(Q,[],H+1)); n0 = size(Qc,2);

                    if isfield(nonlinear_elements{nl}, 'Ind')
                        Indx = nonlinear_elements{nl}.Ind;
                        nx = length(Indx);
                    else
                        Indx = 1:n0;
                        nx = n0;
                    end

                    Qc = [Qc(1,Indx); Qc(2:end,Indx)/2; zeros(N-H-1,nx)];
                    Qc(end-H+1:end,:) = flipud(conj(Qc(2:H+1,:)));
                    q = ifft(Qc)*N;

                    % Exponents and coefficients
                    pp = nonlinear_elements{nl}.exponents;
                    Et = nonlinear_elements{nl}.coefficients;
                    
                    % Evaluate polynimials
                    z = prod(kron(q,ones(size(pp,1),1)) .^ repmat(pp,N,1),2);
                    
                    % Evaluate forces
                    nz = size(Et,2);
                    fnl = (Et*reshape(z,nz,N))';
                    
                    % Forward FFT
                    Fnl_ind = fft(fnl)/N;
                    Fnl_ind = [Fnl_ind(1,:);Fnl_ind(2:H+1,:)*2];
                    
                    Fnl = zeros(H+1,n0);
                    Fnl(:,Indx) = Fnl_ind;

                    % Add forces to global force vector
                    F = F + reshape(transpose(Fnl),[],1);
                    %% Analytical gradients 
                    % (If the double loop bothers you, vectorize!)
                     for l=1:nx
                        % Apply IDFT to Jacobian of Q
                        ndx = size(dQ,2);
                        dQlc = [dQ(Indx(l),:); dQ(n0+Indx(l):n0:end,:)/2; zeros(N-H-1,ndx)];
                        dQlc(end-H+1:end,:) = flipud(conj(dQlc(2:H+1,:)));
                        dql = ifft(dQlc)*N;
                        
                        % Derive
                        notl = setdiff(1:nx,l);
                        dzql_dql = repmat(pp(:,l),N,1).*...
                            (kron(q(:,Indx(l)),ones(size(pp,1),1)).^repmat(pp(:,l)-1,N,1));
                        dzql_dql(isnan(dzql_dql)) = 0;
                        dz_dql = prod([ dzql_dql ...
                            kron(q(:,Indx(notl)),ones(size(pp,1),1)).^...
                            repmat(pp(:,notl),N,1)], 2);
                        df_dql = (Et*reshape(dz_dql,nz,N))';
                        
                        % Derive
                        for i=1:nx
                            dfi = repmat(df_dql(:,i),1,size(dql,2)).*dql;
                            dFi = fft(dfi(end-N+1:end,:))/N;
                            dFi = [dFi(1,:);dFi(2:H+1,:)*2];
                            dF(Indx(i):n0:end,:) = dF(Indx(i):n0:end,:) + dFi;
                        end
                        
                    end
                otherwise
                    error(['Unknown global nonlinearity type ' ...
                        nonlinear_elements{nl}.type '.']);
            end
        end
    end
end
end