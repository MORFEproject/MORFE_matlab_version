function dx = generateOdeFun(t, x, fext, omega, Model, type)
FlagCustom = 0;
if isa(fext, 'function_handle')
    FlagCustom = 1;
end
if (nargin > 5) & strcmp(type,'first order')
    A = Model.A; B = Model.B;
    Ve = Model.Ve; ndof = Model.ndof;
    nonlinear_elements = Model.nonlinear_elements;
    if FlagCustom
        dx = real(A\(B*x + Fn(x, ndof, nonlinear_elements, type) + fext(t)*Ve));
    else
        dx = real(A\(B*x + Fn(x, ndof, nonlinear_elements, type) + fext*Ve*cos(omega*t)));
    end
else
    M = Model.M; D = Model.D; K = Model.K;
    Ve = Model.Ve; ndof = Model.ndof;
    nonlinear_elements = Model.nonlinear_elements;
    
    A = blkdiag(M,M);
    B = [zeros(size(M)) M;
        -K -D];
    if FlagCustom
        dx = real(A\(B*x + Fn(x, ndof, nonlinear_elements) + fext(t)*[zeros(size(Ve));Ve]));
    else
        dx = real(A\(B*x + Fn(x, ndof, nonlinear_elements) + fext*[zeros(size(Ve));Ve]*cos(omega*t)));
    end
end
end

function fn = Fn(x, ndof, nonlinear_elements, type)
if (nargin > 3) & strcmp(type, 'first order')
    fn = zeros(ndof, 1);
    n_ne = length(nonlinear_elements);
    for ni = 1:n_ne
        if nonlinear_elements{ni}.islocal
            switch nonlinear_elements{ni}.type
                case 'DPIM'
                    fn = fn + nonlinear_elements{ni}.Fn(x);
            end
        else
            switch nonlinear_elements{ni}.type
                case 'polynomialstiffness'
                    Ind = nonlinear_elements{ni}.Ind;
                    pp = nonlinear_elements{ni}.exponents;
                    Et = nonlinear_elements{ni}.coefficients;
                    z = prod(kron(x(Ind)',ones(size(pp,1),1)) .^ pp,2);
                    nz = size(Et,2);
                    fn(Ind) = fn(Ind) + (Et*reshape(z,nz,1));
            end
        end
    end
else
    u = x(1:ndof); v = x(ndof+1:end);
    fn = zeros(2*ndof,1);
    n_ne = length(nonlinear_elements);
    for ni = 1:n_ne
        if nonlinear_elements{ni}.islocal
            switch nonlinear_elements{ni}.type
                case 'custom'
                    dn = nonlinear_elements{ni}.force_direction;
                    fn(ndof+1:end) = fn(ndof+1:end) - dn*nonlinear_elements{ni}.NF(dn'*u,dn'*v);
            end
        else
            switch nonlinear_elements{ni}.type
                case 'polynomialstiffness'
                    Ind = nonlinear_elements{ni}.Ind;
                    pp = nonlinear_elements{ni}.exponents;
                    Et = nonlinear_elements{ni}.coefficients;
                    z = prod(kron(u(Ind)',ones(size(pp,1),1)) .^ pp,2);
                    nz = size(Et,2);
                    fn(ndof+Ind) = fn(ndof+Ind) - (Et*reshape(z,nz,1));
            end
        end
    end
end
end