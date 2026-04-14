function RM_DPIM = getDPIMReducedModel(Model_hostbeam, master, parameter_setting)
Model_hostbeam_DPIM = struct('M',Model_hostbeam.M, 'C',Model_hostbeam.D,'K',Model_hostbeam.K,...
    'N', Model_hostbeam.ndof,'nonlinear_force_type','tensor');
nonlinear_force = cell(3,1); nonlinear_force{2} = Model_hostbeam.fnl{1};
nonlinear_force{3} = Model_hostbeam.fnl{2};

Model_hostbeam_DPIM.nonlinear_force = nonlinear_force;
% Manifold setting
Nmax = 10;
errorBar = 0.1;
p = 3;
imod = master;
ParameterisationStyle = 'GP'; % or CNF, RNF, GP
CalculationStyle = 'multi_index'; % index or mutli_index
ModeAnalysisStyle = 'RealMode'; % RealMode or ComplexMode
ManifoldStyle = 'zero'; % zero or lsqminnorm
CalculationMethod = 'Step'; % Total or Step
if nargin > 2
    p = parameter_setting.order;
    ParameterisationStyle = parameter_setting.style;
end
style = struct('PS',ParameterisationStyle,'CS',CalculationStyle,'MAS',ModeAnalysisStyle,...
    'MS',ManifoldStyle, 'CM', CalculationMethod);
Model_hostbeam_DPIM = spectrum_analysis(Model_hostbeam_DPIM,Nmax,style);

Master = [imod imod+Nmax];
Slave = setdiff(1:2*Nmax,Master);

Indp = getIndp(length(Master), p, style);

Model_hostbeam_DPIM.Indp = Indp;
Model_hostbeam_DPIM = checkResonance(Model_hostbeam_DPIM, Master, Slave, Indp, errorBar, style);

tic
Compact_ReducedModel = InvariantManifold(Model_hostbeam_DPIM, style);
toc

% Generate coordinate transform matrix
nm = length(imod);
invT = zeros(2*nm);
for ni = 1:length(imod)
    invT([ni ni+nm], [ni ni+nm]) = [1 1;
        1/1i -1/1i];
end
T = inv(invT);

% Transfer the complex model to real one
tic
RealModel = getRealModel(Compact_ReducedModel,T);
toc

ManifoldFun = getManifoldFun(RealModel, 'real');

RM_DPIM = RealModel; RM_DPIM.Psifun = ManifoldFun.Psi;
RM_DPIM.Upsilonfun = ManifoldFun.Upsilon; RM_DPIM.T = T;
RM_DPIM.Xr = Model_hostbeam_DPIM.spectrum.X(:,Master);
RM_DPIM.Yr = Model_hostbeam_DPIM.spectrum.Y(:,Master);
Fext = [Model_hostbeam.fext;zeros(size(Model_hostbeam.fext))];
RM_DPIM.Fextr = inv(T)*RM_DPIM.Xr'*Fext; RM_DPIM.TX = inv(T)*RM_DPIM.Xr';
RM_DPIM.Ndof = 2*nm;
end