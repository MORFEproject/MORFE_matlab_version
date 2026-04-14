%% A comparative study of DPIM and the modal approach on the stable response of a clamped-clamped von Kármán beam
clear;
close all;
clc

%% Parameters setting
parameter_model = struct('lhb',1,'hhb',5e-3,'bhb',0.12,'Ehb',69e9,'rhohb',2700,...
    'etahb',0.01); % Clamped-clamped von Kármán beam's parameters
parameter_setting = struct('nh_element',10,'style','GP','order',7,'location_excitation',1/2); % Setting 
Master_hb = 1; % Master modes for reduction

% excitation conditions
fextm = 1;
omega = 30*2*pi;

%% Build the dynamic equation for clamped-clamped von Kármán beam
% Obtain the full-order dynamic equations
Model_hostbeam = build_modelhb(parameter_setting.nh_element,...
    parameter_model,parameter_setting.location_excitation);

[Vm, D2] = eigs(Model_hostbeam.K, Model_hostbeam.M, 5, 'smallestabs');
Om_hb2 = diag(D2);
[Om_hb2,ind] = sort(Om_hb2);
Vm = Vm(:,ind);
Fre_hb = sqrt(Om_hb2)/2/pi;
disp(['The first-order natural frequency is ',num2str(Fre_hb(1)),' Hz'])

n_g = length(Model_hostbeam.Mu);
Ind = [];
Vector = [];
for ni = 1:n_g
    Ind = [Ind Model_hostbeam.Mu{ni}.I];
    Vector = [Vector Model_hostbeam.Mu{ni}.vector];
end
pp = full(Ind)';
Et = full(Vector);
Indx = 1:size(pp,2);
nonlinear_elements{1} = struct('type','polynomialstiffness','islocal',0,'ishysteretic',0,...
    'exponents',pp,'coefficients', Et, 'Ind', Indx);
Model_full = struct('M',Model_hostbeam.M,'D',Model_hostbeam.D,'K',Model_hostbeam.K,...
    'Ve',Model_hostbeam.fext,'ndof',Model_hostbeam.ndof);
Model_full.nonlinear_elements = nonlinear_elements;
Model_full.type = 'second order';
Model_full.fextm = fextm;

% Obtain the reduced-order model based on modal approach
% Reduction based on modal approach
RM_hb = ModalApproachRM(Model_hostbeam, Vm(:,Master_hb));

n_g = length(RM_hb.Mur);
Ind = [];
Vector = [];
for ni = 1:n_g
    Ind = [Ind RM_hb.Mur{ni}.I];
    Vector = [Vector RM_hb.Mur{ni}.vector];
end
pp = full(Ind)';
Et = full(Vector);
Indx = 1:size(pp,2);

nonlinear_elements{1} = struct('type','polynomialstiffness','islocal',0,'ishysteretic',0,...
    'Ind',Indx,'exponents',pp,'coefficients',Et);
Model_mam = struct('M',RM_hb.Mr,'D',RM_hb.Dr,'K',RM_hb.Kr,...
    'Ve',RM_hb.fextr,'ndof',length(Master_hb),'Vmc',RM_hb.Vm);
Model_mam.nonlinear_elements = nonlinear_elements;
Model_mam.type = 'second order';
Model_mam.fextm = fextm;

% Obtain the reduced-order model based on DPIM
% Reduction based on DPIM
RM_hb = getDPIMReducedModel(Model_hostbeam, Master_hb, parameter_setting);
Ah = eye(size(RM_hb.f{1}.I)); Bh = RM_hb.f{1}.vector; Vhe = RM_hb.Fextr;

n_g = length(RM_hb.f);
Ind = [];
Vector = [];
for ni = 2:n_g
    Ind = [Ind RM_hb.f{ni}.I];
    Vector = [Vector RM_hb.f{ni}.vector];
end
pp = full(Ind)';
Et = full(Vector);
Indx = 1:size(pp,2);

nonlinear_elements{1} = struct('type','polynomialstiffness','islocal',0,'Ind',Indx,...
    'exponents',pp,'coefficients',Et);
Model_dpim = struct('A',Ah,'B',Bh,'Ve',Vhe,'Psifun',RM_hb.Psifun, ...
    'Upsilonfun',RM_hb.Upsilonfun, 'ndof', size(Ah,1),'Indz',1:size(RM_hb.f{1}.I,2));
Model_dpim.nonlinear_elements = nonlinear_elements;
Model_dpim.type = 'first order';
Model_dpim.fextm = fextm;


%% Perform frequency response calculation with HBM embedded in Nlvib toolbox
% Starting and ending frequency
Om_s = 0.8*omega;
Om_e = 1.2*omega;

% Setting for HBM
analysis = 'frf'; % Analysis type
H = 3; % Truncation order of HBM
N = 2^7; % Sample point per period
dscale = 1;
ds = 2;
fscl = 1;
Indout = 5*3 - 1; % Output dof

% FRC calculation of full-order model
% Step 1: Find suitable initial condition
om = Om_s;
x0 = zeros((2*H+1)*Model_full.ndof,1);
Sopt = struct('flag',1,'fscl',fscl,'stepadapt',1,'stepmax',1,'reversaltolerance',0.1,...
    'Dscale',[1e-2*ones(size(x0));om],'ds',ds,'dsmin',ds/20,'dsmax',10*ds);
Solopt = optimset(optimset(@fsolve),'Display','off',...
        'Jacobian','on','MaxIter',50);%'DerivativeCheck','on');
X_HB_full0 = solve_and_continue(x0,...
@(X) HB_residual_mix(X,Model_full,H,N,analysis,[],Sopt.fscl),...
om,Om_e,ds,Sopt,{},Solopt);

% Step 2: Perform numerical continuation
x0 = X_HB_full0(1:end-1,end);
om = X_HB_full0(end,end);
Sopt = struct('flag',1,'fscl',fscl,'stepadapt',1,'stepmax',1e4,'reversaltolerance',0.1,...
'Dscale',[dscale*max(abs(x0))*ones(size(x0));om],...
'ds',ds,'dsmin',ds/20,'dsmax',10*ds);
X_HB_full = solve_and_continue(x0,...
@(X) HB_residual_mix(X,Model_full,H,N,analysis,[],Sopt.fscl),...
om,Om_e,ds,Sopt,{},Solopt);
Manifold_full.X_HB = X_HB_full; ndof = Model_full.ndof; Npoint = 2^9; type = 'frf';
Manifold_full = getOrbits(Manifold_full,ndof,Npoint,type);
for ni = 1:length(Manifold_full.Orbits)
    Xoutfull(ni) = max(abs(Manifold_full.Orbits{ni}(:,Indout)));
end

% FRC calculation of LNM-ROM
% Step 1: Find suitable initial condition
om = Om_s;
x0_hb = zeros((2*H+1)*Model_mam.ndof,1);
Sopt = struct('flag',1,'fscl',fscl,'stepadapt',1,'stepmax',1,'reversaltolerance',0.2,...
    'Dscale',[1e-2*ones(size(x0_hb));om],'ds',ds,'dsmin',ds/20,'dsmax',10*ds);
Solopt = optimset(optimset(@fsolve),'Display','off',...
        'Jacobian','on','MaxIter',50);%'DerivativeCheck','on');
X_HB_mam0 = solve_and_continue(x0_hb,...
@(X) HB_residual_mix(X,Model_mam,H,N,analysis,[],Sopt.fscl),...
om,Om_e,ds,Sopt,{},Solopt);

% Step 2: Perform numerical continuation
x0_hb = X_HB_mam0(1:end-1,end);
om = X_HB_mam0(end,end);
Sopt = struct('flag',1,'fscl',fscl,'stepadapt',1,'stepmax',1e4,'reversaltolerance',0.2,...
'Dscale',[dscale/2*max(abs(x0_hb))*ones(size(x0_hb));om]...
,'ds',ds,'dsmin',ds/20,'dsmax',10*ds);
X_HB_mam = solve_and_continue(x0_hb,...
@(X) HB_residual_mix(X,Model_mam,H,N,analysis,[],Sopt.fscl),...
om,Om_e,ds,Sopt,{},Solopt);
Manifold_mam.X_HB = X_HB_mam; ndof = Model_mam.ndof; Npoint = 2^9; type = 'frf';
Manifold_mam = getOrbits(Manifold_mam,ndof,Npoint,type);
for ni = 1:length(Manifold_mam.Orbits)
    Xoutmam(ni) = max(abs(Manifold_mam.Orbits{ni}*[Model_mam.Vmc(Indout,:)';zeros(size(Model_mam.Vmc,2),1)]));
end

% FRC calculation of NNM-ROM
% Step 1: Find suitable initial condition
ndof = Model_dpim.ndof;
om = Om_s;
x0_hb = zeros((2*H+1)*ndof,1);

Sopt = struct('flag',1,'fscl',fscl,'stepadapt',1,'stepmax',1,'reversaltolerance',0.1,...
    'Dscale',[1e-2*ones(size(x0_hb));om],'ds',ds,'dsmin',ds/20,'dsmax',10*ds);
Solopt = optimset(optimset(@fsolve),'Display','off',...
        'Jacobian','on','MaxIter',50);%,'DerivativeCheck','on');
t1 = tic;
X_HB_dpim0 = solve_and_continue(x0_hb,...
@(X) HB_residual_mix(X,Model_dpim,H,N,analysis,[],Sopt.fscl),...
om,Om_e,ds,Sopt,{},Solopt);
toc(t1)

% Step 2: Perform numerical continuation
x0_hb = X_HB_dpim0(1:end-1,end);
om = X_HB_dpim0(end,end);
Sopt = struct('flag',1,'fscl',fscl,'stepadapt',1,'stepmax',1e4,'reversaltolerance',0.1,...
'Dscale',[dscale/2*max(abs(x0_hb))*ones(size(x0_hb));om],'ds',ds,'dsmin',ds/20,'dsmax',10*ds);
t1 = tic;
X_HB_dpim = solve_and_continue(x0_hb,...
@(X) HB_residual_mix(X,Model_dpim,H,N,analysis,[],Sopt.fscl),...
om,Om_e,ds,Sopt,{},Solopt);
toc(t1)
Manifold_dpim.X_HB = X_HB_dpim; Npoint = N; type = 'frf';
Manifold_dpim = getOrbits(Manifold_dpim,ndof,Npoint,type);
for ni = 1:length(Manifold_dpim.Orbits)
    Xnf = Model_dpim.Psifun(Manifold_dpim.Orbits{ni}(:,1:2*length(Master_hb)));
    Vnf = Model_dpim.Upsilonfun(Manifold_dpim.Orbits{ni}(:,1:2*length(Master_hb)));
    Xoutdpim(ni) = max(abs(Xnf(:,Indout)));
end

%% Perform FRC comparison
figure
hold on
box on
plot(X_HB_full(end,:)/2/pi,Xoutfull,'-k')
plot(X_HB_mam(end,:)/2/pi,Xoutmam,'--r')
plot(X_HB_dpim(end,:)/2/pi,Xoutdpim,'--b')
le = legend('Full-order model','LNM-ROM', 'NNM-ROM','FontName','Times New Roman','FontSize',10);
set(le,'box','off')
xlabel('Excitation frequency (Hz)','FontName','Times New Roman','FontSize',10)
ylabel('Displacement (m)','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);
xlim([24 34])