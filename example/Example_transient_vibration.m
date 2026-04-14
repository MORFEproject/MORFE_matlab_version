%% A comparative study of DPIM and the modal approach on the transient response of a clamped-clamped von Kármán beam
clear;
close all;
clc

%% Parameters setting
parameter_model = struct('lhb',1,'hhb',5e-3,'bhb',0.12,'Ehb',69e9,'rhohb',2700,...
    'etahb',0.005); % Clamped-clamped von Kármán beam's parameters
parameter_setting = struct('nh_element',10,'style','GP','order',7,'location_excitation',1/2); % Setting 
Master_hb = 1; % Master modes for reduction

% excitation conditions
fextm = 20;
omega = 30*2*pi;

%% Build the dynamic equation for clamped-clamped von Kármán beam
Model_hostbeam = build_modelhb(parameter_setting.nh_element,...
    parameter_model,parameter_setting.location_excitation);

[Vm, D2] = eigs(Model_hostbeam.K, Model_hostbeam.M, 5, 'smallestabs');
Om_hb2 = diag(D2);
[Om_hb2,ind] = sort(Om_hb2);
Vm = Vm(:,ind);
Fre_hb = sqrt(Om_hb2)/2/pi;
disp(['The first-order natural frequency is ',num2str(Fre_hb(1)),' Hz'])

% Obtain the ode function for full-order dynamic equations
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

Model_fullOde = @(t, x, fext, omega) generateOdeFun(t, x, fext, omega, Model_full);

% Reduction based on modal approach
RM_hb = ModalApproachRM(Model_hostbeam, Vm(:,Master_hb));

% Obtain the ode function for reduced-order model based on modal approach
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
    'Ve',RM_hb.fextr,'ndof',length(Master_hb));
Model_mam.nonlinear_elements = nonlinear_elements;

Model_mamOde = @(t, x, fext, omega) generateOdeFun(t, x, fext, omega, Model_mam);


% Reduction based on DPIM
RM_hb = getDPIMReducedModel(Model_hostbeam, Master_hb, parameter_setting);
Ah = eye(size(RM_hb.f{1}.I)); Bh = RM_hb.f{1}.vector; Vhe = RM_hb.Fextr;

% Obtain the ode function for reduced-order model based on DPIM
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

Model_dpimOde = @(t, x, fext, omega) generateOdeFun(t, x, fext, omega, Model_dpim,'first order');


%% Perform time-domain comparison
tspan = [0,1];
x0 = zeros(Model_full.ndof*2,1);
[tnf,xnf] = ode45(@(t,x) Model_fullOde(t,x,fextm, omega), tspan, x0);

x0 = zeros(Model_mam.ndof*2,1);
[tnm,xnm] = ode45(@(t,x) Model_mamOde(t,x,fextm, omega), tspan, x0);

x0 = zeros(Model_dpim.ndof,1);
[tnd,znd] = ode45(@(t,x) Model_dpimOde(t,x,fextm, omega), tspan, x0);
xnd = Model_dpim.Psifun(znd);

figure
hold on
box on
plot(tnf,xnf(:,5*3-1),'-k')
plot(tnm,xnm(:,1)*Vm(5*3-1,Master_hb),'--r')
plot(tnd,xnd(:,5*3-1),'--b')
le = legend('Full-order model','LNM-ROM', 'NNM-ROM','FontName','Times New Roman','FontSize',10);
set(le,'box','off')
xlabel('$t$ (s)','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
ylabel('Displacement (m)','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);