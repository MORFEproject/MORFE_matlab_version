function Model = build_modelhb(nElements,parameter,xe)
%% Finite Element Setup
% Geometry
startLIN = tic;

if nargin <= 1
    l = 1;
    h = 5e-3;
    b = 0.12; 
    % Mesh parameters
    
    % Material properties
    
    E       = 69e9;  % 70e9 % 200e9 % Young's modulus
    rho     = 2700; % 2700 % 7850 % density
    nu      = 0.33;    % nu
    kappa   = 0; % material damping modulus 1e8
    etahb   = 0.005;
else
    l = parameter.lhb;
    h = parameter.hhb;
    b = parameter.bhb;

    E = parameter.Ehb;
    rho = parameter.rhohb;
    nu = 0.33;
    kappa = 0;
    etahb = parameter.etahb;
end

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor = @()BeamElement(b, h, myMaterial); % same element all across the domain

% Meshing the geometry
dx = l/nElements;
x = (0:dx:l).';
nNodes = size(x,1);
nodes = [x, zeros(nNodes,1)];
elements = [1:nNodes-1;2:nNodes].';

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% Plot mesh
%figure('Name','Mesh'); PlotMesh(nodes,elements,0);

%% Assemble linear stiffness, mass and damping
disp('Assembling M,C,K matrices')
% % parallelized assembly
% cluster = parcluster('local');
% cluster.NumWorkers = 4;
% parpool(cluster, 4)
% MyAssembly = Assembly(myMesh,true); 

MyAssembly = Assembly(MyMesh);
K = MyAssembly.stiffness_matrix();
M = MyAssembly.mass_matrix();
%C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition([1 nNodes],[1 2 3],0) % Cantilevered beam
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
%C = MyAssembly.constrain_matrix(C);



%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));

V = MyAssembly.unconstrain_vector(V0);
mod = 1;
v1 = reshape(V(:,mod),3,[]);
%figure;
%PlotFieldonDeformedMesh(nodes,elements,v1(1:2,:).','factor',0.2);
%title(['Mode ' num2str(mod) ', Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )

%% external force assembly
disp('Assembling external force vector')

outnode = MyMesh.nNodes;
if nargin > 2
    outdof = (floor(outnode*xe)+1)*3-1; % 激励点
else
    outdof = (floor(outnode/2)+1)*3-1;
end

outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdofvec = MyAssembly.constrain_vector(outdofvec);
outdof = find(outdofvec);

fext = outdofvec;


%weights = true(nElements,1); 
%fext = MyAssembly.constrain_vector(MyAssembly.uniform_body_force('weights',weights));


computationTimeLIN = toc(startLIN);
%% Tensor Assembly
disp('Getting nonlinearity coefficients')

fnl = cell(1,2);
disp('Assembling Tensors')
startTensors = tic;
fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);       
computationTimeTensors = toc(startTensors);
disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])

% apply boundary conditions
tic
for j = 1:length(fnl)
    fnl{j} = reduce_symmetric_tensor(MyAssembly.constrain_tensor(fnl{j}));
end
toc
Model = struct('M',M,'K',K,'fext',fext,'excitdof',outdof,'Nodes',nodes,'cind',[1 nNodes]);
Model.fnl = fnl;

N_hb = 50;
[Fre, Vm] = SpectrumAnalysis(Model.K, Model.M, N_hb);
disp(Fre(1:10));

D = etahb*pinv(Vm')*diag(sqrt(diag(Vm'*Model.K*Vm)))*pinv(Vm);

Model.D = D;

Mu = cell(length(Model.fnl),1);
for ni = 1:length(Model.fnl)
    Mu{ni} = change_tensor_Mu(Model.fnl{ni});
end
Model.Mu = Mu;
Model.ndof = size(Model.M,1);
end

% Spectral Analysis Function
function [fre, Vm] = SpectrumAnalysis(K, M, n)
    % Compute the eigenspectrum and modal mass
    [V, D] = eigs(full(K), full(M), n, 'smallestabs');
    O2 = diag(D);
    fre = sqrt(O2) / (2 * pi);
    
    % Compute mass-normalized mode shapes
    Vm = V ./ diag(V' * M * V)';
end

function T_reduced = reduce_symmetric_tensor(T)
    % Extract non-zero values and indices of the tensor
    [subs, vals] = find(T);
    
    % Sort the last 2:n+1 modes and normalize the mode patterns
    sorted_modes = [subs(:, 1) sort(subs(:, 2:end), 2)];
    
    % Find unique patterns and sum coefficients of identical patterns
    [unique_modes, ~, idx] = unique(sorted_modes, 'rows');
    new_vals = accumarray(idx, vals);
    
    % Reconstruct indices and tensor
    new_subs = unique_modes;
    T_reduced = sptensor(new_subs, new_vals, size(T));
end