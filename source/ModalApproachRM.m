% Model Order Reduction Function via Modal Approach
function RM = ModalApproachRM(FullModel, Vm)
% Calculate reduced-order mass, damping, stiffness matrices, and external force vector
Mr = diag(diag(Vm' * FullModel.M * Vm));
Kr = diag(diag(Vm' * FullModel.K * Vm));
Dr = diag(diag(Vm' * FullModel.D * Vm));

fextr = Vm' * FullModel.fext;
% Construct the reduced-order model structure
RM = struct('Mr', Mr, 'Dr', Dr, 'Kr', Kr, 'fextr', fextr, 'Vm', Vm, 'FullModel', FullModel);
RM.ndof = size(Mr, 1);
% If the original model contains nonlinear forces, perform order reduction
if isfield(FullModel, 'fnl')
    n_p = length(FullModel.fnl); % Order of nonlinear forces
    fnlr = cell(1, n_p);
    Mur = fnlr;
    for ni = 1:n_p
        VmcT = cell(1, ni + 2);
        [VmcT{:}] = deal(Vm');
        % Perform modal transformation on the nonlinear force tensor
        fnlr{ni} = ttm(FullModel.fnl{ni}, VmcT, 1:ni+2);
        % Change the tensor's Mu (calls an external function with unknown specific functionality)
        Mur{ni} = change_tensor_Mu(fnlr{ni});
    end
    
    RM.fnlr = fnlr;
    RM.Mur = Mur;
end
end