function ManifoldFun = getManifoldFun(ReducedModel,varargin)
if (nargin > 1)&&strcmp(varargin(1),'real')
    Psi = ReducedModel.Psi;
    Upsilon = ReducedModel.Upsilon;
    PsiFun = @(z) RealManifoldFunction(z, Psi);
    Upsilon = @(z) RealManifoldFunction(z, Upsilon);
else
    Psi = ReducedModel.Psi;
    Upsilon = ReducedModel.Upsilon;
    I = ReducedModel.I;
    PsiFun = @(z) ManifoldFunction(z, Psi, I);
    Upsilon = @(z) ManifoldFunction(z, Upsilon, I);
end

ManifoldFun.Psi = PsiFun;
ManifoldFun.Upsilon = Upsilon;
end

function Fun = ManifoldFunction(z, Psi, I)
PsiTotal = [];
ITotal = [];
for ni = 1:length(Psi)
    PsiTotal = [PsiTotal Psi{ni}];
    ITotal = [ITotal I{ni}];
end
ITotal = ITotal';

if size(z,1)>1
    zn = prod(kron(z,ones(size(ITotal,1),1)) .^ repmat(ITotal,size(z,1),1),2);
    Fun = (PsiTotal*reshape(zn,size(PsiTotal,2),[])).';
else
    zn = prod(z.^ITotal,2);
    Fun = PsiTotal*zn;
end
end

function Fun = RealManifoldFunction(z, Psi)
PsiTotal = [];
ITotal = [];
for ni = 1:length(Psi)
    PsiTotal = [PsiTotal Psi{ni}.vector];
    ITotal = [ITotal Psi{ni}.I];
end
ITotal = ITotal';
if size(z,1)>1
    zn = prod(kron(z,ones(size(ITotal,1),1)) .^ repmat(ITotal,size(z,1),1),2);
    Fun = (PsiTotal*reshape(zn,size(PsiTotal,2),[])).';
else
    zn = prod(z.^ITotal,2);
    Fun = PsiTotal*zn;
end
end