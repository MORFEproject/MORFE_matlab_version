function Manifold = getOrbits(Manifold, ndof, Npoint, type)

if nargin < 4
    type = 'nma';
end

X_HB = Manifold.X_HB;
norbits = size(X_HB, 2);
Orbits = cell(1, norbits);
Tn = Orbits;

for ni = 1:norbits
    switch type
        case 'nma'
            tn = linspace(0, 2*pi/X_HB(end-2, ni), Npoint);
        case 'frf'
            tn = linspace(0, 2*pi/X_HB(end, ni), Npoint);
    end
    [Yn, dYn] = yfunction(X_HB, ni, ndof, tn, zeros(ni,1), type);
    Orbits{ni} = [Yn dYn];
    Tn{ni} = tn;
end
Manifold.Orbits = Orbits;
Manifold.Tn = Tn;
end