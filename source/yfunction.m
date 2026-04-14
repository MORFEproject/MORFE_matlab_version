function [yn,dyn] = yfunction(X_HB,ind,n,t,Phi_HB,type)
if nargin <= 5
    type = 'nma';
end

switch type
    case 'nma'
        if nargin <= 4
            phi = 0;
        else
            phi = Phi_HB(ind);
        end
        
        w = X_HB(end-2,ind);
        a = 10.^X_HB(end,ind);
        Qh = a*IVec(X_HB(1:end-3,ind),n);
        
        H = (size(Qh,2) - 1)/2;
        
        Hamonic = zeros(2*H+1,size(t,2));
        Hamonic(1,:) = ones(size(t));
        for i = 1:H
            Hamonic(2*i,:) = cos(i*w*(t - phi/w));
            Hamonic(2*i+1,:) = sin(i*w*(t - phi/w));
        end
        yn = (Qh*Hamonic)';
        dyn = (w*Qh*Delta(H)*Hamonic)';
    case 'frf'
        w = X_HB(end,ind);
        Qh = IVec(X_HB(1:end-1,ind),n);

        H = (size(Qh,2) - 1)/2;

        Hamonic = zeros(2*H+1,size(t,2));
        Hamonic(1,:) = ones(size(t));
        for i = 1:H
            Hamonic(2*i,:) = cos(i*w*t);
            Hamonic(2*i+1,:) = sin(i*w*t);
        end
        yn = (Qh*Hamonic)';
        dyn = (w*Qh*Delta(H)*Hamonic)';
end
end