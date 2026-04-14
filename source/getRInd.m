function R = getRInd(Psi, f, Indp, I, CS)

if nargin <= 4
    CS = 'index';
end

N = size(Psi{1},1);
n = size(f{1},1);

R = zeros(N, 1);
Ind = Indp.I;

switch CS
    case 'index'
        p = length(I);
        for s = 1:n
            for k = 2:p-1
                for l = 0:p-k
                    p_1 = p - k + 1;
                    p_2 = k;
                    I1 = [I(1:l);s;I(l+k+1:p)];
                    I2 = I(l+1:l+k);
        
                    j1 = find_index_VinA(I1,Ind{p_1});
                    j2 = find_index_VinA(I2,Ind{p_2});
        
                    Rskl = Psi{p_1}(:, j1)*f{p_2}(s, j2);
                    R = R + Rskl;
                end
            end
        end

    case 'multi_index'
        p = sum(I);
        for s = 1:n
            for k = 2:p-1
                p_1 = p - k + 1;
                p_2 = k;

                mp2 = size(Ind{p_2},2);
                for kni = 1:mp2
                    I2 = Ind{p_2}(:,kni);
                    I1 = I-I2; I1(s) = I1(s)+1;
                    if min(I1) >= 0
                        j1 = find_index_VinA(I1, Ind{p_1});
                        Rskn = I1(s)*Psi{p_1}(:,j1)*f{p_2}(s,kni);
                        R = R + Rskn;
                    else
                        continue
                    end
                end
            end
        end
end
end
