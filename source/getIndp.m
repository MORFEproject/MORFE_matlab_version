function Indp = getIndp(n, n_p, style)
    I = cell(n_p,1);
    CS = style.CS;
    switch CS
        case 'index'
            for p_i = 1:n_p
                I{p_i} = recurse([], n, p_i)';
            end
        case 'multi_index'
            for p_i = 1:n_p
                I{p_i} = flip(findAllVectors(n, p_i), 2);
            end
    end
    Indp.I = I;
    Indp.type = CS;
end

function vectors = recurse(prefix, n, p_i)
    if p_i == 0
        vectors = prefix;
    else
        vectors = [];
        for i = 1:n
            newPrefix = [prefix i]; 
            newVectors = recurse(newPrefix, n, p_i-1);
            vectors = [vectors; newVectors];
        end
    end
end
