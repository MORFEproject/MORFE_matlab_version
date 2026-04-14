function I_k = generateCombinations(k, I)
    n = numel(I);  % Ensure using numel to get the number of elements in I
    m = numel(k);  % Ensure using numel to get the number of elements in k

    % Precompute and store all possible column vectors
    combCache = cell(m, 1);
    for ni = 1:m
        combCache{ni} = findAllVectors(n, k(ni));
    end

    % Initialize output matrix
    I_k = [];
    % Call recursive function to generate all valid combinations
    I_k = recursiveCombine(1, zeros(n, m), combCache, I_k, I, m);

    function I_k = recursiveCombine(dim, currentMatrix, combCache, I_k, I, m)
        if dim > m
            % Check if the current matrix's row sums equal the vector I
            if all(sum(currentMatrix, 2) == I)
                I_k = cat(3, I_k, currentMatrix);  % Add the current matrix to the output
            end
            return;
        end

        % Get all possible vectors for the current dimension
        currentVectors = combCache{dim};
        for i = 1:size(currentVectors, 2)
            % Set the vector for the current dimension
            newMatrix = currentMatrix;
            newMatrix(:, dim) = currentVectors(:, i);
            % Recursively construct vectors for the remaining dimensions
            I_k = recursiveCombine(dim + 1, newMatrix, combCache, I_k, I, m);
        end
    end
end