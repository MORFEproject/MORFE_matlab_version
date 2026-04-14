function R = findAllVectors(n, p)
    R = [];
    R = recursiveFind(0, n, p, [], R);
end

function R = recursiveFind(currentDim, totalDim, remainingSum, currentVector, R)
    if currentDim == totalDim - 1
        R = [R, [currentVector; remainingSum]];
    else
        for i = 0:remainingSum
            R = recursiveFind(currentDim + 1, totalDim, remainingSum - i, [currentVector; i], R);
        end
    end
end
