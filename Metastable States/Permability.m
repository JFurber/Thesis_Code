function [newMatrix, countOffDiagonalNonZero] = Permability(transitionMatrix)
    % Get the size of the matrix
    [k, ~] = size(transitionMatrix);
    
    % Initialize the new matrix with NaNs
    newMatrix = NaN(k, k);
    
    % Extract the off-diagonal elements that are not zero
    offDiagonalElements = transitionMatrix(~eye(k) & transitionMatrix ~= 0);
    
    % Calculate the average of the off-diagonal elements
    avgOffDiagonal = mean(offDiagonalElements);
    
    % Count the number of off-diagonal non-zero elements
    countOffDiagonalNonZero = numel(offDiagonalElements);
    
    % Iterate over the matrix to fill the new matrix
    for i = 1:k
        for j = 1:k
            if i ~= j && transitionMatrix(i, j) ~= 0
                newMatrix(i, j) = transitionMatrix(i, j) - avgOffDiagonal;
            end
        end
    end
end

