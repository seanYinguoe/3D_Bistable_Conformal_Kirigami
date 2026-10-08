function intersection = findBisectorIntersection(A1, A2, B1, B2)
    % Calculate midpoints of the vectors (works for column vectors)
    midA = (A1 + A2) / 2;
    midB = (B1 + B2) / 2;
    
    % Directions of original vectors (column format preserved)
    dirA = A2 - A1;
    dirB = B2 - B1;
    
    % Perpendicular directions (column vectors, 90-degree rotation)
    perpA = [-dirA(2); dirA(1)];  % Semicolon creates column vector
    perpB = [-dirB(2); dirB(1)];
    
    % Solve linear system A*[t; s] = B
    A = [perpA, -perpB];  % Combine column vectors
    B = midB - midA;
    
    if det(A) == 0
        error('Perpendicular bisectors are parallel.');
    else
        ts = A \ B;  % Solve system
        intersection = midA + ts(1)*perpA;
    end
end
