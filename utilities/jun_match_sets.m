function [match, nomatch_A, nomatch_B] = jun_match_sets(A, B)
% Matches the elements of A and B.
%
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

nA = length(A);
nB = length(B);

if length(unique(A)) ~= nA
    disp('A must be unique');
    return
elseif length(unique(B)) ~= nB
    disp('B must be unique');
    return
end

validA = true(nA,1);
validB = true(nB,1);

match = zeros(min(nA,nB), 2);
counter = 0;
for a = 1:nA
    ind = find(A(a) == B(validB));
    
    if ind
        counter = counter + 1;
        
        bs = find(validB);
        b = bs(ind);
        
        match(counter,:) = [a b];
        validA(a) = false;
        validB(b) = false;
    end
end

match = match(1:counter,:);
nomatch_A = find(validA);
nomatch_B = find(validB);

end
