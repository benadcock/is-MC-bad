%--- Description ---%
%
% Filename: find_margin.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: finds the reduced margin R(S) of a set S (see Sec SM2.3)
%
% Inputs:
% S - an d x n set of multi-indices
%
% Output:
% RS - the reduced margin of S

function RS = find_margin(S)

%%% find margin %%%

[d,n] = size(S);

MS = [];
for i = 1:n
    
    for j = 1:d
        
        % add e_j to the ith entry of S
        z = S(:,i);
        z(j) = z(j) + 1;
        
        % test for membership of S
        if ismember(z',S','rows') == 0
            MS = [MS z];
        end
        
    end
    
end

MS = (unique(MS','rows'))'; % eliminate repeats

%%% compute reduced margin %%%

RS = [];

for i = 1:size(MS,2)
    
    z = MS(:,i);
    
    flg = 1;
    
    for j = 1:d
        
        if z(j) > 0
            w = z; w(j) = z(j) - 1;
            if ismember(w',S','rows') == 0
                flg = 0;
            end
        end
        
    end
    
    if flg == 1
        RS = [RS z];
    end
    
end

end

