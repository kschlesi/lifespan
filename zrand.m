function [zrand,zdist] = zrand(partitions)
% calculates the Rand z-score for p partitions of an n-node network
% INPUTS:  partitions: pxn matrix of community assignments
% OUTPUTS: zrand: the average Rand z-score over all pairs of partitions
%          zdist: a symmetrix pxp matrix containing the z-score for each
%                 pair of partitions

[p,n] = size(partitions);
zdist = zeros(p);
for a=1:p   % loops over all unique pairs of partitions
    aa=a
    for b=a:p
%         Nab = zeros(max(partitions(a,:)),max(partitions(b,:)));
%         for i=1:size(Nab,1) % loops over communities in partition a
%             for j=1:size(Nab,2) % ...and those in partition b
%                 Nab(i,j) = sum((partitions(a,:)==i).*(partitions(b,:)==j));
%             end
%         end
%         % mu and sigma, computed from null model
%         [mu,sigma] = statsAB(n,Nab);
%         % ...and Wab (number of pairs grouped together in both a and b)
%         Wab = Nab(Nab>1);
%         Wab = sum(factorial(Wab)./(2.*factorial(Wab-2)));
%         % calculate Zab and enter into zdist
%         zdist(a,b) = (Wab - mu)/sigma;  
        
        zdist(a,b) = zrand_pair(partitions(a,:),partitions(b,:));
    end
end

zdist = triu(zdist)+triu(zdist,1)'; % fill in lower triangle for symmetry
% zrand = average over only the upper triangle (NOT the diagonal)
zrand = 2*sum(sum(triu(zdist,1)))/(p*(p-1));

end