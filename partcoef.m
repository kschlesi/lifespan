function [slicecoefs] = partcoef(C,A)
% identifies the participation coefficient of each node in a multislice
% network. 
% input: C, an nxt matrix of multislice community assignments
%        A, a tx1 cell array of nxn weighted, undirected adjacency matrices
% output: slicecoefs, an nxt matrix containing the coefficient for each of
%                     the n nodes in the current slice only
% ADAPTED from Brain Connectivity Toolbox (please cite),
%              https://sites.google.com/site/bctnet/measures/list

[n,t] = size(C);
k = cell(t,1);
slicecoefs = zeros(size(C));
for T=1:t
    k{T} = sum(A{T},2);
    Gc = (A{T}~=0)*diag(C(:,T));
    kin2 = zeros(n,1);
    for i=1:max(C(:,T))
        kin2 = kin2 + sum(A{T}.*(Gc==i),2).^2;
    end
    kdiv2 = kin2./(k{T}.^2);
    slicecoefs(:,T) = ones(n,1) - kdiv2;
    slicecoefs(~k{T},T) = 0;
    %slicecoefs(:,T) = participation_coef(A{T},C(:,T));
end

end
