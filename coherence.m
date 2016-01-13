function [maxcohr, maxcomm, cohrs] = coherence(idREG,partn)
% COHERENCE returns cohrs, a measure of the coherence of a subset of m nodes
% in a partition of n nodes.
% inputs: idREG, an mx1 vector containing the IDs of the m nodes in subset
%         partn, an nxt matrix of c cluster assignments for t layers of n nodes
% output: cohrs, a cxt vector containing the fraction of the m subset nodes
%                contained in each cluster in each layer
%         maxcohr, the maximum fraction contained in cohrs
%         maxcomm, the ID of the cluster containing the most idREG nodes

idREG = idREG(idREG>0);
m = length(idREG);
t = size(partn,2);
c = numel(unique(partn));
cohrs = zeros(c,1);
for T=1:t
    for i=1:c
        cohrs(i,T) = sum(partn(idREG,T)==i)./m;
    end
end
[maxcohr,maxcomm] = max(cohrs',[],2);