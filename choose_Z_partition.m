function [choice,zrands] = choose_Z_partition(C)
% input 'C' is a pxn matrix of community partitions (p partitions of n-node network)
% output 'choice' is the nx1 partition that has the highest Z-score with
% all the other partitions.
% ******currently IGNORES identical partitions....

[p,n] = size(C);

% compute zrand between all pairs
% keep track of partitions that are repeats
zrands = zeros(p);
repeats = zeros(p,1);
for i=1:p
    disp(i);
    for j=i+1:p
        zrands(i,j) = zrand_pair(C(i,:),C(j,:));
        [~,Cjnew,~] = comm_overlap(C(j,:)',C(i,:)');
        if all(C(i,:)==Cjnew') % if they map exactly
            repeats(j) = 1;
        end
    end
end

% remove repeat partitions from zrand
[zrands,newix] = remove_missings(zrands,find(repeats>0));

% make zrands symmetric
zrands = zrands + zrands';

% find partition with highest overall Z
ix = newix(find(sum(zrands)==max(sum(zrands)),1,'first'));
choice = C(ix,:)';




