function surety = sureties(origs,fpartn)
% origs = original p partitions (pxn)
% fprtn = final partitions (nx1)

[~,n] = size(origs);
fpartn = fpartn(:);
assert(numel(fpartn)==n);

% find pairwise percentages
assoc = mod_allegiance(origs);
assocf = mod_allegiance(fpartn');

% calculate sureties wrt final partition
surety = (sum(assoc.*assocf)./sum(assocf))';

end