function [submat,newix] = remove_missings(origmat,missings)

[u,n] = size(origmat);
if u~=n; error('input matrix must be square!'); end;

newix = removeval(1:n,missings)';

submat = origmat(newix,:);
submat = submat(:,newix);

end