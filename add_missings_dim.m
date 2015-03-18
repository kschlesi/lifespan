function origdim = add_missings_dim(inmat,newix,n,dim,pad)

if size(inmat,dim)~=numel(newix)
    error('the desired dimension is the wrong size');
end

if ~isnumeric(pad) || ~isscalar(pad)
    error('pad value must be numerical scalar');
end

newsize = size(inmat); % size of inmat
newsize(dim) = n; % set newsize to size of output
tdims = ndims(inmat); % total number of dimensions of input
if length(inmat)==numel(inmat); tdims = 1; end; % special case if 1-dim
origdim = shiftdim(zeros(newsize),dim-1); % put shortened dimension in front
inmat = shiftdim(inmat,dim-1); % same with inmat
evalstring = repmat(',:',1,tdims-1);
for ix=1:n
    if ismember(ix,newix)
        eval(['origdim(ix' evalstring ') = inmat(newix==ix' evalstring ');']);
    else
        eval(['origdim(ix' evalstring ') = pad.*ones(size(origdim(ix' evalstring ')));']);
    end
end

if tdims~=1
    origdim = shiftdim(origdim,tdims-dim+1); % shift back to original dimensions
end

end