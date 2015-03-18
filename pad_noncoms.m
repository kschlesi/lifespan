function padvals = pad_noncoms(vals,assignments)
% this is a specific helper function for choking_results; for multislice
% layers in which not all communities are represented, it takes "sigs" or
% a column vector representing the p-values of each extant community, of 
% size % (# extant communities)x1, and pads it with ones (p=1) in the
% correct slots for missing communities, as shown by the ordered list of 
% extant communities, "assignments," of size (# extant communities)x1.
% returns the padded value vector padvals, of size (max(assigmnents))x1.
    
    padvals = ones(max(assignments),1);
    padvals(assignments) = vals;

end