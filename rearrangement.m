function remat = rearrangement(orig, new, newix)
% this function computes the rearrangement matrix given orig partitions,
% new partitions (both partition matrices pxn -- runs x nodes), and
% "newix", which maps new to old indices: newix(new_index) = old_index

orig_small = orig(:,newix);
assert(all(size(orig_small) == size(new)));
remat = mod_allegiance(orig_small) - mod_allegiance(new);

end