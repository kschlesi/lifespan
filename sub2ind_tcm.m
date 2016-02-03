function ind = sub2ind_tcm(s1_,s2_,layer)

ind = zeros(size(layer));

for i=1:numel(layer)

if layer(i)==1
    ind(i) = sub2ind([11,11],s1_(i),s2_(i));
end

if layer(i)==2
    ind(i) = (sub2ind([8,7],s1_(i),s2_(i)) + 121);
end

end

end