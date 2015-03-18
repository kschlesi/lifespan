function [rr,pp] = correlate(sample1,sample2)

sample1 = sample1(:);
sample2 = sample2(:);
if (length(sample1)~=length(sample2))
    error('samples must have the same number of elements');
end

[rr,pp] = corrcoef(sample1,sample2);
rr = rr(1,2);
pp = pp(1,2);

figure;
scatter(sample1,sample2);
legend(['Pearson''s r = ' num2str(rr) ', p = ' num2str(pp)]);

end