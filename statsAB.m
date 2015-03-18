function [mu,sigma] = statsAB(n,Nab)
% for use with Rand z-score calculation from Traud, Kelsic, Mucha, & Porter 2010.
% gives mean and standard deviation of Wab under null model used in paper 
% (which preserves row&column sums, i.e. number and size of communities in
% each partition.
% INPUTS: n (number of nodes)
%         Nab (contingency matrix between partitions a and b)
% OUTPUTS: mu (expected number of pairs grouped together in both partitions)
%          sigma (standard deviation of # of pairs grouped together in 
%                 both partitions under the null model)

    % MU 
        % calculate total number of pairs M (consider missings???)
        M = n*(n-1)/2;
        % calculate Ma (number of pairs grouped together in a)
        Ma = sum(Nab,2); % gives a groupings   
        Ma = Ma(Ma>1);
        Ma = sum(factorial(Ma)./(2.*factorial(Ma-2)));
        % ...and Mb
        Mb = sum(Nab,1)'; % gives b groupings
        Mb = Mb(Mb>1);
        Mb = sum(factorial(Mb)./(2.*factorial(Mb-2)));
        mu = Ma*Mb/M;
        
    % ...and sigma (standard deviation of Wab (from null model))
        Ca = n*(n^2-3*n-2) - 8*(n+1)*Ma + 4*sum(sum(Nab,2).^3);
        Cb = n*(n^2-3*n-2) - 8*(n+1)*Mb + 4*sum(sum(Nab,1)'.^3);
        sigma = M/16 - (4*Ma-2*M)^2*(4*Mb-2*M)^2/(256*M^2) + Ca*Cb/(16*n*(n-1)*(n-2)) + ...
            ((4*Ma-2*M)^2-4*Ca-4*M)*((4*Mb-2*M)^2-4*Cb-4*M)/(64*n*(n-1)*(n-2)*(n-3));
        sigma = sqrt(sigma);

end