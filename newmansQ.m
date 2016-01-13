function Q = newmansQ(atest,Ctest)
% computes Newman's Q-value for modularity of a network with nxn adjacency
% matrix A, given a partition C.

P = full(sparse(Ctest,1:length(Ctest),1));
Q = sum(sum((P*atest).*P))/sum(sum(atest));

% construct matrix e
%c = unique(Ctest);
% disp(c);
% e = zeros(numel(c));
% for i=1:numel(c)
%   e(i,i) = sum(sum(atest(Ctest==i,Ctest==i)));
%   for j=i+1:numel(c)
%     % count total weight between comm i and comm j
%     e(i,j) = sum(sum(atest(Ctest==i,Ctest==j)));
%   end
% end
% figure; bcolor(e);
% %assert(sum(sum(e)) == sum(sum(triu(atest))));
% e = e + triu(e,1)';
% Q = trace(e) - sum(sum(e^2));
% Q = Q/sum(sum(atest));

end