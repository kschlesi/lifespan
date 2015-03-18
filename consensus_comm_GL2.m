function [S2,Q2,X_new3,qpc] = consensus_comm_GL2(C)

% FUNCTION: To identify consistent communities in an ensemble
% 
% OUTPUT:
% S2 is a pxn matrix of new community assignments
% Q2 is the associated modularity value
% X_new3 is the thresholded nodal association matrix
% qpc is the quality of the consensus (the lower the value the better)
%
% INPUT:
% C is a pxn matrix of community assignments where p is the number of optimizations
% that have been performed. n is the number of nodes in the system. C gives the real partitions. 

% find size of network & number of optimizations
npart = numel(C(:,1));
m = numel(C(1,:));

% initialize
X = zeros(m,m);
X_rand3 = X;
C_rand3 = zeros(size(C));

% try a random permutation approach
for i = 1:npart;
    pr = randperm(m);
    C_rand3(i,:) = C(i,pr);  % C_rand3, npart random assignments of m nodes
end
for i = 1:npart; % creates X (an mxm count of how many times nodes i and j
                 %            were assigned to the same community fo realz)
    ii=i         % & X_rand3 (an mxm count of how many times nodes i and j
                 %            were assigned to the same community by chance)
    tic
    for k = 1:m
        for p = 1:m;
            if isequal(C(i,k),C(i,p))
                X(k,p) = X(k,p)+1;
            else
                X(k,p) = X(k,p)+ 0;
            end         
            if isequal(C_rand3(i,k),C_rand3(i,p))
                X_rand3(k,p) = X_rand3(k,p)+1;
            else
                X_rand3(k,p) = X_rand3(k,p)+ 0;
            end
        end
    end
    toc
end

% keep only associated assignments that occur more often than expected in
% the random data (assign them to X_new3, an mxm count with nonsigs zeroed)
X_new3 = zeros(m,m);
X_new3(X>max(max(triu(X_rand3,1)))) = X(X>max(max(triu(X_rand3,1))));

neqs = (X_new3~=X_new3');
if sum(sum(neqs))
    error('nodal association matrix not symmetric!!');
end


% recompute optimal partition on this new matrix of kept community
% association assignments (output: S2, the new pxn assignment matrix)
[S2,Q2] = genlouvainREPs(X_new3,npart,1);
%X_new3 = X_new3./p;

% define the quality of the consensus
qpc = sum(sum(abs(diff(S2))));


% visualize
% [MATreordered] = reorderMAT(X_new3,50000,'line');
% figure; imagesc(MATreordered);
% does this 


% 