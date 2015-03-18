function [C,Q] = genlouvainREPs(A,p,gamma,omega)
% this code does p community detection optimizations on adjacency matrix b
% using generalized Louvain algorithm from netwiki package.
% genlouvain returns S, an nx1 vector of assignments, and a 1x1 Q-value.

% inputs: A (original NxN (potentially sparse) adj matrix for static)
%           OR (Tx1 cell array of T NxN matrices for categorical multislice)
%         gamma (value of gamma for use in community optimization penalty)
%         omega (value of omega for multislice; omit or set = 0 for static)
%         p (number of genlouvain optimizations to perform)
%
% outputs: C (pxN matrix of community assignments)
%            OR (pxNxT tensor of multislice community assignments)
%          Q (px1 matrix of associated quality values)

if nargin<4 || ~omega  % static community detection, A is nxn matrix
    N = size(A,1); % number of nodes in the system
    k = full(sum(A)); % vector of node degrees
    twom = sum(k); % m = number of edges in system
    B = full(A - gamma*(k'*k)/twom);% or define B as fcn handle to improve memory:
    %B = @(v) A(:,v) - gamma*k'*k(v)/twom; % (see 'genlouvain.m' example)
    C = zeros(p,N); % pxN matrix of assignments  
    Q = zeros(p,1);
    for i=1:p  % performing p optimizations
        disp(i);
        [S,Q1] = genlouvain(B,10000,0);
        C(i,:) = S';
        Q(i) = Q1/twom;
    end

else % multislice community detection, A is Tx1 array of NxN matrices
    N=length(A{1});
    T=length(A);
    B=spalloc(N*T,N*T,(N+T)*N*T);
    twomu=0;
    for s=1:T
        k=sum(A{s});
        twom=sum(k);
        twomu=twomu+twom;
        indx=(1:N)+(s-1)*N;
        B(indx,indx)=A{s}-gamma*(k'*k)/twom;
    end
    twomu=twomu+T*omega*N*(T-1);
    all2all = N*[(-T+1):-1,1:(T-1)];
    B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    C = zeros(p,N,T);
    Q = zeros(p,1);
    for i=1:p  % performing p optimizations
        disp(i);
        [S,Q1] = genlouvainmex(B);
        Q(i) = Q1/twomu;
        C(i,:,:) = reshape(S,N,T);
    end
end

end