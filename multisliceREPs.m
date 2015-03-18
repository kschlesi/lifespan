function [S,Q] = multisliceREPs(A,p,gamma,omega)

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
S = zeros(p,N,T);
Q = zeros(p,1);
for i=1:p
[S1,Q1] = genlouvainmex(B);
Q(i) = Q1/twomu;
S(i,:,:) = reshape(S1,N,T);
end
