% mle example

% general expression to be minimized
    
   % mu = @(Phit,k,n,L) L.*Phit.^n./(Phit.^n+k); % form of function to be fit
    
% load data and create plot of rate function v. Phit

    theta = [0.005,0.005,0.01,0.005,0.008,0.01,0.03,0.12,0.24,0.19,0.27];
    H = 100.*[50,48,48,48,45,40,38,35,20,8,2];
    J = H.*theta;
    Phit = 0:0.1:1;
    %n=1; L=0.5;
    negLogLo = @(k) -1*sum((H-J).*log(1-(L.*Phit.^n./(Phit.^n+k))) ...
                                          + J.*log(L.*Phit.^n./(Phit.^n+k)));
    params0 = [0.5;1;0.5]; % starting parameters for mu function
    %f = @(k,n,L) negLogL(mu,Phit,k,n,L,J,H);
    params = fminunc(negLogLo, 0.5);

myfun = @(x,c) sum(x.^2+c);
c = 3;                              % define parameter first
x = fminunc(@(x) myfun(x,c),[1;1]);