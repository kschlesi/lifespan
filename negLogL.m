function y = negLogL(k,n,L,mu,Phit,J,H)

    y = -1*sum((H-J).*log(1-mu(Phit,k,n,L)) + J.*log(mu(Phit,k,n,L)));

end