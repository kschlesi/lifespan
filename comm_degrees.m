function [cSize,inND,outND] = comm_degrees(A,partnS,nodeix,n,t,nruns,nsubs)

% calculate node degrees
%overallND = zeros(n,nruns,nsubs);
inND = zeros(n,nruns,nsubs);
outND = zeros(n,nruns,nsubs);
cSize = zeros(n,nruns,nsubs);
for k=1:nsubs
  for T=1:t
    weights = A(k).adj{T};
    sameC = mod_allegiance([partnS(~~nodeix,T,k),partnS(~~nodeix,T,k)]',0);
    %overallND(~~nodeix,T,k) = sum(weights)';
    inND(~~nodeix,T,k) = sum(weights.*sameC)';
    outND(~~nodeix,T,k) = sum(weights.*~sameC)';
    cSize(~~nodeix,T,k) = sum(sameC)';
  end   
end

end