function flex = flexibility(partns,iscat)
% input: a txn matrix (t partitions, one per slice, for n nodes)
% input: iscat, a boolean; 1 if categorical
% output: a nx1 vector (flexibility of each node)

if nargin<2
    iscat = 0; % default = time-ordered
end

n = size(partns,2);
t = size(partns,1);
flex = zeros(n,1);

if iscat % categorical

   for T=1:t-1
       for T2=T+1:t
           switches = (partns(T,:)~=partns(T2,:)); 
           flex = flex + switches';
       end
   end
   flex = flex*2/(t*(t-1));

else % time-ordered
   
    for T=1:t-1
        switches = (partns(T,:)~=partns(T+1,:));
        flex = flex + switches';
    end
    flex = flex/(t-1);
end