function [fits,fix,nonfits,nix] = intervalSort(X,div,scheme,pct)
% given input intervals and input div, returns a list of intevals
% that fit in the divisions and a list of those that do not, based on scheme
% INPUTS: X, a Tx2 list of intervals 
%            (X(:,1)=beginning times) and X(:,2)=ending times)
%         div, a Dx2 list of divisionsinto which to sort
%            (div(:,1))=beginning times and div(:,2)=ending times)
%         scheme = a string naming the desired sorting scheme
%            ('within' = intervals entirely enclosed in divisions (default), 
%             'pct' = intervals at least pct% within a division (default pct=0.5),
%             'begin' = intervals begin within a division)
% OUTPUTS: fits, an Fx2 list of the entries of 'X' that fit sorting criteria
%          fix, an Fx2 list of the row of 'div' that each interval fit into
%          nonfits, a (T-F)x2 list of the entries of 'X' that were not sorted
%          nix, a list of the rows of 'div' that had no intervals in them

% test whether intervals are all positive
if sum((X(:,2)-X(:,1))<0)
    error('All input intervals must have positive length!');
end

if nargin<3
    warning('Using default sorting scheme ''within''')
    scheme = 'within';
end

if (strcmp(scheme,'within')+strcmp(scheme,'pct')+strcmp(scheme,'begin'))==0
    warning('Using default sorting scheme ''within''')
    scheme = 'within';
end

if nargin<4
    pct = 0.5;
end

if pct<0 || pct>1
    error('Percent input must be between 0 and 1!');
end
    
% sorting!!!
fits = zeros(size(X));
nonfits = zeros(size(X));
fix = zeros(size(X));
nix = ones(size(div,1),1);
f = 0; % counts number of fits
nf = 0; % counts number of nonfits
switch scheme
    
    case 'within',
        for i=1:size(X,1)
            didfit = 0;
            for j=1:size(div,1)
                if X(i,1)>=div(j,1) && X(i,2)<=div(j,2)
                    f = f+1;
                    fits(f,:) = X(i,:);
                    fix(f,:) = [i,j];
                    nix(j) = 0;
                    didfit = 1;
                end
            end  % end loop over divisions
            if ~didfit
                nf = nf+1;
                nonfits(nf,:) = X(i,:);
            end                
        end % end loop over intervals
        
    case 'begin',
        for i=1:size(X,1)
            didfit = 0;
            for j=1:size(div,1)
                if X(i,1)>=div(j,1) && X(i,1)<=div(j,2)
                    f = f+1;
                    fits(f,:) = X(i,:);
                    fix(f,:) = [i,j];
                    nix(j) = 0;
                    didfit = 1;
                end
            end % end loop over divisions
            if ~didfit
                nf = nf+1;
                nonfits(nf,:) = X(i,:);
            end
        end % end loop over intervals
        
    case 'pct',
        for i=1:size(X,1)
            didfit = 0;
            for j=1:size(div,1)
                pctsize = (X(i,2)-X(i,1))*pct; % pct% the length of X interval i
                % length of X interval i that falls within division j
                amtindiv = min(X(i,2),div(j,2))-max(X(i,1),div(j,1));
                if amtindiv>=pctsize
                    f = f+1;
                    fits(f,:) = X(i,:);
                    fix(f,:) = [i,j];
                    nix(j) = 0;
                    didfit = 1;
                end
            end  % end loop over divisions
            if ~didfit
                nf = nf+1;
                nonfits(nf,:) = X(i,:);
            end
        end % end loop over intervals
        
end % end switch on scheme

% appropriately size outputs for return
fits = fits(1:f,:);
nonfits = nonfits(1:nf,:);
fix = fix(1:f,:);
nix = find(nix);

end