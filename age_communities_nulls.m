%%%% GenLouvain! for community detection on age data
function output = age_communities_nulls(k,Ain,inpath,outpath,p,t,ts,run_ts,...
            tosave,remove_miss,goRange,missing,tag)

%%%%%%%%%%%%%%%%%%%%%%%%% MODIFY THIS PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inpath = '/Users/kimberly/Documents/lifespan/NNs/'; % should end with /
% outpath = '/Users/kimberly/Documents/lifespan/CD/'; % should end with /
% subjects = 108;
% missing_subjs = [];
% p = 100;    % number of genlouvain optimizations to perform
% t = 3;      % number of slices in multislice (total)
% ts = 316;   % number of time samples in each slice
% run_ts = 316; % number of time samples in each functional run
% tosave = 1;
% remove_miss = 0;
% 
% % list of gamma,omega pairs to test (omega~=0 gives categorical multislice)
% goRange = [1.05, 0;
%           ];                                    
% 
% % ib = list of subjects to include in analysis      
% ib = removeval((1:subjects)',missing_subjs);
% ib = 1;
% disp(ib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set number of nodes (brain regions) & slices per run
    try 
        A = csvread([inpath 'ageNN' num2str(k) '_run' num2str(1) ...
                 '-' num2str(1) '_tw' num2str(ts) '.txt']);
    catch 
            output = [num2str(k) ' file missing'];
            return
    end
    n = size(A,1);
    oldN = n;
% this is not needed since "missing" is passed in
%     if remove_miss
%         missing = []; % list of missing regions
%         %missing = find_missing_files('/Users/kimberly/Google Drive/choking/',nregs,ib,nruns);
%     end
    TperR = floor(run_ts/ts);
   
    
for goix = size(goRange,1)   % loop over gamma/omega pairs
    gamma = goRange(goix,1);
    omega = goRange(goix,2);

%for k=1:numel(ib)  %loop over subjects
    %person = ib(k);
    person = k;
    disp(['SUBJECT ' num2str(person)]);
    disp(['GAMMA = ' num2str(gamma)]);
	disp(['OMEGA = ' num2str(omega)]);

    % ensure doesn't already exist
    if exist([outpath 'ageCD' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_run' num2str(1) '-' num2str(1) '_tw' ...
                  num2str(ts) tag '_Call.txt'],'file')
        output = [num2str(k) ' file exists'];
        return
    end

% this is not needed since we want returned matrices with size oldN
%     if remove_miss   % if removing nodes, switch to appropriate n
%       if missing(person,:,:)
%         n = oldN - numel(unique([missing(person,:,:),0])) + 1;
%       end
%     end
    
if ~omega   % does static community detection on each slice separately

    Call = zeros(t*p,oldN);  % orig. brain region ordering, (t*p)xn
    Cnewall = Call;
    Qall = zeros(t*p,1);
    Qnewall = Qall;
    Qconsall = zeros(t,1);
    Zall = zeros(t*p*(p-1)/2,1);
    Znewall = zeros(t*p*(p-1)/2,1);
    for T=1:t  % loops over slices
%         r = floor((T-1)/TperR)+1;
%         i = mod(T-1,TperR)+1;
%         slicefile = [inpath 'ageNN' num2str(person) '_run' num2str(r) ...
%                      '-' num2str(i) '_tw' num2str(ts) '.txt'];
        
        % original nxn (potentially sparse) adjacency matrix for genlouvain
        %A = csvread(slicefile);
        A = Ain{T};
        if remove_miss
            [A,new_ix] = remove_missings(A,squeeze(missing(T,:)));
            n = length(A);
        else
            A(isnan(A)) = 0;
        end     
        
        % does p optimizations of genlouvain optimization on A with value gamma
        [C,Q] = genlouvainREPs(A,p,gamma);
        [Zallmean,zdist] = zrand(C);  % calculates rand-Z values (avg and dist)
        m = reshape(triu(zdist,1),p*p,1);
        m = m(triu(ones(p),1)>0);
        Zall((T-1)*p*(p-1)/2+1:T*p*(p-1)/2) = m; % saves rand-Z distribution
                
        % community consensus
        [Cnew,Qnew,~,Qcons] = consensus_comm_GL2(C);  % uses "genlouvainREPs"
        if min(isnan(Zallmean))
            Znewall((T-1)*p*(p-1)/2+1:T*p*(p-1)/2) = NaN;
        else
            [~,zdistc] = zrand(Cnew); % calculates rand-Z values (avg and dist)
            m = reshape(triu(zdistc,1),p*p,1);
            m = m(triu(ones(p),1)>0);
            Znewall((T-1)*p*(p-1)/2+1:T*p*(p-1)/2) = m; % saves rand-Z distribution
        end
        
        % add missing rows of 0s back in
        if remove_miss
            C = add_missings_dim(C,new_ix,oldN,2,0);
            Cnew = add_missings_dim(Cnew,new_ix,oldN,2,0);
        end

        Call(p*(T-1)+1:p*T,:) = C;
        Cnewall(p*(T-1)+1:p*T,:) = Cnew;
        Qall(p*(T-1)+1:p*T) = Q;
        Qnewall(p*(T-1)+1:p*T) = Qnew;
        Qconsall(T) = Qcons;
                
    end  % end loop over slices
    
else 
    % multislice!!! 
    %A = cell(t,1);
    A = Ain;
    for T=1:t
%         r = floor((T-1)/TperR)+1;
%         i = mod(T-1,TperR)+1;
%         slicefile = [inpath 'ageNN' num2str(person) '_run' num2str(r) ...
%                      '-' num2str(i) '_tw' num2str(ts) '.txt'];
%         A{T} = csvread(slicefile);
        if remove_miss
            [A{T},new_ix] = remove_missings(A{T},unique(missing));
            n = length(A{T});
        else
            A{T}(isnan(A{T})) = 0;
        end
    end
    [Cmult,Qall] = genlouvainREPs(A,p,gamma,omega); % output is pxnxt
    [~,zdist] = zrand(reshape(Cmult,p,n*t));
    m = reshape(triu(zdist,1),p*p,1);
    m = m(triu(ones(p),1)>0);
    Zall = m;
    % reshape into (t*p)xn for plotting
    Call = zeros(t*p,n);
    for T=1:t
        Call((T-1)*p+1:T*p,:) = Cmult(:,:,T);
    end
    
    % community consensus: input matrix must be px(n*t)
    [Cmultnew,Qnewall,~,Qconsall] = consensus_comm_GL2(reshape(Cmult,p,n*t));  % uses "genlouvainREPs"
    if isnan(Zall)
        Znewall(1:p*(p-1)/2) = NaN;
    else
        [~,zdistc] = zrand(Cmultnew);
        m = reshape(triu(zdistc,1),p*p,1);
        m = m(triu(ones(p),1)>0);
        Znewall = m;
    end
    Cmultnew = reshape(Cmultnew,p,n,t); % this is pxnxt
    % passes out (t*p)xn matrix for plotting
    for T=1:t
        Cnewall((T-1)*p+1:T*p,:) = Cmultnew(:,:,T);
    end
    clear Cmult Cmultnew;
    
    % add missing rows with 0s into Call, Cnewall
    if remove_miss
        Call = add_missings_dim(Call,new_ix,oldN,2,0);
        Cnewall = add_missings_dim(Cnewall,new_ix,oldN,2,0);
    end
    
end  % end if-else switching between static and multislice

    if tosave
        % save orig. p partitions in '_Call'
        dlmwrite([outpath 'ageCD' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_tw' num2str(ts) tag '_Call.txt'],Call);
        
        % save community consensus partitions in '_Cnewall'
        dlmwrite([outpath 'ageCD' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_tw' num2str(ts) tag '_Cnewall.txt'],Cnewall);
        
        % save Qs and CC Qs in '_Qs'
        dlmwrite([outpath 'ageCD' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_tw' num2str(ts) tag '_Qs.txt'],[Qall,Qnewall]);
        
        % save quality of community consensus in '_Qconsall'
        dlmwrite([outpath 'ageCD' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_tw' num2str(ts) tag '_Qconsall.txt'],Qconsall);
        
        % save rand-Z dists and CC rand-Z dists in '_Zs'
        dlmwrite([outpath 'ageCD' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_tw' num2str(ts) tag '_Zs.txt'],[Zall,Znewall]);  
    end
    
% end  % end loop over individuals

end % end loop over gamma/omega pairs

output = [num2str(k) ' ' num2str(gamma) ' ' num2str(omega) ' matrix created'];

end %function

% example A from genlouvain netwiki site
% A = [0,1,1,0,0,0,0,0,0,1;
%      1,0,1,0,0,0,0,0,0,0;
%      1,1,0,0,0,0,0,0,0,0;
%      0,0,0,0,1,1,0,0,0,1;
%      0,0,0,1,0,1,0,0,0,0;
%      0,0,0,1,1,0,0,0,0,0;
%      0,0,0,0,0,0,0,1,1,1;
%      0,0,0,0,0,0,1,0,1,0;
%      0,0,0,0,0,0,1,1,0,0;
%      1,0,0,1,0,0,1,0,0,0];

% example A: small world network (n=100, k=3, phi_=0.2, isperiodic)
% A = smallworld(100,3,0.2,1);