%%%% GenLouvain! for community detection on age data

%%%%%%%%%%%%%%%%%%%%%%%%% MODIFY THIS PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inpath = '/Users/kimberly/Documents/lifespan/NNs/'; % should end with /
outpath = '/Users/kimberly/Documents/lifespan/CD/'; % should end with /
subjects = 108;
missing_subjs = [];
p = 100;    % number of genlouvain optimizations to perform
t = 3;      % number of slices in multislice (per run)
ts = 316;   % number of time samples in each slice
tosave = 1;
remove_miss = 0;

% list of gamma,omega pairs to test (omega~=0 gives categorical multislice)
goRange = [1.05, 0;
          ];                                    

ib = removeval((1:subjects)',missing_subjs);
ib = 1;
disp(ib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set number of nodes (brain regions) & timesteps per run
    A = csvread([inpath 'ageNN' num2str(ib(1)) '_run' num2str(1) ...
                 '-' num2str(1) '_tw' num2str(ts) '.txt']);
    n = size(A,1);
    if remove_miss
        missing = []; % list of missing regions
        %missing = find_missing_files('/Users/kimberly/Google Drive/choking/',nregs,ib,nruns);
        oldN = n;
    end
   
    
for goix = size(goRange,1)
    gamma = goRange(goix,1);
    omega = goRange(goix,2);

for k=1:numel(ib)
    person = ib(k);
    disp(['SUBJECT ' num2str(person)]);
    disp(['GAMMA = ' num2str(gamma)]);
	disp(['OMEGA = ' num2str(omega)]);
    
    if remove_miss
      if missing(person,:,:)
        n = oldN - numel(unique([missing(person,:,:),0])) + 1;
      end
    end
    
if ~omega   % does static community detection on each slice separately

    Call = zeros(t*p,n);  % orig. brain region ordering, (t*p)xn
    Cnewall = Call;
    Qall = zeros(t*p,1);
    Qnewall = Qall;
    Qconsall = zeros(t,1);
    Zall = zeros(t*p*(p-1)/2,1);
    Znewall = zeros(t*p*(p-1)/2,1);
    for T=1:t  % loops over slices (fix this!!!!)
        %figure out how to open A 
        % should be a function of T and ts
        
        % original nxn (potentially sparse) adjacency matrix for genlouvain
        A = csvread(['/Users/kimberly/Desktop/choking_NNs/NN' num2str(person+100) ...
            '_' num2str(nregs) scheme num2str(runnum) code '.txt']);
        if remove_miss
            A = remove_missings(A,squeeze(missing(person,runnum,:)));
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
            [Znewallmean,zdistc] = zrand(Cnew); % calculates rand-Z values (avg and dist)
            m = reshape(triu(zdistc,1),p*p,1);
            m = m(triu(ones(p),1)>0);
            Znewall((T-1)*p*(p-1)/2+1:T*p*(p-1)/2) = m; % saves rand-Z distribution
        end

        Call(p*(T-1)+1:p*T,:) = C;
        Cnewall(p*(T-1)+1:p*T,:) = Cnew;
        Qall(p*(T-1)+1:p*T) = Q;
        Qnewall(p*(T-1)+1:p*T) = Qnew;
        Qconsall(T) = Qcons;
                
    end  % end loop over slices
    
else 
    % multislice!!! 
    A = cell(t,1);
    for T=1:t
        runnum = T;
        A{T} = csvread(['/Users/kimberly/Desktop/choking_NNs/NN' num2str(person+100) ...
            '_' num2str(nregs) scheme num2str(runnum) code '.txt']);
        if remove_miss
            A{T} = remove_missings(A{T},squeeze(missing(person,runnum,:)));
        else
            A{T}(isnan(A{T})) = 0;
        end
    end
    [Cmult,Qall] = genlouvainREPs(A,p,gamma,omega); % output is pxnxt
    [Zall1,zdist] = zrand(reshape(Cmult,p,n*t));
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
        [Znewall1,zdistc] = zrand(Cmultnew);
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
    
end  % end if-else switching between static and multislice

    if tosave
        % save orig. p partitions in '_Call'
        dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
            num2str(person+100) '/run' num2str(runnum) '/' scheme '/CDr' num2str(gamma) num2str(omega) ...
            '_' num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code '_Call.txt'],Call);
        
        % save community consensus partitions in '_Cnewall'
        dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
            num2str(person+100) '/run' num2str(runnum) '/' scheme '/CDr' num2str(gamma) num2str(omega) ...
            '_' num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code '_Cnewall.txt'],Cnewall);
        
        % save Qs and CC Qs in '_Qs'
        dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
            num2str(person+100) '/run' num2str(runnum) '/' scheme '/CDr' num2str(gamma) num2str(omega) ...
            '_' num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code '_Qs.txt'],[Qall,Qnewall]);
        
        % save quality of community consensus in '_Qconsall'
        dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
            num2str(person+100) '/run' num2str(runnum) '/' scheme '/CDr' num2str(gamma) num2str(omega) ...
            '_' num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code '_Qconsall.txt'],Qconsall);
        
        % save rand-Z dists and CC rand-Z dists in '_Zs'
        dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
            num2str(person+100) '/run' num2str(runnum) '/' scheme '/CDr' num2str(gamma) num2str(omega) ...
            '_' num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code '_Zs.txt'],[Zall,Znewall]);  
    end
    
end  % end loop over individuals

end % end loop over gamma/omega pairs

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