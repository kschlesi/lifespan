%%%% GenLouvain! for community detection on age data

%%%%%%%%%%%%%%%%%%%%%%%%% MODIFY THIS PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inpath = '/Users/kimberly/Documents/lifespan/NNs/'; % should end with /
outpath = '/Users/kimberly/Documents/lifespan/CD/'; % should end with /
subjects = 108;
missing_subjs = [46;64;81;82];
p = 100;    % number of genlouvain optimizations to perform
t = 3;      % number of slices in multislice (total)
ts = 316;   % number of time samples in each slice
run_ts = 316; % number of time samples in each functional run
tosave = 1;
remove_miss = 0;

% list of gamma,omega pairs to test (omega~=0 gives categorical multislice)
goRange = [1.15, 0.001;
          ];                                    

% ib = list of subjects to include in analysis      
ib = removeval((1:subjects)',missing_subjs);
ib = ib(1);
disp(ib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set number of nodes (brain regions) & slices per run
    A = csvread([inpath 'ageNN' num2str(ib(1)) '_run' num2str(1) ...
                 '-' num2str(1) '_tw' num2str(ts) '.txt']);
    n = size(A,1);
    if remove_miss
        missing = []; % list of missing regions
        %missing = find_missing_files('/Users/kimberly/Google Drive/choking/',nregs,ib,nruns);
        oldN = n;
    end
    TperR = floor(run_ts/ts);
   
    
for goix = size(goRange,1)   % loop over gamma/omega pairs
    gamma = goRange(goix,1);
    omega = goRange(goix,2);

for k=1:numel(ib)  %loop over subjects
    person = ib(k);
    disp(['SUBJECT ' num2str(person)]);
    disp(['GAMMA = ' num2str(gamma)]);
	disp(['OMEGA = ' num2str(omega)]);
    
    if remove_miss   % if removing nodes, switch to appropriate n
      if missing(person,:,:)
        n = oldN - numel(unique([missing(person,:,:),0])) + 1;
      end
    end
    
    % create cell array A for multislice CD
    A = cell(t,1);
    for T=1:t
        r = floor((T-1)/TperR)+1;
        i = mod(T-1,TperR)+1;
        slicefile = [inpath 'ageNN' num2str(person) '_run' num2str(r) ...
                     '-' num2str(i) '_tw' num2str(ts) '.txt'];
        A{T} = csvread(slicefile);
        if remove_miss
            A{T} = remove_missings(A{T},squeeze(missing(person,runnum,:)));
        else
            A{T}(isnan(A{T})) = 0;
        end
    end
    
    % multislice!!! and rand z score
    [Cmult,Qall] = genlouvainREPs_f(A,p,gamma,omega); % output is pxnxt
    [Zall1,zdist] = zrand(reshape(Cmult,p,n*t));
    m = reshape(triu(zdist,1),p*p,1);
    m = m(triu(ones(p),1)>0);
    Zall = m;
    
    % reshape into (t*p)xn for plotting
    Call = zeros(t*p,n);
    for T=1:t
        Call((T-1)*p+1:T*p,:) = Cmult(:,:,T);
    end
    
    % community consensus: input matrix must be pxnxt, output Cmultnew is pxnxt
    [Cmultnew,Qnewall,~,Qconsall] = consensus_comm_GLf(Cmult,gamma,omega);
    if isnan(Zall)
        Znewall(1:p*(p-1)/2) = NaN;
    else
        [Znewall1,zdistc] = zrand(reshape(Cmultnew,[p,n*t]));
        m = reshape(triu(zdistc,1),p*p,1);
        m = m(triu(ones(p),1)>0);
        Znewall = m;
    end
    
    % passes out (t*p)xn matrix for plotting
    for T=1:t
        Cnewall((T-1)*p+1:T*p,:) = Cmultnew(:,:,T);
    end
    clear Cmult Cmultnew;
    

    if tosave
        % save orig. p partitions in '_Call'
        dlmwrite([outpath '/ageCDf' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_run' num2str(r) '-' num2str(i) '_tw' ...
                  num2str(ts) '_Call.txt'],Call);
        
        % save community consensus partitions in '_Cnewall'
        dlmwrite([outpath '/ageCDf' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_run' num2str(r) '-' num2str(i) '_tw' ...
                  num2str(ts) '_Cnewall.txt'],Cnewall);
        
        % save Qs and CC Qs in '_Qs'
        dlmwrite([outpath '/ageCDf' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_run' num2str(r) '-' num2str(i) '_tw' ...
                  num2str(ts) '_Qs.txt'],[Qall,Qnewall]);
        
        % save quality of community consensus in '_Qconsall'
        dlmwrite([outpath '/ageCDf' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_run' num2str(r) '-' num2str(i) '_tw' ...
                  num2str(ts) '_Qconsall.txt'],Qconsall);
        
        % save rand-Z dists and CC rand-Z dists in '_Zs'
        dlmwrite([outpath '/ageCDf' num2str(gamma) num2str(omega) '_' ...
                  num2str(person) '_run' num2str(r) '-' num2str(i) '_tw' ...
                  num2str(ts) '_Zs.txt'],[Zall,Znewall]);  
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