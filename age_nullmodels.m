% code for null model tests...

inpath = '/Users/kimberly/Documents/lifespan/';

%goRange = [1.15, 0.005]; ts = 316;  stats_tag = 316;   % [gamma, omega] values and timesteps per window
%     goRange = [1.15,0.001]; ts = 316; stats_tag = '316alt1';
%     goRange = [1   ,0.005]; ts = 316; stats_tag = '316alt2';
%     goRange = [1   ,0.001]; ts = 316; stats_tag = '316alt3';
goRange = [1.2, 0.05]; ts = 52; stats_tag = '52';
s_plot = 0;             % list of subjects to plot (out of all subjs with both CD and behavior)

corr_type = 'Spearman'; % type of correlation (Spearman, Pearson, or Kendall)
is_runwise = 0;         % should correlations be done separately within each run?
inc_leg = 0;            % should a legend of all subjs be included?
remove_vis = 0;         % should tagged nodes be removed?
rm_tag = '_nv';         % sets tag for removing nodes (visual tag = '_nv')
reg_motion = 1;         % check correlations with motion partialled out?
extra_tag = '_nullspace';%'_norm1';   % if density normalization is used in CD

null_test = 1;
nullp = 100;
time_cont = 0;
netdiag = 0;

totalsubjs = 108;
missing_subjs = [46;64;81;82]; % final list
missing = []; remove_miss = 0;
tosave = 1;
% subject 64 is missing behavioral info only
ts_run = 316;
nruns = 3;

%% load subjects and demographics

% load subject IDs
%load('lifespan_subjs.mat');
load('age_timeseries.mat');
subj_IDs = cell(length(age_timeseries),1);
for k=1:length(age_timeseries)
    subj_IDs{k} = age_timeseries(k).ID;
end
clear age_timeseries;

% load file paths
NNinpath = [inpath 'NNs/'];
CDinpath = [inpath 'CD/'];
if remove_vis
    tag = rm_tag;
else
    tag = [];
end

% import age, perormance, and motion data
ages = xlsread([inpath 'Lifespan Behavioral Data for Sharing.xls'],'Sheet1','F2:F127');
sexs = xlsread([inpath 'Lifespan Behavioral Data for Sharing.xls'],'Sheet1','E2:E127');
dprimes = xlsread([inpath 'Lifespan Behavioral Data for Sharing.xls'],'Sheet1','EE2:EE127');
critss = xlsread([inpath 'Lifespan Behavioral Data for Sharing.xls'],'Sheet1','EF2:EF127');
[~,ageIDs,~] = xlsread([inpath 'Lifespan Behavioral Data for Sharing.xls'],'Sheet1','C2:C127');
vars = load('motion_extras.mat');
sAges = zeros(length(subj_IDs),1); sSexs = sAges;
sDPrime = sAges; rDPrime = zeros(length(subj_IDs),3);
sCSwitch = sAges; rCSwitch = zeros(length(subj_IDs),3);
sMotion = sAges; sGmin = sAges;
for i=1:length(subj_IDs)
    sAges(i) = ages(strcmp(ageIDs,subj_IDs(i)));
    sSexs(i) = sexs(strcmp(ageIDs,subj_IDs(i)));
    sDPrime(i) = dprimes(strcmp(ageIDs,subj_IDs(i)));
    sCSwitch(i) = critss(strcmp(ageIDs,subj_IDs(i)));
  if any(strcmp(vars.vars.IDs,subj_IDs(i)))
    sMotion(i) = vars.vars.motion(strcmp(vars.vars.IDs,subj_IDs(i)));
    sGmin(i) = vars.vars.gmin(strcmp(vars.vars.IDs,subj_IDs(i)));
  else
    sMotion(i) = nan;
    sGmin(i) = nan;
  end
end

%% Load atlas regions and known brain systems/classes
% load atlas classifications for n=194 hybrid atlas
n = 194;
fid = fopen('hybrid_labels_long.txt');
regName = textscan(fid,'%s',n,'Delimiter','\n');
regName = regName{:};
fclose(fid);
[regClass,regCN] = region_class_map(regName,inpath);
legClass = unique(regClass);
nC = numel(legClass);

%% NULL MODEL 2: load and shuffle CD results
% real distributions
[t,n,ib,Cplot,Cplotall,~,~,Ncvg] = load_CD_results(totalsubjs, missing_subjs,...
                      nruns, ts, ts_run, goRange, tag, CDinpath);
goix = 1;
p = size(Cplotall,2)/numel(ib);
TperR = floor(ts_run/ts);
% compute partitions, flexibilities, ncomms
partnS = zeros(n,t,numel(ib)); % find final partitions
flexS = zeros(n,numel(ib));   % compute flexibilities
ncommS = zeros(t+1,numel(ib));  % number of communities in each slice & overall
%ncomm1 = zeros(t+1,numel(ib));  % number of communities in each slice & overall
sings = zeros(t+1,numel(ib));   % number of singletons
scommS = cell(numel(ib),1);     % SIZE of each numbered community in each slice
for k=1:numel(ib)
    origs = squeeze(Cplotall(goix,p*(k-1)+1:p*k,:,:));
    partnS(:,:,k) = squeeze(mode(origs,1));
    flexS(:,k) = flexibility(partnS(:,:,k)',1);
	ncommS(end,k) = max(max(partnS(:,:,k)));
    alls = partnS(:,:,k);
    sings(k) = numel(find(hist(alls(:),1:1:ncommS(end,k))==1));
    %ncomm1(end,k) = ncommS(end,k) - sings;
    for T=1:t
        ncommS(T,k) = numel(removeval(unique(partnS(:,T,k)),0));
        sings(T,k) = numel(find(hist(partnS(:,T,k),1:1:ncommS(T,k))==1));
        %ncomm1(T,k) = ncommS(T,k) - sings;
    end
    scommS{k} = hist(partnS(:,:,k),1:ncommS(end,k));
end
ncomm1 = ncommS - sings;

maxSize = max(cellfun(@max,cellfun(@max,scommS,'UniformOutput',false)));
maxNum = max(ncommS(end,:));
slicehistarrayS = cellfun(@(x)hist(x,1:maxSize),scommS,'UniformOutput',false);
slicehistarrayI = cellfun(@(x)hist(x.*(x<=maxSize),1:maxSize),scommS,'UniformOutput',false);
slicehistS = cell2mat(cellfun(@(x)sum(x,2),slicehistarrayS,'UniformOutput',false)');
slicehistI = cell2mat(cellfun(@(x)sum(x,2),slicehistarrayI,'UniformOutput',false)');
shistS = cellfun(@(x)sum(x,2),scommS,'UniformOutput',false);
maxTotSize = max(cell2mat(cellfun(@max,shistS,'UniformOutput',false)));
skipA = 50;
shistS = cell2mat(cellfun(@(x)hist(x,1:skipA:maxTotSize),shistS,'UniformOutput',false));

mean_drop = @(x,dr,varargin) nanmean(nanify(x,dr),varargin{:});
std_drop = @(x,dr,varargin) nanstd(nanify(x,dr),varargin{:});
meanCommSize = cell2mat(cellfun(@mean,cellfun(@mean,scommS,'UniformOutput',false),'UniformOutput',false));
meanCommSize1 = cell2mat(cellfun(@mean,cellfun(@(x)mean_drop(x,1),scommS,'UniformOutput',false),'UniformOutput',false));
stdCommSize = cell2mat(cellfun(@(x)std(x,[],1),scommS,'UniformOutput',false));
stdCommSize1 = cell2mat(cellfun(@(x)std_drop(x,1,[],1),scommS,'UniformOutput',false));

% compute recruitments
nodeix = ones(n,1);
recruS = zeros(n,numel(ib));
    for k=1:numel(ib)
        partnsToUse = partnS(~~nodeix,:,k)';
        MA = mod_allegiance(partnsToUse,0);
        recruS(~~nodeix,k) = corecruitment(MA,regCN(~~nodeix));
    end
    % class-wise recruitments 
    clRecru = zeros(nC,numel(ib));
    for c=1:nC
        clRecru(c,:) = mean(recruS((~~(nodeix.*(regCN==c))),:));
    end
    

% null distributions & flexibilities
    disp('null distributions');
    partnRand = zeros(n,t,numel(ib),nullp);
    flexRand = zeros(n,numel(ib),nullp);   % compute flexibilities
    for j=1:nullp
     disp(j);
      for k=1:numel(ib)
        pr = randperm(n); % create random partitions
        for T=1:t
            if time_cont == 0;
                pr = randperm(n);       % create random partitions
            end
            partnRand(:,T,k,j) = partnS(pr,T,k);  % p random assignments of n nodes
        end
        flexRand(:,k,j) = flexibility(partnRand(:,:,k,j)',1);
      end
    end 
    
% null recruitments
recruRand = zeros(n,numel(ib),nullp);
clRecruRand = zeros(nC,numel(ib),nullp);
  for j=1:nullp
    disp(['null recruitments ' num2str(j) '/' num2str(nullp)]);
    for k=1:numel(ib)
        partnsToUse = partnRand(~~nodeix,:,k,j)';
        MA = mod_allegiance(partnsToUse,0);
        recruRand(~~nodeix,k,j) = corecruitment(MA,regCN(~~nodeix));
    end
    % class-wise recruitments 
    for c=1:nC
        clRecruRand(c,:,j) = mean(recruRand((~~(nodeix.*(regCN==c))),:,j));
    end
  end

%% NULL MODEL 2: plots.

 % we shuffle all partitions, keeping total community size and number
      % but destroying all else... 
      % now we ask, how much do community size and number control flex?
      aa=squeeze(mean(flexRand,1));
      pvals_randflex = zeros(size(aa,1),1);
      pvals_allflex = zeros(size(aa,1),1);
      for i=1:104; 
          pvals_randflex(i)=numel(find(mean(flexS(:,i))>aa(i,:)))/size(aa,2); 
          pvals_allflex(i)=numel(find(mean(flexS(:,i))>aa))/numel(aa); 
      end
      
   % plot rand flexes by individual
   figure; plot(aa');
           title('null flexibilities by subject');
           xlabel('random iteration'); ylabel('flexibility');
           legend('single individual');
   figure; plot(mean(flexS)); hold all;
           plot(aa);
           title('true and null flexibilities');
           xlabel('subject'); ylabel('flexibility');
           legend('true flexibility','null distribution');
           hold off;
   [flexsort,fix] = sort(mean(flexS));
   figure; plot(flexsort); hold all;
           plot(aa(fix,:));
           title('true and null flexibilities (sorted)');
           xlabel('subject'); ylabel('flexibility');
           legend('true flexibility','null distribution');
           hold off;
   
   % make and hist these.      
      corrRs = zeros(nullp,1);
      corrPs = zeros(nullp,1);
      for i=1:nullp
      [corrRs(i),corrPs(i)] = correlate(sAges(ib),aa(:,i),'NoFigure',...
          'type','Spearman','partial',sMotion(ib));
      end      
    
      figure; hist(corrRs,nullp/2); % true = 0.3908 Pearson, 0.53 SpearmanMP
      hold on; 
              plot([0.53,0.53],[0,max(hist(corrRs,nullp/2))],'-r');
              title('distribution of null age-flexibility correlations (true \rho = 0.53)');
              xlabel('Spearman''s \rho'); ylabel('number of null instances');
      figure; hist(corrPs,nullp/2); % true = 4.1e-5 Pearson, 8.7e-9 SpearmanMP
      hold on; 
              plot([8.7e-9,8.7e-9],[0,max(hist(corrPs,nullp/2))+2],'-r');
              title('distribution of null age-flexibility p-values (true p = 8.9e-9)');
              xlabel('Spearman''s \rho'); ylabel('number of null instances');
    figure; hist(corrP1s(corrP1s<0.05),sum(corrP1s<0.05)/2); % true = 0.003 SpearmanMP
      hold on; 
              plot([8.7e-9,8.7e-9],[0,max(hist(corrP1s(corrP1s<0.05),sum(corrP1s<0.05)/2))+2],'-r');
              title('distribution of null age-flexibility p-values (true p = 0.003)');
              xlabel('Spearman''s \rho'); ylabel('number of null instances');
                            
              
              
%% still with null model 2:
% how much do community number and size control recruitments?

% SUBJECT recruitment
[sorted_rec,recix] = sort(mean(recruS));
figure; plot(sorted_rec,'-'); hold on;
        plot(squeeze(mean(recruRand(:,recix,:),1)),'--'); hold on;
        %plot(mean(recruS),'-k'); hold on;
        %plot(squeeze(mean(recruRand,1)),'--'); hold on;
        xlabel('subject'); ylabel('mean recruitment over nodes');
        title('subject-wise recruitment v. null');
        legend('recruitment','null distribution');
        
figure; plot(sorted_rec,'-'); hold on;
        nullsort = sort(squeeze(mean(recruRand(:,recix,:),1)),2);
        plot(nullsort(:,floor(nullp*0.95)));
        xlabel('subject'); ylabel('mean recruitment over nodes');
        title('subject-wise recruitment v. null');
        legend('recruitment','null distribution');
        
% NODE recruitment
figure; plot(mean(recruS,2),'-k'); hold on;
        plot(squeeze(mean(recruRand,2)),'--'); hold on;
        xlabel('node'); ylabel('mean recruitment over subjects');
        title('node-wise recruitment v. null');
        legend('recruitment','null distribution');
        
% CLASSWISE recruitment
figure; for i=1:nC
            subplot(3,4,i);
            [~,six] = sort(clRecru(i,:));
            plot(clRecru(i,six),'-k'); hold on;
            plot(squeeze(clRecruRand(i,six,:)),'--'); hold on;
            xlabel('subject'); ylabel('class recruitment');
            %legend('recruitment','null distribution');
            title(legClass(i));
        end
        
figure; for i=1:nC
            subplot(3,4,i);
            [~,six] = sort(clRecru(i,:));
            plot(clRecru(i,six)); hold on;
            nullsort = sort(squeeze(clRecruRand(i,six,:)),2);
            plot(nullsort(:,floor(nullp*0.95)));
            legend('recruitment','95% confidence','Location','NorthWest');
            title(legClass(i));
        end

%% NULL MODEL 1: load and shuffle connectivities

% load adjacency matrices
A = load_adjs(t,ib,subj_IDs,ts,ts_run,NNinpath);
total_conn = zeros(numel(ib),1);
      for ii=1:numel(ib)
          adjs = A(ii).adj;
          total_conn(ii) = 0;
          for T=1:size(adjs,1)
              total_conn(ii) = total_conn(ii) + sum(sum(adjs{T}));
          end
      end

% load all connectivity matrices
p = 100;
[t,n,ib] = load_CD_results(totalsubjs,missing_subjs,nruns,ts,ts_run,goRange,tag,CDinpath);
A = load_adjs(t,ib,subj_IDs,ts,ts_run,NNinpath);
ix = find(triu(ones(n),1));     % find triu indices

% for each matrix, shuffle all connection weights uniformly at random; this
% preserves the overall connectivity but destroys other community structure
% FIRST, shuffle ONLY in space, not time.
for P=1:nullp
  disp(['iteration ' num2str(P)]);
  Anull = A;  
  for k=1:numel(ib)
    % shuffle, ENFORCE SYMMETRY (i.e. shuffle triu, then reflect over diag)
    shix = ix(randperm(n*(n-1)/2)); % shuffle triu indices
    for T=1:t
        adjmat = zeros(n);
        adjmat(ix) = A(k).adj{T}(shix);     % shuffle adj triu entries
        Anull(k).adj{T} = adjmat + adjmat'; % place in new Anull struct
        assert(all(diag(adjmat)==0));
    end
  end
  %save([NNinpath 'Anull' tag extra_tag num2str(P) '.mat'],'Anull');
  clear Anull;
end

%%
% for each shuffled subject, re-compute the community structure from the
% multislice shuffled connectivity matrix... compute stats... and save.
%%%%%%%%%%%%%%%%%%% RUN ON CLUSTER!! %%%%%%%%%%%%%%%%%%%%%

%for P=1:nullp
  disp(['iteration ' num2str(P)]);
  tic;
  load([NNinpath 'Anull' tag extra_tag num2str(P) '.mat']);
  for k=1:numel(ib)
    output = age_communities_nulls(ib(k),Anull(k).adj,...
                            NNinpath,CDinpath,p,t,ts,ts_run,...
                            tosave,remove_miss,goRange,missing,...
                            [tag extra_tag num2str(P)]);
    disp(output);
  end
  [~,~,~,~,Cplotall] = load_CD_results(1,[],nruns,...%totalsubjs,missing_subjs,nruns,...
                       ts,ts_run,goRange,[tag extra_tag num2str(P)],CDinpath);
  if size(Cplotall,1)==1
      goix = 1;
  end                      
  Astats = struct('partnRand',cell(numel(ib),1),...
                  'flexRand',cell(numel(ib),1),...
                  'ncommRand',cell(numel(ib),1)...
                  );
  partnRand = zeros(n,t,numel(ib)); % find final partitions
  flexRand = zeros(n,numel(ib));    % compute flexibilities
  ncommRand = zeros(t+1,numel(ib)); % number of communities in each slice & overall
  for k=1:numel(ib)
      origs = squeeze(Cplotall(goix,p*(k-1)+1:p*k,:,:));
      partnRand(:,:,k) = squeeze(mode(origs,1));
      flexRand(:,k) = flexibility(partnRand(:,:,k)',1);
      for T=1:t
          ncommRand(T,k) = numel(removeval(unique(partnRand(:,T,k)),0));
      end
      ncommRand(end,k) = max(max(partnRand(:,:,k)));
      Astats(k).partnRand = partnRand;
      Astats(k).flexRand = flexRand;
      Astats(k).ncommRand = ncommRand;
  end
  %save([NNinpath 'ageRandStats' tag extra_tag num2str(P) '.mat'],'Astats');
  clear Astats;
  toc;



%% Null Model 1: Analysis
% load in data
% ib_orig = ib;
% msg = load('msg_subj_runs.txt');
% msg_subj = unique(msg(:,1));
% ib = removeval(ib,msg_subj);
% sub_ix = ismember(ib_orig,ib);
sub_ix = 1:numel(ib);
sub_find = find(sub_ix);
ncomm_Rand = zeros(numel(ib),p);
ncomm1_Rand = zeros(numel(ib),p);
flex_Rand = zeros(numel(ib),p);
scomm_Rand = cell(numel(ib),p);
msg = [];
for k=1:numel(ib)
  disp(['subject ' num2str(ib(k))]);
  for P=1:p
    try
        % load random attempt P
        partnR = load([CDinpath 'ageRandStats/ageRandStats' num2str(ib(k)) tag extra_tag num2str(P) '_partn.txt']);
        %ncommR = load([CDinpath 'ageRandStats/ageRandStats' num2str(ib(k)) tag extra_tag num2str(P) '_ncomm.txt']);
        flexR = load([CDinpath 'ageRandStats/ageRandStats' num2str(ib(k)) tag extra_tag num2str(P) '_flex.txt']);
        ncomm_Rand(k,P) = max(max(partnR));
        sings = numel(find(hist(partnR(:),1:1:ncomm_Rand(k,P))==1));
        ncomm1_Rand(k,P) = ncomm_Rand(k,P) - sings;
        %ncomm_Rand(k,P) = ncommR(end);
        flex_Rand(k,P) = mean(flexR);
        scomm_Rand{k,P} = hist(partnR,1:ncomm1_Rand(k,P));
    catch err    
       msg = [msg;ib(k),P];
    end
  end
end % end loop over ib

% shist_Rand = zeros(max(max(ncomm_Rand)),numel(ib));
% for k=1:numel(ib)
%   shist_hold = zeros(max(ncomm_Rand(k,:)),p);
%   for P=1:p  
%     shist_hold(:,P) = sum([scomm_Rand{k,P};zeros(max(ncomm_Rand(k,:))-ncomm_Rand(k,P),t)],2);
%   end
%   shist_Rand(:,k) = [sum(shist_hold,2);zeros(max(max(ncomm_Rand))-max(ncomm_Rand(k,:)),1)];
% end

maxSizeR = max(max(cellfun(@max,cellfun(@max,scomm_Rand,'UniformOutput',false))));
maxNumR = max(max(ncomm1_Rand));
slicehistarray_Rand = cellfun(@(x)hist(x,1:maxSizeR),scomm_Rand,'UniformOutput',false);
slicehist_Rand = cell2mat(  cellfun(@(y)shiftdim(y,-2),  cellfun(@(x)sum(x,2),slicehistarray_Rand,'UniformOutput',false) ,'UniformOutput',false) );
shist_Rand = cellfun(@(x)sum(x,2),scomm_Rand,'UniformOutput',false);
maxTotSizeR = max(max(cell2mat(cellfun(@max,shist_Rand,'UniformOutput',false))));
skipR = 20;
shist_Rand = cell2mat( cellfun(@(y)shiftdim(y,-1), cellfun(@(x)hist(x,1:skipR:maxTotSizeR),shist_Rand,'UniformOutput',false),'UniformOutput',false) ) ;

mean_drop = @(x,dr,varargin) nanmean(nanify(x,dr),varargin{:});
std_drop = @(x,dr,varargin) nanstd(nanify(x,dr),varargin{:});
meanCommSize_Rand = cell2mat(cellfun(@mean,cellfun(@mean,scomm_Rand,'UniformOutput',false),'UniformOutput',false));
meanCommSize1_Rand = cell2mat(cellfun(@mean,cellfun(@(x)mean_drop(x,1),scomm_Rand,'UniformOutput',false),'UniformOutput',false));
stdCommSize_Rand = cell2mat(cellfun(@(x)std(x,[],1),scomm_Rand,'UniformOutput',false));
stdCommSize1_Rand = cell2mat(cellfun(@(x)std_drop(x,1,[],1),scomm_Rand,'UniformOutput',false));

%% null model 1: PLots
% null distributions for community number
figure; plot(ncomm1(end,sub_ix),'k'); hold all;
        plot(ncomm1_Rand,'--');
        title('number of communities');
        xlabel('subject ID'); ylabel('total number of communities');
        legend('real data','null distribution','(preserving only','overall connectivity)');
[sorted_ncommS,ncomix] = sort(ncomm1(end,sub_ix));
figure; plot(sorted_ncommS,'k'); hold all;
        plot(ncomm1_Rand(ncomix,:),'--');        
        title('number of communities (sorted by subject total)');
        xlabel('subject ID'); ylabel('total number of communities');
        legend('real data','null distribution','(preserving only','overall connectivity)');
[~,tconix] = sort(total_conn(sub_ix));
figure; plot(ncomm1(end,tconix),'k'); hold all;
        plot(ncomm1_Rand(tconix,:),'--');        
        title('number of communities (sorted by subject total)');
        xlabel('subject ID'); ylabel('total number of communities');
        legend('real data','null distribution','(preserving only','overall connectivity)');

figure; plot(total_conn(tconix),ncomm1(end,tconix),'k'); hold all;
        plot(total_conn(tconix),ncomm1_Rand(tconix,:),'--');        
        title('number of communities (sorted by overall connectivity)');
        xlabel('OC'); ylabel('total number of communities');
        legend('real data','null distribution','(preserving only','overall connectivity)');

        
correlate(ncomm1(end,sub_ix),mean(ncomm1_Rand,2));
correlate(ncomm1(end,sub_find([1:27,29:34,36:end])),mean(ncomm1_Rand([1:27,29:34,36:end],:),2));
title('true number communities v. mean null value')
xlabel('true number of communities'); ylabel('mean number of communities over 100 null networks');

correlate(ncommS(end,sub_find([1:27,29:34,36:end])),var(ncomm1_Rand([1:27,29:34,36:end],:),[],2));
title('true var number communities v. var null value')
xlabel('true number of communities'); ylabel('var in number of communities over 100 null networks');


% correlations of ncomm with age (null and real)
agevec = sAges(ib_orig);
motionvec = sMotion(ib_orig);
correlate(agevec,ncommS(end,:),'type','Spearman','partial',motionvec);
title('age v. number of communities');
xlabel('age'); ylabel('number of communities');
correlate(agevec(sub_ix),mean(ncomm1_Rand,2),'type','Spearman','partial',motionvec(sub_ix));
title('age v. number of communities (NULL)');
xlabel('age'); ylabel('number of communities (null)');
correlate(agevec(sub_find([1:27,29:34,36:end])),mean(ncomm1_Rand([1:27,29:34,36:end],:),2),...
    'type','Spearman','partial',motionvec(sub_find([1:27,29:34,36:end])));
title('age v. number of communities (NULL), no outliers');
xlabel('age'); ylabel('number of communities (null)');

   % make and hist these.      
      corrR1s = zeros(nullp,1);
      corrP1s = zeros(nullp,1);
      for i=1:nullp
      [corrR1s(i),corrP1s(i)] = correlate(sAges(ib),ncomm1_Rand(:,i),'NoFigure',...
          'type','Spearman','partial',sMotion(ib));
      end      
    
      figure; hist(corrR1s,nullp/2); % true = 0.29 SpearmanMP
      hold on; 
              plot([0.2885,0.2885],[0,max(hist(corrR1s,nullp/2))],'-r');
              title('distribution of null age-flexibility correlations (true \rho = 0.29)');
              xlabel('Spearman''s \rho'); ylabel('number of null instances');
      figure; hist(corrP1s,nullp/2); % true = 0.003 SpearmanMP
      hold on; 
              plot([0.003,0.003],[0,max(hist(corrP1s,nullp/2))+2],'-r');
              title('distribution of null age-flexibility p-values (true p = 0.003)');
              xlabel('Spearman''s \rho'); ylabel('number of null instances');
      figure; hist(corrP1s(corrP1s<0.05),sum(corrP1s<0.05)/2); % true = 0.003 SpearmanMP
      hold on; 
              plot([0.003,0.003],[0,max(hist(corrP1s(corrP1s<0.05),sum(corrP1s<0.05)/2))+2],'-r');
              title('distribution of null age-flexibility p-values (true p = 0.003)');
              xlabel('Spearman''s \rho'); ylabel('number of null instances');
              


% random and actual bargraphs of community size distributions

% TOTAL community size over all windows
figure; bar(squeeze(mean(shist_Rand,2))');
        xlabel('community size S'); ylabel('number of communities with size S');
                set(gca,'XTick',0:1:numel(0:skipR:maxTotSizeR)-1,'XTickLabel',0:skipR:maxTotSizeR);
        title('total dynamic community size, NULL');
figure; bar(shistS');
        set(gca,'XTick',0:1:numel(0:skipA:maxTotSize)-1,'XTickLabel',0:skipA:maxTotSize);
        title('total dynamic community size');
        
% single window community size distributions, summed over windows
figure; bar(squeeze(mean(slicehist_Rand,2))');
        xlabel('community size S'); ylabel('number of communities with size S');
        title('single-window community size distributions, NULL');
figure; bar(slicehistS);
        xlabel('community size S'); ylabel('number of communities with size S');
        title('single-window community size distributions');
figure; bar(slicehistI); ylim([0,50]);
        xlabel('community size S'); ylabel('number of communities with size S');
        title('single-window community size distributions');