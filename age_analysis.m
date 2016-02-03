%%% lifespan data analysis.
%%% TO GENERATE FIGURES for lifespan paper, select desired timewindow
%%% from line 9/10, and include 'figureN' as an entry in p_plot on line 13,
%%% where N is the figure number (e.g. 'figure2' or 'figure6a'); then run.
%%
%%%%%%%%%%%%%%%%%%%%%%%%   modify parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
inpath = '/Users/kimberly/Documents/lifespan/';

%goRange = [1.15, 0.005]; ts = 316;  stats_tag = 316;   % [gamma, omega] values and timesteps per window
%     goRange = [1.15,0.001]; ts = 316; stats_tag = '316alt1';
%     goRange = [1   ,0.005]; ts = 316; stats_tag = '316alt2';
%     goRange = [1   ,0.001]; ts = 316; stats_tag = '316alt3';
goRange = [1.2, 0.05]; ts = 52; stats_tag = '52';
s_plot = 0;%,7,9];                   % list of subjects to plot (out of all subjs with both CD and behavior)
c_plot = {'Fronto-parietal'};          % list of classes to plot (when necessary)
p_plot = {};
% possible member functions of p_plot:
% comm flex ncomm corr_perfage cohs recru perfcohs adjs
% corr_age... _flex _ncomm _cohsrecru _cflex _cncomm _ccR
% corr_perf... _flex _ncomm _cohsrecru _cflex _cncomm _ccR
% class ccomm cncomm cflex crecru ccohs ccR
% figure1 figure2 figure3 figure4b figure5 figure6 figure7
% save_stats

corr_type = 'Spearman'; % type of correlation (Spearman, Pearson, or Kendall)
is_runwise = 0;         % should correlations be done separately within each run?
inc_leg = 0;            % should a legend of all subjs be included?
remove_vis = 0;         % should tagged nodes be removed?
rm_tag = '_nv';         % sets tag for removing nodes (visual tag = '_nv')
reg_motion = 1;         % check correlations with motion partialled out?

null_test = 0;
nullp = 100;
netdiag = 0;

totalsubjs = 108;
missing_subjs = [46;64;81;82]; % final list
% subject 64 is missing behavioral info only
ts_run = 316;
nruns = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   LOAD EVERYTHING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load subject IDs
%load('lifespan_subjs.mat');
load('age_timeseries.mat');
subj_IDs = cell(length(age_timeseries),1);
for k=1:length(age_timeseries)
    subj_IDs{k} = age_timeseries(k).ID;
end
totalBOLD = zeros(size(age_timeseries));
for k=1:length(age_timeseries)
    totalBOLD(k) = sum(sum(age_timeseries(k).run1))./size(age_timeseries(k).run1,1) + ...
                    sum(sum(age_timeseries(k).run2))./size(age_timeseries(k).run2,1) + ...
                    sum(sum(age_timeseries(k).run3))./size(age_timeseries(k).run3,1) ;
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
critr1 = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','H2:H127');
critr2 = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','M2:M127');
critr3 = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','R2:R127');
dprimer1 = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','D2:F127');
dprimer2 = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','I2:K127');
dprimer3 = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','N2:P127');
[~,ageIDsR,~] = xlsread([inpath 'Lifespan Behavioral Data for 3 Runs.xls'],'MAIN','C2:C127');
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
    rDPrime(i,1) = mean(dprimer1(strcmp(ageIDsR,subj_IDs(i)),[1,3]));
    rDPrime(i,2) = mean(dprimer2(strcmp(ageIDsR,subj_IDs(i)),[1,3]));
    rDPrime(i,3) = mean(dprimer3(strcmp(ageIDsR,subj_IDs(i)),[1,3]));
    rCSwitch(i,1) = critr1(strcmp(ageIDsR,subj_IDs(i)));
    rCSwitch(i,2) = critr2(strcmp(ageIDsR,subj_IDs(i)));
    rCSwitch(i,3) = critr3(strcmp(ageIDsR,subj_IDs(i)));
  if any(strcmp(vars.vars.IDs,subj_IDs(i)))
    sMotion(i) = vars.vars.motion(strcmp(vars.vars.IDs,subj_IDs(i)));
    sGmin(i) = vars.vars.gmin(strcmp(vars.vars.IDs,subj_IDs(i)));
  else
    sMotion(i) = nan;
    sGmin(i) = nan;
  end
end

% load community detection results
[t,n,ib,Cplot,Cplotall,~,~,Ncvg] = load_CD_results(totalsubjs, missing_subjs,...
                      nruns, ts, ts_run, goRange, tag, CDinpath);
goix = 1;
p = size(Cplotall,2)/numel(ib);
TperR = floor(ts_run/ts);
%ib = removeval(ib,[find(isnan(sDPrime)),find(isnan(sCSwitch))]);
nsubs = length(ib);
toUse = zeros(totalsubjs,1); toUse(ib) = 1; 
toUse(isnan(sDPrime)) = 0; toUse(isnan(sCSwitch)) = 0;
toUse1 = ones(nsubs,1); toUse1(isnan(sDPrime(ib))) = 0; toUse1(isnan(sCSwitch(ib))) = 0;

% load adjacency matrices
A = load_adjs(t,ib,subj_IDs,ts,ts_run,NNinpath);
%pplot = ib(kplot);

% compute partitions, flexibilities, ncomms
partnS = zeros(n,t,numel(ib)); % find final partitions
sureS = zeros(n,t,numel(ib),2);  % sureties (CHANGE THIS DUMB MEASURE)
flexS = zeros(n,numel(ib));   % compute flexibilities
ncommS = zeros(t+1,numel(ib));  % number of communities in each slice & overall
RflexS = zeros(n,nruns,nsubs);
RncommS = zeros(nruns+1,nsubs);
for k=1:numel(ib)
    origs = squeeze(Cplotall(goix,p*(k-1)+1:p*k,:,:));
    origsO = squeeze(Cplot(goix,p*(k-1)+1:p*k,:,:));
    partnS(:,:,k) = squeeze(mode(origs,1));
    flexS(:,k) = flexibility(partnS(:,:,k)',1);
    for T=1:t
        ncommS(T,k) = numel(removeval(unique(partnS(:,T,k)),0));
        %sureS(:,T,k,1) = sureties(origsO(:,:,T),partnS(:,T,k));
        %sureS(:,T,k,2) = sureties(origs(:,:,T),partnS(:,T,k));
    end
    for run=1:nruns
        RflexS(:,run,k) = flexibility(partnS(:,TperR*(run-1)+1:TperR*run,k)',1);
        RncommS(run,k) = numel(removeval(unique(partnS(:,TperR*(run-1)+1:TperR*run,k)),0));
    end
    ncommS(end,k) = max(max(partnS(:,:,k)));
    RncommS(end,k) = ncommS(end,k);
end

% generate random null distribution of partitions
if null_test
    disp('null distributions');
    partnRand = zeros(n,t,numel(ib),nullp);
    flexRand = zeros(n,numel(ib),nullp);   % compute flexibilities
    for j=1:nullp
     disp(j);
      for k=1:numel(ib)
        for T=1:t
            pr = randperm(n);       % create random partitions
            partnRand(:,T,k,j) = partnS(pr,T,k);  % p random assignments of n nodes
        end
        flexRand(:,k,j) = flexibility(partnRand(:,:,k,j)',1);
      end
    end 
end


%%% compute size of each community
cMaxMin = zeros(t+1,2,nsubs);
cSizesA = zeros(max(max(ncommS)),nsubs);
cSizesT = zeros(max(max(ncommS)),t,nsubs);
for j = 1:nsubs
    cvals = removeval(unique(partnS(:,:,j)),0);
    maxc = zeros(1,t);
    minc = n.*ones(1,t);
    maxcA = 0;
    mincA = n*t;
    for c = cvals'
        cSizes = sum(partnS(:,:,j)==c);
        cSizesT(c,:,j) = cSizes;
        cSizesA(c,j) = sum(cSizes);
        maxc = max([cSizes;maxc]);
        minc = min([cSizes;minc]);
        maxcA = max([sum(cSizes),maxcA]);
        mincA = min([sum(cSizes),mincA]);
    end
    cMaxMin(:,:,j) = [maxc, maxcA; minc, mincA]';            
end

% load atlas classifications for n=194 hybrid atlas
fid = fopen('hybrid_labels_long.txt');
regName = textscan(fid,'%s',n,'Delimiter','\n');
regName = regName{:};
fclose(fid);
[regClass,regCN] = region_class_map(regName,inpath);
legClass = unique(regClass);

% load physical positions of nodes
%%% note: these are for officers, not lifespan subjs
%%% use for visualization purposes ONLY
pos = load('hybrid_centroids.mat');  
pos = squeeze(mean(pos.hybrid_c,2));
pos = pos(ib,:,:);

% set node IDs to consider
if remove_vis
    nodeix = ~ismemvar(regClass,'Visual');
else
    nodeix = ones(n,1);
end

% compute cohesions (nee coherences)
nC = numel(unique(regCN(~~nodeix)));
are_cohs = 0;
if ismember('cohs',p_plot) || ismember('perfcohs',p_plot) ||...
   ismember('ccohs',p_plot) || ismember('corr_age_cohsrecru',p_plot) %||...
   %ismember('corr_perf_cohsrecru',p_plot)
    reps = 1000;
    subjmcohr = zeros(nsubs,t,nC);
    subjmcomm = zeros(nsubs,t,nC);
    avgmcohr = zeros(nsubs,t,nC);
    allmcohr = zeros(nsubs,t,reps,nC);
    mcohr_ps = zeros(nsubs,t,3,nC);
    %fskew_ps = zeros(nsubs,t,3,nC);
    for k=1:nsubs
        disp(['cohesion for ' num2str(ib(k))]);
      for c=unique(regCN(~~nodeix))'
        disp(legClass(c));
        % find cohesion per classified region
        [maxcohr,maxcomm,cc] = coherence(find(regCN(~~nodeix,:)==c),partnS(~~nodeix,:,k));
        [p1,~,p3,mcohr,~,fskew] = coherence_null(find(regCN(~~nodeix,:)==c),partnS(~~nodeix,:,k),reps,0);
        subjmcohr(k,:,c) = maxcohr';
        subjmcomm(k,:,c) = maxcomm';
        avgmcohr(k,:,c) = mean(mcohr,2);
        allmcohr(k,:,:,c) = shiftdim(mcohr,-1);
        mcohr_ps(k,:,:,c) = shiftdim(p1(1:t,:),-1);
      end
    end
    are_cohs = 1;
end

% compute recruitments
%if ismember('recru',p_plot) || ismember('crecru',p_plot)
    recruS = zeros(n,nsubs,nruns);
    corecruS = zeros(nC,nC,nsubs,nruns);
    for k=1:nsubs
      for r=1:nruns
        partnsToUse = partnS(~~nodeix,(r-1)*TperR+1:r*TperR,k)';
        MA = mod_allegiance(partnsToUse,0);
        recruS(~~nodeix,k,r) = corecruitment(MA,regCN(~~nodeix));
        corecruS(:,:,k,r) = class_corecruitment(MA,regCN(~~nodeix));
      end
    end
    % class-wise recruitments 
    clRecru = zeros(nC,nsubs,nruns);
    for c=1:nC
        clRecru(c,:,:) = mean(recruS((~~(nodeix.*(regCN==c))),:,:),1);
    end
%end
if null_test
  disp('null recruitments');
  recruRand = zeros(n,nsubs,nruns,nullp);
  clRecruRand = zeros(nC,nsubs,nruns,j);
  for j=1:nullp
    for k=1:nsubs
      for r=1:nruns
        partnsToUse = partnRand(~~nodeix,(r-1)*TperR+1:r*TperR,k,j)';
        MA = mod_allegiance(partnsToUse,0);
        recruRand(~~nodeix,k,r,j) = corecruitment(MA,regCN(~~nodeix));  
      end
    end
    % class-wise recruitments 
    disp(j);
    for c=1:nC
        clRecruRand(c,:,:,j) = mean(recruRand((~~(nodeix.*(regCN==c))),:,:,j),1);
    end
  end
end

% compute average flexibility per class
nCl = numel(unique(regCN(~~nodeix)));
clSize = zeros(nCl,1);
clFlex = zeros(nCl,nsubs);
RclFlex = zeros(nCl,nruns,nsubs);
for c=1:nC
    clSize(c) = sum(regCN(~~nodeix)==c);
    clFlex(c,:) = mean(flexS(~~(nodeix.*(regCN==c)),:));
    for run=1:nruns
        RclFlex(c,run,:) = mean(RflexS(~~(nodeix.*(regCN==c)),run,:),1);
    end
end
if null_test
  disp('null class flexibilities');
  clFlexRand = zeros(nCl,nsubs,nullp);
  for j=1:nullp
    for c=1:nC
      clFlexRand(c,:,j) = mean(flexRand(~~(nodeix.*(regCN==c)),:,j));
%       for run=1:nruns
%           RclFlex(c,run,:,j) = mean(RflexS(~~(nodeix.*(regCN==c)),run,:,j),1);
%       end
    end
  end
end


% compute number of communities in each class
cncommS = zeros(nC,nsubs,t+1);
for c=1:nC
  for k1=1:nsubs
    for T=1:t
      ncms = numel(removeval(unique(partnS(~~(nodeix.*(regCN==c)),T,k1)),0));
      cncommS(c,k1,T) = ncms;
    end
    ncall = numel(removeval(unique(partnS(~~(nodeix.*(regCN==c)),:,k1)),0));
    cncommS(c,k1,t+1) = ncall;
  end
end

% compute basic network diagnostics from adj
% calculate node degrees
if netdiag
overallND = zeros(n,nruns,nsubs);
inND = zeros(n,nruns,nsubs);
outND = zeros(n,nruns,nsubs);
cSize = zeros(n,nruns,nsubs);
for k=1:nsubs
  for T=1:t
    weights = A(k).adj{T}(~~nodeix,~~nodeix);
    sameC = mod_allegiance([partnS(~~nodeix,T,k),partnS(~~nodeix,T,k)]',0);
    overallND(~~nodeix,T,k) = sum(weights)';
    inND(~~nodeix,T,k) = sum(weights.*sameC)';
    outND(~~nodeix,T,k) = sum(weights.*~sameC)';
    cSize(~~nodeix,T,k) = sum(sameC)';
  end   
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE ALL PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=s_plot
  %if k~=0
    %posPass = squeeze(pos(k,:,~~nodeix))'; % Nx3 positions
  %else
    posPass = squeeze(pos(1,:,~~nodeix))';
  %end
  if ismember('save_stats',p_plot)
      rpsAgeAll = zeros(33,4);
      rpsDPAll = zeros(33,4);
      rpsCSAll = zeros(33,4);
      reg_motion = 1;
      fig_tag = 'NoFigure';
  else
      fig_tag = [];
  end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% PLOTS from chosen keywords %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % partitions, both layout and BrainView
        if ismember('comm',p_plot) && k~=0
        figure;
        bcolor(partnS(~~nodeix,:,k)); colorbar; 
        xlabel('run'); set(gca,'XTick',1.5:1:t+1.5-1,'XTickLabel',1:t);
        ylabel('brain region ID');
        title(['task-associated communities, subject ' subj_IDs(ib(k))]);
%         for T=1:t
%             figure;
%             title(['subject ' subj_IDs{ib(k)} ' communities, slice ' num2str(T)]);
%             BrainView(posPass,partnS(~~nodeix,T,k),0,'View','sagittal','ExtFigureCmd');
%         end      
        end
        
        % number of communities
        if (ismember('ncomm',p_plot) || ismember('figure3',p_plot)) && k==0
            if ~ismember('figure3',p_plot)
            figure; bcolor(ncommS'); colorbar;
            title('number of communities');
            xlabel('run'); set(gca,'XTick',1.5:1:t+1.5,'XTickLabel',1:t);
            ylabel('subject')
            end
            
            xstep = 10;
            [sortedcomms,six] = sortentry(ncommS,'row',0,sAges(ib));
            %figure; bcolor(sortedcomms'); colorbar;
            figure; bcolor(sortedcomms(1:end-1,:)'); colorbar;
            title('number of communities sorted by subject age');
            xlabel('run'); set(gca,'XTick',1.5:1:t+1.5,'XTickLabel',1:t);
            ylabel('subject age'); set(gca,'YTick',1:xstep:nsubs,'YTickLabel',sAges(ib(six(1:xstep:nsubs))));
        end
            
        if (ismember('corr_age_ncomm',p_plot) && k==0 ) || ismember('save_stats',p_plot)
            [rpsAgeAll(12,1),rpsAgeAll(12,2)] = correlate(sAges(ib),ncommS(end,:),'type','Spearman');
            title('total community number v. age');
            xlabel('age')
            ylabel('total community number')
            if reg_motion %%%%% regress motion
                [rpsAgeAll(12,3),rpsAgeAll(12,4)] = correlate(sAges(ib),ncommS(end,:),'partial',sMotion(ib),'type','Spearman','DispResiduals','Raw');
                title('total community number v. age, MC');
                xlabel('age')
                ylabel('total community number')
            end
            if is_runwise && t~=nruns
              for run=1:nruns
                correlate(sAges(ib),RncommS(run,:),'type',corr_type)
                title(['community number v. age, run ' num2str(run)]);
                xlabel('age')
                ylabel('total community number')
                if reg_motion
                correlate(sAges(ib),RncommS(run,:),'partial',sMotion(ib),'type',corr_type)
                title(['community number v. age, run ' num2str(run) ', MC']);
                xlabel('age')
                ylabel('total community number')
                end
              end
            end
        end
        
        % hmmm
        if ismember('scomm',p_plot) && k~=0
            cscale1 = max(max(partnS(:,:,k)));
            cmin = 0;
            if T==1
                figure;
                side = ceil(sqrt(length(taskcell)));
                suptitle(['subject ' subj_IDs{k} ' communities']);
              for T1 = 1:length(taskcell)
                max(partnS(:,T1,k))
                subplot(side,side,T1);
                BrainView(posPass,partnS(~~nodeix,T1,k),0,'View','sagittal','ExtFigureCmd');
                title(char(taskcell(T1)));
              end
            end
        end

        % example of particular classes (size) v. communities (color)
        if ismember('ccomm',p_plot) && k~=0
            for w=0
                for c=1:9
                    figure;
                    side=ceil(sqrt(length(taskcell)));
                    suptitle([char(subj_IDs(ib(k))) ' communities and ' char(unique(regClass(regCN==c))) ...
                         ' regions...'])
                    for T1=1:length(taskcell)
                        subplot(side,side,T1);
                        BrainView(posPass,partnS(~~nodeix,T1,k),0,'NodeSize',50*((regCN==c)*3+1),'View','sagittal','ExtFigureCmd');
                        title(['...in ' char(taskcell(T1)) ', w = ' num2str(w)]);
                    end
                end
            end
        end
        
        % general plot of colored classes
        if ismember('class',p_plot) && k~=0
            if k~=0
            figure;
            % pick out one thing by indices
           % pickix = strcmp(regClass,'Default mode');  % should be nx1, binary
           % pickix = ismemvar(strfind(regName,'amygdala'),3);
            pickix = zeros(n,1);
            BrainView(posPass,regCN(~~nodeix),0,'NodeSize',100*(pickix(~~nodeix)*3+1),'View','sagittal','ExtFigureCmd');
            title('locations of classified regions');
            theC = colorbar;
            set(theC,'YTicklabel',legClass);
            colormap jet
            end
        end
        
        % example of particular classes (size) v. communities (color)
        if ismember('cclas',p_plot) && k~=0
            for w=0
                for T1=1:length(taskcell)
                    figure;
                    side=ceil(sqrt(numel(unique(partnS(~~nodeix,T1,k)))));
                    ax = zeros(numel(unique(partnS(~~nodeix,T1,k))),1);
                for cm=1:numel(unique(partnS(~~nodeix,T1,k)))
                        if numel(unique(partnS(~~nodeix,T1,k)))<=6; ax(cm) = subplot(side,side-1,cm);
                        else ax(cm) = subplot(side,side,cm); end;
                        BrainView(posPass,regCN(~~nodeix),0,'NodeSize',50*((partnS(~~nodeix,T1,k)==cm)*3+1),'View','sagittal','ExtFigureCmd');
                        title(['...in comm. ' num2str(cm)]);
                        spp = get(ax(cm), 'Position');
                        set(ax(cm), 'Position', [spp(1)*0.85 spp(2) 0.85*spp(3) spp(4)]);
                end
                suptitle([char(subj_IDs(ib(k))) ' classified regions in ' char(taskcell(T1)) ', w = ' num2str(w)]);
                theC = colorbar;
                set(theC,'Position', [.8314 .11 .0581 .8150]);%,'Ytick',log(L));
                set(theC,'YTicklabel',legClass);
                colormap jet
                end
            end
        end
             
        
        % flexibility maps
        if ismember('flex',p_plot) || ismember('figure4b',p_plot)
          if k==0 
              if ~ismember('figure4b',p_plot)
              % create a mean flexibility per node figure
              BrainView(posPass,mean(flexS(~~nodeix,:),2),0,'View','sagittal');
              colorbar('Location','EastOutside');
              title('Mean Node Flexibility across Subjects');
              % create a variance in flexibility per node figure
              BrainView(posPass,var(flexS(~~nodeix,:),0,2),0,'View','sagittal');
              colorbar('Location','EastOutside');
              title('Variance in Node Flexibility across Subjects');
              end
              correlate(mean(flexS(~~nodeix,:),2),var(flexS(~~nodeix,:),0,2));
              title('mean v. variance in flexibility across subjects')
              xlabel('mean flexibility across subjects, for each node')
              ylabel('variance in flexibility across subjects')              
          else
              if ~ismember('figure4b',p_plot)
              BrainView(posPass,flexS(~~nodeix,k),0,'View','sagittal');
              colorbar('Location','EastOutside');
              title(['Node Flexibility, Subject ' char(subj_IDs(ib(k)))]);
              end
          end
        end
        
        % flexibility correlation with age
        if ismember('corr_age_flex',p_plot) || ismember('figure5',p_plot) || ismember('save_stats',p_plot)
          if k==0 
              [rpsAgeAll(1,1),rpsAgeAll(1,2)] = correlate(sAges(ib),mean(flexS(~~nodeix,:),1),'type','Spearman');
              title('mean flexibility over all nodes v. age')
              ylabel('mean flexibility over all nodes')
              xlabel('subject age')
              ttest2(mean(flexS(~~nodeix,sAges(ib)>40),1),mean(flexS(~~nodeix,sAges(ib)<40),1))
              legend('Location','SouthEast');
              if reg_motion %%%%%% regress motion              
                  [rpsAgeAll(1,3),rpsAgeAll(1,4)] = correlate(sAges(ib),mean(flexS(~~nodeix,:),1),'partial',sMotion(ib),'type','Spearman','DispResiduals','Raw');
                  title('mean flexibility over all nodes v. age, MC')
                  ylabel('mean flexibility over all nodes')
                  xlabel('subject age')
                  legend('Location','SouthEast');
              end
              if ~ismember('figure5',p_plot) && ~ismember('save_stats',p_plot)
              if is_runwise && t~=nruns
                for run=1:nruns
                  correlate(sAges(ib),mean(RflexS(~~nodeix,run,:),1),'type',corr_type);
                  title(['mean flexibility over all nodes v. age, run ' num2str(run)])
                  ylabel('mean flexibility over all nodes')
                  xlabel('subject age')
                  legend('Location','SouthEast');
                  if reg_motion %%%%%% regress motion              
                    correlate(sAges(ib),mean(RflexS(~~nodeix,run,:),1),...
                        'partial',sMotion(ib),'type',corr_type);
                      title(['mean flexibility over all nodes v. age, run ' num2str(run) ', MC'])
                      ylabel('mean flexibility over all nodes')
                      xlabel('subject age')
                      legend('Location','SouthEast');
                  end
                end
              end
              end
          end
        end
        
        % correlations of performance with age
        if ismember('corr_perfage',p_plot) && k==0
            % age
            figure;
            suptitle('Performance v. Age');
            subplot(1,2,1);
            correlate(sAges(~isnan(sDPrime)),sDPrime(~isnan(sDPrime)),'type',corr_type,'ExtFigureCmd');
            xlabel('age'); ylabel('dprime'); title('dprime v. age');
            subplot(1,2,2);
            correlate(sAges(~isnan(sCSwitch)),sCSwitch(~isnan(sCSwitch)),'type',corr_type,'ExtFigureCmd');
            xlabel('age'); ylabel('crit switch'); title('criterion switch score v. age');
        end
        
        if (ismember('corr_perf_flex',p_plot) && k==0) || ismember('save_stats',p_plot)
            %corr_type = 'Pearson'
            % flexibility
            flexToUse = mean(flexS(~~nodeix,:),1);
            dpToUse = sDPrime(ib);
            csToUse = sCSwitch(ib);
            motToUse = sMotion(ib);
            
            figure;
            suptitle('Performance v. Flexibility');
            subplot(1,2,1);
            [rpsDPAll(1,1),rpsDPAll(1,2)] = correlate(flexToUse,dpToUse,'type','Pearson','ExtFigureCmd');
            xlabel('flexibility'); ylabel('dprime'); title('dprime v. flexibility');
            subplot(1,2,2);
            [rpsCSAll(1,1),rpsCSAll(1,2)] = correlate(flexToUse,csToUse,'type','Pearson','ExtFigureCmd');
            xlabel('flexibility'); ylabel('crit switch'); title('crit shift score v. flexibility');
            if reg_motion
                figure;
                suptitle('Performance v. Flexibility, motion regressed');
                subplot(1,2,1);
                [rpsDPAll(1,3),rpsDPAll(1,4)] = correlate(flexToUse,dpToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd');
                xlabel('age'); ylabel('dprime'); title('dprime v. age, MC');
                subplot(1,2,2);
                [rpsCSAll(1,3),rpsCSAll(1,4)] = correlate(flexToUse,csToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd');
                xlabel('age'); ylabel('crit shift'); title('criterion shift score v. age, MC');
            end
            
            if is_runwise && t~=nruns
              for run=1:nruns
                flexToUse = mean(RflexS(~~nodeix,run,:),1);
                dpToUse = rDPrime(ib,run);
                csToUse = rCSwitch(ib,run);
                
                figure;
                suptitle(['Performance v. Flexibility, run ' num2str(run)]);
                subplot(1,2,1);
                correlate(flexToUse(~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                    'type',corr_type,'ExtFigureCmd');
                xlabel('flexibility'); ylabel('dprime'); title('dprime v. flexibility');
                subplot(1,2,2);
                correlate(flexToUse(~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                    'type',corr_type,'ExtFigureCmd');
                xlabel('flexibility'); ylabel('crit switch'); title('crit switch score v. flexibility');
                if reg_motion
                    figure;
                    suptitle(['Performance v. Flexibility, run ' num2str(run) ', MC']);
                    subplot(1,2,1);
                    correlate(flexToUse(~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                        'partial',motToUse(~isnan(dpToUse)),'type',corr_type,'ExtFigureCmd');
                    %correlate(sAges(~isnan(sDPrime)),stats.r,'ExtFigureCmd');
                    xlabel('age'); ylabel('dprime'); title('dprime v. age');
                    subplot(1,2,2);
                    correlate(flexToUse(~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                        'partial',motToUse(~isnan(csToUse)),'type',corr_type,'ExtFigureCmd');
                    %correlate(sAges(~isnan(sCSwitch)),stats2.r,'ExtFigureCmd');
                    xlabel('age'); ylabel('crit switch'); title('criterion switch score v. age');
                end    
              end
            end
        end
            
        if (ismember('corr_perf_ncomm',p_plot) && k==0  ) || ismember('save_stats',p_plot)
            % ncomms
            figure;
            suptitle('Performance v. Number of Communities');
            subplot(1,2,1);
            %ncToUse = mean(ncommS(1:3,:),1);
            ncToUse = ncommS(end,:);
            dpToUse = sDPrime(ib);
            csToUse = sCSwitch(ib);
            [rpsDPAll(12,1),rpsDPAll(12,2)] = correlate(ncToUse,dpToUse,'type','Pearson','ExtFigureCmd');
            xlabel('no. comms'); ylabel('dprime'); title('dprime v. no. comms');
            subplot(1,2,2);
            [rpsCSAll(12,1),rpsCSAll(12,2)] = correlate(ncToUse,csToUse,'type','Pearson','ExtFigureCmd');
            xlabel('no. comms'); ylabel('crit switch'); title('crit switch score v. no. comms');
            if reg_motion
                figure;
                suptitle('Performance v. ncomms, motion regressed');
                subplot(1,2,1);
                [rpsDPAll(12,3),rpsDPAll(12,4)] = correlate(ncToUse,dpToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd');
                xlabel('age'); ylabel('dprime'); title('dprime v. nc, MC');
                subplot(1,2,2);
                [rpsCSAll(12,3),rpsCSAll(12,4)] = correlate(ncToUse,csToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd');
                xlabel('age'); ylabel('crit shift'); title('criterion shift score v. nc, MC');
            end
            
            if is_runwise && t~=nruns
              for run=1:nruns
                ncToUse = RncommS(run,:);
                dpToUse = rDPrime(ib,run);
                csToUse = rCSwitch(ib,run);  
                figure;
                suptitle(['Performance v. Number of Communities, run ' num2str(run)]);
                subplot(1,2,1);
                correlate(ncToUse(~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),'ExtFigureCmd');
                xlabel('no. comms'); ylabel('dprime'); title('dprime v. no. comms');
                subplot(1,2,2);
                correlate(ncToUse(~isnan(csToUse)),csToUse(~isnan(csToUse)),'ExtFigureCmd');
                xlabel('no. comms'); ylabel('crit switch'); title('crit switch score v. no. comms');
              end
            end
        end
        
        if ismember('cohs',p_plot)
          for c=1:nC
            if ismember(legClass(c),c_plot) && k==0
            % plot subject-wise cohesions & null cohesions, sorted by average
            figure; xspace = 5;
            % sort subjects by mean coherence over tasks
            [sortedcohr,six] = sortentry(subjmcohr(:,:,c),'col',0,mean(subjmcohr(:,:,c),2));
            sortedavg = sortentry(avgmcohr(:,:,c),'col',0,mean(subjmcohr(:,:,c),2));
            plot(sortedcohr);
            hold on;
            plot(sortedavg,':');
            title(['cohesion of ' legClass(c) ' regions during each slice']);
            ylabel(['cohesion of ' legClass(c) ' regions']);
            xlabel('subject'); set(gca,'XTick',1:xspace:numel(ib),'XTickLabel',ib(six(1:xspace:end)));
            legend(numvec2cell(1:t),'Location','Northwest');
            
            % plot subject-wise p-values for cohesion, sorted by avg cohesion (same as above)
            figure;
            toplot = sortentry(mcohr_ps(:,:,2,c),'col',0,mean(subjmcohr(:,:,c),2));
            semilogy(toplot,'o'); hold on; semilogy(1:1:nsubs,0.05*ones(nsubs,1),'k:');
            legend(numvec2cell(1:t),'Location','NorthEast');
            title(['p-values for significant cohesion of ' legClass(c) 'regions']);
            xlabel('subject'); ylabel('p-value');
            set(gca,'XTick',1:xspace:numel(ib),'XTickLabel',ib(six(1:xspace:end)));
            end
          end
          
          % plot cohesion of each class, bare and then with class size
          % regressed out
          if k==0
              figure;
              plot(squeeze(mean(subjmcohr,2))',':');
              hold on;
              plot(squeeze(mean(mean(subjmcohr,2),1)),'-k');
              title('cohesion in each node class');
              xlabel('node class'); set(gca,'XTickLabel',legClass);
              ylabel('cohesion, avg. over runs');
              
              szs = histc(regCN(~~nodeix),(1:1:nC));
              stats = regstats(squeeze(mean(mean(subjmcohr,2),1)),szs);
              figure;
              plot(stats.r);
              title('node class cohesion, with effect of class size removed');
              xlabel('node class'); set(gca,'XTickLabel',legClass);
          end

        end
        
        if ((ismember('corr_age_cohsrecru',p_plot) || ismember('figure7',p_plot)) && k==0) || ismember('save_stats',p_plot)
            [rpsAgeAll(23,1),rpsAgeAll(23,2)] = correlate(sAges(ib),mean(mean(recruS(~~nodeix,:,:),1),3),'type','Spearman');
            title('age v. recruitment overall')
            [rpsAgeAll(23,3),rpsAgeAll(23,4)] = correlate(sAges(ib),mean(mean(recruS(~~nodeix,:,:),1),3),'partial',sMotion(ib),'DispResiduals','Raw','type','Spearman');
            title('age v. recruitment overall, MC')
            crs = zeros(nC,2); 
            cps = zeros(nC,2);
        for c=1:nC
            if ~ismember('figure7',p_plot) || ~ismember('save_stats',p_plot)
            if are_cohs
            % cohesion
            correlate(sAges(ib),squeeze(mean(subjmcohr(:,:,c),2)),'type',corr_type);
            title(legClass(c))
%%%%%%%%%%%%%%%%%%%%%%% regress motion
            if reg_motion
                correlate(sAges(ib),squeeze(mean(subjmcohr(:,:,c),2)),...
                    'partial',sMotion(ib),'type',corr_type,'DispResiduals','Raw');
                % stats = regstats(squeeze(mean(subjmcohr(:,:,c),2)),sMotion(ib));
                % correlate(sAges(ib),stats.r)
                xlabel('ages'); ylabel('class cohesion');
                title([legClass{c}, ' MC']);
            end
            end
            end
            
            % recruitment
            [crs(c,1),cps(c,1)] = correlate(sAges(ib),squeeze(mean(clRecru(c,:,:),3)),'type',corr_type);
            xlabel('ages'); ylabel('class recruitment');
            title(legClass(c))
            if reg_motion   %%%%% regress motion
                [crs(c,2),cps(c,2)] = correlate(sAges(ib),mean(clRecru(c,:,:),3),...
                    'partial',sMotion(ib),'type',corr_type,'DispResiduals','Raw');
                % stats = regstats(squeeze(mean(subjmcohr(:,:,c),2)),sMotion(ib));
                % correlate(sAges(ib),stats.r)
                xlabel('ages'); ylabel('class recruitment');
                title([legClass{c}, ' MC']);
            end
            % play around with multiple comparisons tests
%             [fdrPC,fdrNA] = FDR(cps(:,1),0.05);
%             [fdrPC2,fdrNA2] = FDR(cps(:,2),0.05);
%             fdrs = repmat([fdrPC,fdrPC2],nC,1);
%             fdrsNA = repmat([fdrNA,fdrNA2],nC,1);
%             legClass(cps(:,1)<fdrs(:,1))
%             legClass(cps(:,2)<fdrs(:,2))
%             legClass(cps(:,1)<fdrsNA(:,1))
%             legClass(cps(:,2)<fdrsNA(:,2))
            
            if ~ismember('figure7',p_plot) || ~ismember('save_stats',p_plot)
            if is_runwise && t~=nruns
              for run=1:nruns
                % recruitment
                correlate(sAges(ib),clRecru(c,:,run),'type',corr_type);
                xlabel('ages'); ylabel('class recruitment');
                title([legClass{c} ', run ' num2str(run)])
                if reg_motion %%%%%%%% regress motion
                    correlate(sAges(ib),clRecru(c,:,run),...
                        'partial',sMotion(ib),'type',corr_type);
                    % stats = regstats(squeeze(mean(subjmcohr(:,:,c),2)),sMotion(ib));
                    % correlate(sAges(ib),stats.r)
                    xlabel('ages'); ylabel('class recruitment');
                    title([legClass{c} ', run ' num2str(run) ', MC']);
                end  
              end
            end
            end
            
        end
            rpsAgeAll(24:33,1) = crs(:,1);
            rpsAgeAll(24:33,2) = cps(:,1);
            rpsAgeAll(24:33,3) = crs(:,2);
            rpsAgeAll(24:33,4) = cps(:,2);
        end
        
        % are motion affected classes smaller?
        rest = zeros(nC,1); pest=rest;
        rest1 = zeros(nC,1); pest1=rest;
        rest2 = zeros(nC,1); pest2=rest;
        for c=1:nC; 
            toCorrrr = clFlex(c,:);
            toCorrr = mean(clRecru(c,:,:),3);
            toCorr = cncommS(c,:,end);
            [rest(c),pest(c)]=correlate(sMotion(ib),toCorrrr,'NoFigure'); 
            [rest1(c),pest1(c)]=correlate(sMotion(ib),toCorrr,'NoFigure'); 
            [rest2(c),pest2(c)]=correlate(sMotion(ib),toCorr,'NoFigure'); 
        end
        figure;
        [rrr,ppp] = correlate(clSize,rest,'ExtFigureCmd'); hold all;
        ylabel('correlation: motion v. system-specific measure')
        xlabel('system size')
        [rrr1,ppp1] = correlate(clSize,rest1,'ExtFigureCmd');
        [rrr2,ppp2] = correlate(clSize,rest2,'ExtFigureCmd');
        legend(['flexibility, r = ' num2str(rrr) ', N.S.'],...
               ['recruitment, r = ' num2str(rrr1) ', p < 0.05'],...
               ['n comms, r = ' num2str(rrr2) ', N.S.']);
        
        
        
        if (ismember('corr_perf_cohsrecru',p_plot) && k==0) || ismember('save_stats',p_plot)
            figure;
              recoToUse = mean(mean(recruS(~~nodeix,:,:),1),3);
             % overall recruitment
              subplot(1,2,1);
              [rpsDPAll(23,1),rpsDPAll(23,2)] = correlate(recoToUse,dpToUse,'type','Pearson','ExtFigureCmd',fig_tag);
              xlabel('recruitment'); ylabel('dprime'); title('dprime v. recruitment');
              subplot(1,2,2);
              [rpsCSAll(23,1),rpsCSAll(23,2)] = correlate(recoToUse,csToUse,'type','Pearson','ExtFigureCmd',fig_tag);
              xlabel('recruitment'); ylabel('crit shift'); title('crit shift score v. recruitment');
              if reg_motion
                  figure;
                  suptitle('Performance v. ncomms, motion regressed');
                  subplot(1,2,1);
                  [rpsDPAll(23,3),rpsDPAll(23,4)] = correlate(recoToUse,dpToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd',fig_tag);
                  xlabel('age'); ylabel('dprime'); title('dprime v. recru, MC');
                  subplot(1,2,2);
                  [rpsCSAll(23,3),rpsCSAll(23,4)] = correlate(recoToUse,csToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd',fig_tag);
                  xlabel('age'); ylabel('crit shift'); title('criterion shift score v. recru, MC');
              end
            
            rpsDP = zeros(nC,4);  
            rpsCS = zeros(nC,4);
            for c=1:nC
                if are_cohs
                % cohesion
                dpToUse = sDPrime(ib);
                csToUse = sCSwitch(ib);
                cohsToUse = squeeze(mean(subjmcohr(:,:,c),2));
                figure;
                suptitle(legClass(c));
                subplot(1,2,1);
                correlate(cohsToUse(~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                    'type',corr_type,'ExtFigureCmd');
                xlabel('cohesion'); ylabel('dprime'); title('dprime v. cohesion');
                subplot(1,2,2);
                correlate(cohsToUse(~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                    'type',corr_type,'ExtFigureCmd');
                xlabel('cohesion'); ylabel('crit shift'); title('crit shift score v. cohesion');
                end
                
                % recruitment
                recToUse = mean(clRecru(c,:,:),3);
                figure;
                suptitle(legClass(c));
                subplot(1,2,1);
                [rpsDP(c,1),rpsDP(c,2)] = correlate(recToUse(~isnan(dpToUse)),dpToUse,'type','Pearson','ExtFigureCmd');
                xlabel('recruitment'); ylabel('dprime'); title('dprime v. recruitment');
                subplot(1,2,2);
                [rpsCS(c,1),rpsCS(c,2)] = correlate(recToUse,csToUse,'type','Pearson','ExtFigureCmd');
                xlabel('recruitment'); ylabel('crit shift'); title('crit shift score v. recruitment');
                if reg_motion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% crecru perfs typo in rpsCS
                    figure;
                    suptitle('Performance v. ncomms, motion regressed');
                    subplot(1,2,1);
                    [rpsDP(c,3),rpsDP(c,4)] = correlate(recToUse,dpToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd');
                    xlabel('age'); ylabel('dprime'); title('dprime v. recru, MC');
                    subplot(1,2,2);
                    [rpsCS(c,3),rpsCS(c,4)] = correlate(recToUse,csToUse,'partial',motToUse,'type','Pearson','ExtFigureCmd');
                    xlabel('age'); ylabel('crit shift'); title('criterion shift score v. recru, MC');
                end
                
                
                if is_runwise && t~=nruns
                  for run=1:nruns
                    dpToUse = rDPrime(ib,run);
                    csToUse = rCSwitch(ib,run);
                    recToUse = clRecru(c,:,run);
                    figure;
                    suptitle([legClass{c} ', run ' num2str(run)]);
                    subplot(1,2,1);
                    correlate(recToUse(~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                        'type',corr_type,'ExtFigureCmd');
                    xlabel('recruitment'); ylabel('dprime'); title('dprime v. recruitment');
                    subplot(1,2,2);
                    correlate(recToUse(~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                        'type',corr_type,'ExtFigureCmd');
                    xlabel('recruitment'); ylabel('crit shift'); title('crit shift score v. recruitment');
                  end
                end
            end 
            rpsDPAll(24:33,:) = rpsDP;
            rpsCSAll(24:33,:) = rpsCS;
                
        end
        
        
        % analysis of coherences in relation to null model
        if ismember('perfcohs',p_plot)
          for c=1:nC
            if ismember(legClass(c),c_plot) && k==0
            sigs = (mcohr_ps(:,:,2,c)<0.05);
            allx = all(sigs,2);
            somex = any(sigs,2) - all(sigs,2);
            nonex = 1 - any(sigs,2);
            dpToUse = sDPrime(ib);
            csToUse = sCSwitch(ib);
            %figure;
            %scatter(dp)
            
            %correlate(sum(sigs(~isnan(dpToUse),:),2),dpToUse(~isnan(dpToUse)));
            %legend('');
            %correlate(sum(sigs(~isnan(csToUse),:),2),csToUse(~isnan(csToUse)));
            %legend('');
            
            semilogy(toplot,'o'); hold on; semilogy(1:1:nsubs,0.05*ones(nsubs,1),'k:');
            legend(numvec2cell(1:t),'Location','NorthEast');
            title(['p-values for significant cohesion of ' legClass(c) 'regions']);
            xlabel('subject'); ylabel('p-value');
            set(gca,'XTick',1:xspace:numel(ib),'XTickLabel',ib(six(1:xspace:end)));
            end
          end
          
        end
        
        if ismember('recru',p_plot) || ismember('figure6c',p_plot)
          if k~=0 && ~ismember('figure6c',p_plot)
            % plot recruitment for all those people.
            BrainView(posPass,mean(recruS(~~nodeix,k,:),3),0,'View','sagittal');
            title(['recruitment coefficient over whole experiment, subject ' subj_IDs{k}]);
            colorbar;
          else
            BrainView(posPass,mean(mean(recruS(~~nodeix,:,:),3),2),0,'View','sagittal');
            title('recruitment coefficient over whole experiment, all subjects');
            colorbar;
          end
        end
        
        if ismember('crecru',p_plot) || ismember('figure6d',p_plot)
          if k==0 %&& ismember(legClass(c),c_plot)
            % for each class, plot mean recruitment for nodes in that class   
            % for each subject (as well as the overall subject mean)
%             figure;
%             plot(mean(clRecru,3),':'); hold on;
%             plot(diag(mean(mean(corecruS,4),3)),'-k');
%             xlabel('node class'); set(gca,'XTickLabel',legClass);
%             ylabel('mean recruitment coefficient');
%             title('node recruitment by class');
            figure;
            boxplot(mean(clRecru,3)'); hold on;
            plot(diag(mean(mean(corecruS,4),3)),'ok');
            xlabel('node system'); set(gca,'XTickLabel',legClass);
            ylabel('system recruitment coefficient')
            title('node recruitment by system')
            
            if ~ismember('figure6d',p_plot)
            % now plot for each run as well as the overall mean
            figure;
            plot(squeeze(mean(clRecru,2)),':'); hold on;
            plot(diag(mean(mean(corecruS,3),4)),'-k');
            xlabel('node class'); set(gca,'XTickLabel',legClass);
            ylabel('mean recruitment coefficient');
            title('node recruitment by class');
            legend(numvec2cell(1:nruns));
            end
          end
        end
           
        if ismember('ccR',p_plot)
            if k==0
                % co-recruitments
                figure; bcolor(mean(mean(corecruS,3),4));
                xlabel('node class'); set(gca,'XTickLabel',legClass);
                ylabel('node class'); set(gca,'YTickLabel',legClass);
                title('class co-recruitments, averaged over all subjects');
                colorbar;
            else
                % individual co-recruitments
                figure; bcolor(mean(corecruS(:,:,k,:),4));
                xlabel('node class'); set(gca,'XTickLabel',legClass);
                ylabel('node class'); set(gca,'YTickLabel',legClass);
                title(['average class co-recruitments over slices, subject ' subj_IDs{k}]);
                colorbar;
            end
        end


        if ismember('cflex',p_plot) || ismember('figure6b',p_plot)
            % plot flexibility of each class
%             figure;
%             plot(clFlex,':'); hold on;
%             plot(mean(clFlex,2),'-k');
%             xlabel('node class'); set(gca,'XTickLabel',legClass);
%             ylabel('mean flexibility over nodes in class')
%             title('flexibility of individual node classes')
            
            figure;
            boxplot(clFlex'); hold on;
            plot(mean(clFlex,2),'ok');
            xlabel('node class'); set(gca,'XTickLabel',legClass);
            ylabel('mean flexibility over nodes in class')
            title('flexibility of individual node classes')
            
            if ~ismember('figure6b',p_plot)
            % check sig differences between classes
            [~,ix] = sort(mean(clFlex,2));
            disp(mean(clFlex(ix,:),2));
            ttest_flex = zeros(nCl);
            for i=1:nC 
              for j=1:nC 
                [~,ttest_flex(i,j)] = ttest(clFlex(i,:),clFlex(j,:)); 
              end
            end
            figure;
            bcolor(ttest_flex); colorbar;
            figure;
            bcolor(ttest_flex<0.05/(nC*2));colorbar;
            figure;
            bcolor(ttest_flex<0.05);colorbar;
            %sum(ttest_flex<0.05/(nC*2));
            end
        end
        
        if (ismember('corr_age_cflex',p_plot) && k==0) || ismember('save_stats',p_plot)
            % correlation of age with class flexibility
            rps = zeros(nC,4);
            for c=1:nC
            [rps(c,1),rps(c,2)] = correlate(sAges(ib),clFlex(c,:),'type','Spearman',fig_tag);
            title(legClass(c))
            xlabel('ages'); ylabel('class flexibility');
            if reg_motion  %%%%% regress motion
                [rps(c,3),rps(c,4)] = correlate(sAges(ib),clFlex(c,:),...
                    'partial',sMotion(ib),'type','Spearman','DispResiduals','Raw',fig_tag);
                xlabel('ages'); ylabel('class flexibility');
                title([legClass{c}, ' motion regressed']);
            end
            
            if is_runwise && t~=nruns
              for run=1:nruns
                correlate(sAges(ib),RclFlex(c,run,:),'type',corr_type);
                title([legClass{c} ', run ' num2str(run)]);
                xlabel('ages'); ylabel('class flexibility');
                if reg_motion  %%%%% regress motion
                    correlate(sAges(ib),RclFlex(c,run,:),...
                        'partial',sMotion(ib),'type',corr_type);
                    xlabel('ages'); ylabel('class flexibility');
                    title([legClass{c} ', run ' num2str(run) ', MC']);
                end  
              end
            end
            
            end
            rpsAgeAll(2:11,:) = rps;
        end

        if (ismember('corr_perf_cflex',p_plot) && k==0) || ismember('save_stats',p_plot)
            % correlation of performance with class flexibility
            rpsDP = zeros(nC,4);
            rpsCS = zeros(nC,4);
            for c=1:nC
                dpToUse = sDPrime(ib);
                csToUse = sCSwitch(ib);
                figure;
                suptitle(legClass(c));
                subplot(1,2,1);
                [rpsDP(c,1),rpsDP(c,2)] = correlate(clFlex(c,~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                    'type',corr_type,'ExtFigureCmd',fig_tag);
                xlabel('class flexibility'); ylabel('dprime'); title('dprime v. class flexibility');
                subplot(1,2,2);
                [rpsCS(c,1),rpsCS(c,2)] = correlate(clFlex(c,~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                    'type',corr_type,'ExtFigureCmd',fig_tag);
                xlabel('class flexibility'); ylabel('crit shift'); title('crit shift score v. class flexibility');
            if reg_motion  %%%%%% regress motion
                motToUse = sMotion(ib);
                figure;
                suptitle([legClass{c}, ' motion regressed']);
                subplot(1,2,1);
                [rpsDP(c,3),rpsDP(c,4)] = correlate(clFlex(c,~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                    'partial',motToUse(~isnan(dpToUse)),'type',corr_type,'ExtFigureCmd',fig_tag);
                xlabel('class flexibility'); ylabel('dprime'); title('dprime v. class flexibility');
                subplot(1,2,2);
                [rpsCS(c,3),rpsCS(c,4)] = correlate(clFlex(c,~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                    'partial',motToUse(~isnan(csToUse)),'type',corr_type,'ExtFigureCmd',fig_tag);
                xlabel('class flexibility'); ylabel('crit shift'); title('crit shift score v. class flexibility');
            end
            
            if is_runwise && t~=nruns
              for run=1:nruns
                dpToUse = rDPrime(ib,run);
                csToUse = rCSwitch(ib,run);
                figure;
                suptitle([legClass{c} ', run ' num2str(run)]);
                subplot(1,2,1);
                correlate(RclFlex(c,run,~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                    'type','Pearson','ExtFigureCmd');
                xlabel('class flexibility'); ylabel('dprime'); title('dprime v. class flexibility');
                subplot(1,2,2);
                correlate(RclFlex(c,run,~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                    'type','Pearson','ExtFigureCmd');
                xlabel('class flexibility'); ylabel('crit shift'); title('crit shift score v. class flexibility');
                if reg_motion  %%%%%% regress motion
                    figure;
                    suptitle([legClass{c} ',run ' num2str(run) ', MC']);
                    subplot(1,2,1);
                    correlate(RclFlex(c,run,~isnan(dpToUse)),dpToUse(~isnan(dpToUse)),...
                        'partial',motToUse(~isnan(dpToUse)),'type','Pearson','ExtFigureCmd');
                    xlabel('class flexibility'); ylabel('dprime'); title('dprime v. class flexibility');
                    subplot(1,2,2);
                    correlate(RclFlex(c,run,~isnan(csToUse)),csToUse(~isnan(csToUse)),...
                        'partial',motToUse(~isnan(csToUse)),'type','Pearson','ExtFigureCmd');
                    xlabel('class flexibility'); ylabel('crit shift'); title('crit shift score v. class flexibility');
                end  
              end
            end
            
            end
            rpsDPAll(2:11,:) = rpsDP;
            rpsCSAll(2:11,:) = rpsCS;
        end
        
        if ismember('adjs',p_plot) && k~=0
            figure; 
            for T=1:t
            subplot(1,t,T);
            bcolor(A(k).adj{T});
            title(['slice ' num2str(T)])
            end
            %colorbar('Location','SouthOutside');
            suptitle(['adjacency matrices, subject ' subj_IDs{k}])
        end
        
        if ismember('flex_noise_diag',p_plot) && k~=0
            % mean flex v. overall node degree, PER NODE
            correlate(mean(flexS(~~nodeix,:),2),mean(mean(overallND,2),3));
            title('overall node degree v flex, mean over subjs')
            
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(overallND,2),1));
            title('overall node degree v flex, mean over nodes')  
            
            % mean flex v. overall node degree, PER NODE
            correlate(mean(flexS(~~nodeix,:),2),mean(mean(inND,2),3));
            title('in-comm node degree v flex, mean over subjs')
            
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(inND,2),1));
            title('in-comm node degree v flex, mean over nodes') 
            
            % mean flex v. overall node degree, PER NODE
            correlate(mean(flexS(~~nodeix,:),2),mean(mean(outND,2),3));
            title('out-comm node degree v flex, mean over subjs')
            
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(outND,2),1));
            title('out-comm node degree v flex, mean over nodes') 
            
            % ratio of in to out connections
            inDout = mean(mean(inND,2),3)./mean(mean(outND,2),3);
            correlate(mean(flexS(~~nodeix,:),2),inDout);
            title('in/out connection ratio v flex, mean over subjs')
            
        end
        
        if ismember('flex_noise_diag_norm',p_plot) && k~=0
            % mean flex v. overall node degree, PER NODE
            correlate(mean(flexS(~~nodeix,:),2),mean(mean(overallND,2),3));
            title('overall node degree v flex, mean over subjs')
            
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(overallND,2),1));
            title('overall node degree v flex, mean over nodes')  
            
            % mean flex v. overall node degree, PER NODE
            correlate(mean(flexS(~~nodeix,:),2),mean(mean(inND./cSize,2),3));
            title('in-comm node degree v flex, mean over subjs')
            
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(inND./cSize,2),1));
            title('in-comm node degree v flex, mean over nodes') 
            
            % mean flex v. overall node degree, PER NODE
            correlate(mean(flexS(~~nodeix,:),2),mean(mean(outND./(n-cSize),2),3));
            title('out-comm node degree v flex, mean over subjs')
            
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(outND./(n-cSize),2),1));
            title('out-comm node degree v flex, mean over nodes') 
            
            % 3D plot
            figure;
            scatter3(mean(flexS(~~nodeix,:),2),mean(mean(inND./cSize,2),3),mean(mean(overallND,2),3));
            xlabel('flexibility'); ylabel('in-degree'); zlabel('overall degree');
        end
        
        if ismember('flex_noise_diag_comm',p_plot) && k~=0
            
          for c=1:nC
            if c==1; figure; end;
            % mean flex v. overall node degree, PER NODE
            inC = (regCN(~~nodeix)==c);
            inCx = (regCN==c).*(~~nodeix);
            correlate(mean(flexS(~~inCx,:),2),...
                mean(mean(overallND(~~inC,:,:),2),3),'ExtFigureCmd');
            hold all;
            title('overall node degree v flex, mean over subjs')
          end
          legend(legClass);

            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(overallND,2),1));
            title('overall node degree v flex, mean over nodes')  
            
         for c=1:nC  
            if c==1; figure; end;
            % mean flex v. overall node degree, PER NODE
            inC = (regCN(~~nodeix)==c);
            inCx = (regCN==c).*(~~nodeix);
            correlate(mean(flexS(~~inCx,:),2),...
                mean(mean(inND(~~inC,:,:),2),3),'ExtFigureCmd');
            hold all;
            title('in-comm node degree v flex, mean over subjs')
         end
         legend(legClass);

            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(inND,2),1));
            title('in-comm node degree v flex, mean over nodes') 
            
         for c=1:nC   
            % mean flex v. overall node degree, PER NODE
            inC = (regCN(~~nodeix)==c);
            inCx = (regCN==c).*(~~nodeix);
            if c==1; figure; end;
            correlate(mean(flexS(~~inCx,:),2),...
                mean(mean(outND(~~inC,:,:),2),3),'ExtFigureCmd');
            hold all;
            title('out-comm node degree v flex, mean over subjs')
         end
         legend(legClass);
         
            % mean flex v. overall node degree, PER SUBJ
            correlate(mean(flexS(~~nodeix,:),1),mean(mean(outND,2),1));
            title('out-comm node degree v flex, mean over nodes')
            
          for c=1:nC   
            % ratio of in to out connections
            if c==1; figure; end;
            inC = (regCN(~~nodeix)==c);
            inCx = (regCN==c).*(~~nodeix);
            inDout = mean(mean(inND,2),3)./mean(mean(outND,2),3);
            correlate(mean(flexS(~~inCx,:),2),...
                inDout(~~inC,:,:),'ExtFigureCmd');
            hold all;
            title('in/out connection ratio v flex, mean over subjs')
          end  
          legend(legClass);
        end
        
        if ismember('cncomm',p_plot) && k==0
            figure;
            plot(cncommS(:,:,end),':'); hold on;
            plot(mean(cncommS(:,:,end),2),'-k')
            xlabel('node class'); set(gca,'XTickLabel',legClass);
            ylabel('number of communities in node class');
            title('number of communities by node class');
        end
        
        
        if (ismember('corr_age_cncomm',p_plot) && k==0) || ismember('save_stats',p_plot)
            % correlation of age with class ncomms
            rps = zeros(nC,4);
            for c=1:nC
            [rps(c,1),rps(c,2)] = correlate(sAges(ib),cncommS(c,:,end),'type','Spearman',fig_tag);
            title(legClass(c))
            xlabel('ages'); ylabel('class ncomms');
            if reg_motion %%%%%% regress motion
                [rps(c,3),rps(c,4)] = correlate(sAges(ib),cncommS(c,:,end),...
                    'partial',sMotion(ib),'type','Spearman','DispResiduals','Raw',fig_tag);
                xlabel('ages'); ylabel('class ncomms');
                title([legClass{c}, ' motion regressed']);
            end
            end
            
            rpsAgeAll(13:22,:) = rps;            
%             if is_runwise && t~=nruns
%               for run=1:nruns
%                 correlate(sAges(ib),RcncommS(c,:,end),'type',corr_type);
%                 title([legClass{c} ',run ' num2str(run)])
%                 xlabel('ages'); ylabel('class ncomms');
%                 if reg_motion %%%%%% regress motion
%                     correlate(sAges(ib),RcncommS(c,:,end),...
%                         'partial',sMotion(ib),'type',corr_type);
%                     xlabel('ages'); ylabel('class ncomms');
%                     title([legClass{c} ',run ' num2str(run) ', MC']);
%                 end  
%               end
%             end
            
        end
        
        if (ismember('corr_perf_cncomm',p_plot) && k==0) || ismember('save_stats',p_plot)
            rpsDP = zeros(nC,4);
            rpsCS = zeros(nC,4);
            % correlation of performance with class ncomms
            dpToUse = sDPrime(ib);
            csToUse = sCSwitch(ib);
            for c=1:nC
                figure;
                suptitle(legClass(c));
                subplot(1,2,1);
                [rpsDP(c,1),rpsDP(c,2)] = correlate(cncommS(c,~isnan(dpToUse),end),dpToUse(~isnan(dpToUse)),...
                    'type','Pearson','ExtFigureCmd',fig_tag);
                xlabel('class ncomms'); ylabel('dprime'); title('dprime v. class ncomms');
                subplot(1,2,2);
                [rpsCS(c,1),rpsCS(c,2)] = correlate(cncommS(c,~isnan(csToUse),end),csToUse(~isnan(csToUse)),...
                    'type','Pearson','ExtFigureCmd',fig_tag);
                xlabel('class ncomms'); ylabel('crit shift'); title('crit shift score v. class ncomms');
            if reg_motion  %%%%%% regress motion
                motToUse = sMotion(ib);
                figure;
                suptitle([legClass{c}, ' motion regressed']);
                subplot(1,2,1);
                [rpsDP(c,3),rpsDP(c,4)] = correlate(cncommS(c,~isnan(dpToUse),end),dpToUse(~isnan(dpToUse)),...
                    'partial',sMotion(ib),'type','Pearson','ExtFigureCmd',fig_tag);
                xlabel('class ncomms'); ylabel('dprime'); title('dprime v. class ncomms');
                subplot(1,2,2);
                [rpsCS(c,3),rpsCS(c,4)] = correlate(cncommS(c,~isnan(csToUse),end),csToUse(~isnan(csToUse)),...
                    'partial',sMotion(ib),'type','Pearson','ExtFigureCmd',fig_tag);
                xlabel('class ncomms'); ylabel('crit shift'); title('crit shift score v. class ncomms');
            end
            end
            rpsDPAll(13:22,:) = rpsDP;
            rpsCSAll(13:22,:) = rpsCS;
        end
        
      if ismember('figure2',p_plot)
        %%%% performance correlation
        % using only N=104... subjects with behavioral info and brain scans
        ageToUse = sAges(ib); 
        dpToUse = sDPrime(ib); 
        csToUse = sCSwitch(ib);        
        figure;
        scatter(dpToUse(ageToUse>=60),csToUse(ageToUse>=60),'ob')
        hold on;
        scatter(dpToUse(ageToUse==18),csToUse(ageToUse==18),'ok')
        scatter(dpToUse(~~((ageToUse>18).*(ageToUse<60))),csToUse(~~((ageToUse>18).*(ageToUse<60))),'or')
        title('Word Memory Performance Measures')
        xlabel('dprime'); ylabel('criterion shift score')
        legend('age 18','ages 25-33','ages 65-70')
        % should give Pearson's r = -0.06, p > 0.5
        correlate(dpToUse,csToUse,'type','Pearson');
      end
        
        % community sizes
      if ismember('comm_sizes',p_plot) && k==0
        figure; subplot(2,2,1); bcolor(squeeze(cMaxMin(1:t,1,:))'); hold on;
        colorbar('Location','EastOutside');
        subplot(2,2,2); bcolor(squeeze(cMaxMin(1:t,2,:))'); hold on;
        colorbar('Location','EastOutside');
        subplot(2,2,3); bcolor(squeeze(cMaxMin(end,1,:))); hold on;
        colorbar('Location','EastOutside');
        subplot(2,2,4); bcolor(squeeze(cMaxMin(end,2,:))); hold on;
        colorbar('Location','EastOutside');
        
        % histograms
        % all & slicewise
        cSizeHistA = [];
        cSizeHistT = [];
        for j=1:nsubs
            cSizeHistA = [cSizeHistA; cSizesA(1:ncommS(end,j),j)];
            for T=1:t
                cSizeHistT = [cSizeHistT; cSizesT(1:ncommS(end,j),T,j)];
            end
        end
        % figures
        figure; hist(cSizeHistA,50);
        figure; hist(cSizeHistT,50);
      end
      
      if ismember('removevis',p_plot)
        if k==0
             for i=randi(nsubs,[10,1])'; 
                 figure; bcolor(partnS(~~nodeix,:,i)); colorbar; 
             end
            flex_rv = zeros(n,nsubs);
            for i=1:nsubs
                flex_rv(:,i) = add_missings_dim(...
                    flexibility(partnS(~~nodeix,:,i)',1),find(~~nodeix),n,1,0 ...
                                                ); % categorical
            end
            figure; plot(mean(flex_rv(~~nodeix,:),1));
            [sortedflex,sortix] = sort(mean(flex_rv(~~nodeix,:),1));
            figure; plot(sortedflex);
            %hist(mean(flex_rv(~~nodeix,:),1),nsubs);                     
            for i=[sortix(1:10),sortix(end-9:end)]; 
                figure; bcolor(partnS(~~nodeix,:,i)); colorbar; 
            end
        end
      end
      
      if ismember('surety',p_plot)
        if k==0
           % correlate surety with flexibility
             f1 = squeeze(mean(sureS(~~nodeix,:,:,1),2));
             f2 = flexS(~~nodeix,:);
             rs = zeros(nsubs,1);
             ps = zeros(nsubs,1);
           figure;
             for i=1:nsubs
                 [rs(i),ps(i)] = correlate(f1(:,i),f2(:,i),'ExtFigureCmd');
                 hold on;
             end
             hold off;
           figure;
             [rall,pall] = correlate(f1(:),f2(:),'ExtFigureCmd');    
           figure; hist(rs,100);
           figure; hist(ps,100);
           figure; hist(rs(ps<0.05/nsubs),sum(ps<0.05/nsubs));
           [rt,pt]=ttest(rs);
           % there are huge individual differences in flexibility
           % distributions, but relatively little among surety
           % distributions.
           figure;
             suptitle('flexibility v. surety: significant correlations');
             j=1;
             for i=find(ps<0.05/nsubs)'
                 subplot(4,5,j);
                 correlate(f1(:,i),f2(:,i),'ExtFigureCmd','RemoveLeg');
                 title(['subject ' num2str(ib(i))]);
                 j = j + 1;
             end
             %xlabel('surety');
             %ylabel('flexibility');
        end
      end
      
      if ismember('null_tests',p_plot) %&& ts==316
          rage = zeros(2+2*nC,nullp,3);
          page = zeros(2+2*nC,nullp,3);
          for j=1:nullp
            % correlations with age:
              % flexibility overall
              [rage(1,j,1),page(1,j,1)] = correlate(sAges(ib),mean(flexRand(~~nodeix,:,j),1),...
                  'type','Spearman','NoFigure');
              % recruitment overall
              [rage(2,j,1),page(2,j,1)] = correlate(sAges(ib),mean(mean(recruRand(~~nodeix,:,:,j),3),1),...
                  'type','Spearman','NoFigure');
            for cl=1:nC  
              % recruitment classwise
              [rage(2+cl,j,1),page(2+cl,j,1)] = correlate(sAges(ib),mean(clRecruRand(cl,:,:,j),3),...
                  'type','Spearman','NoFigure');
            end
              
            % correlations with dprime:
              % flexibility overall
              [rage(1,j,2),page(1,j,2)] = correlate(sAges(ib),mean(flexRand(~~nodeix,:,j),1),...
                  'type','Pearson','NoFigure');
              % recruitment overall
              [rage(2,j,2),page(2,j,2)] = correlate(sAges(ib),mean(mean(recruRand(~~nodeix,:,:,j),3),1),...
                  'type','Pearson','NoFigure');
            for cl=1:nC  
              % recruitment classwise
              [rage(2+cl,j,2),page(2+cl,j,2)] = correlate(sAges(ib),mean(clRecruRand(cl,:,:,j),3),...
                  'type','Pearson','NoFigure');
            end

            % correlations with cshift:
              % flexibility overall
              [rage(1,j,3),page(1,j,3)] = correlate(sAges(ib),mean(flexRand(~~nodeix,:,j),1),...
                  'type','Pearson','NoFigure');
              % recruitment overall
              [rage(2,j,3),page(2,j,3)] = correlate(sAges(ib),mean(mean(recruRand(~~nodeix,:,:,j),3),1),...
                  'type','Pearson','NoFigure');
            for cl=1:nC  
              % recruitment classwise
              [rage(2+cl,j,3),page(2+cl,j,3)] = correlate(sAges(ib),mean(clRecruRand(cl,:,:,j),3),...
                  'type','Pearson','NoFigure');
            end
          end
          
          % null checks and visualization
        figure;
          [rr,pp] = correlate(sAges(ib),mean(flexS(~~nodeix,:),1),...
                  'type','Spearman','ExtFigureCmd','RemoveLeg'); hold all;
          correlate(sAges(ib),mean(mean(flexRand(~~nodeix,:,:),1),3),...
                  'type','Spearman','ExtFigureCmd','RemoveLeg');
          title('age v. flexibility overall and null');
          legend(['p = ' num2str(sum(rage(1,:,1)>rr)./nullp)]); hold off;
          sum(page(1,:,1)<0.05);
        figure; hist(rage(1,:,1));
          legend(['r = ' num2str(rr) ', p = ' num2str(sum(rage(1,:,1)>rr)./nullp)]);
          
        figure;
          [rr,pp] = correlate(sAges(ib),mean(mean(recruS(~~nodeix,:,:),3),1),...
                  'type','Spearman','ExtFigureCmd','RemoveLeg'); hold all;
          correlate(sAges(ib),mean(mean(mean(recruRand(~~nodeix,:,:,:),3),1),4),...
                  'type','Spearman','ExtFigureCmd','RemoveLeg');
          title('age v. recruitment overall and null');
          legend(['p = ' num2str(sum(rage(2,:,1)<rr)./nullp)]); hold off;
        figure; hist(rage(2,:,1));
          legend(['r = ' num2str(rr) ', p = ' num2str(sum(rage(2,:,1)<rr)./nullp)]);

        for cl=1:nC
        figure;
          [rr,pp] = correlate(sAges(ib),mean(clRecru(cl,:,:),3),...
                  'type','Spearman','ExtFigureCmd','RemoveLeg'); hold all;
          correlate(sAges(ib),mean(mean(clRecruRand(cl,:,:,:),3),4),...
                  'type','Spearman','ExtFigureCmd','RemoveLeg');
          title(['age v. recruitment ' legClass{cl} ' and null']);
          legend(['p = ' num2str(sum(rage(2+cl,:,1)<rr)./nullp)]); hold off;
        figure; hist(rage(2+cl,:,1));
          legend(['r = ' num2str(rr) ', p = ' num2str(sum(rage(2+cl,:,1)<rr)./nullp)]);
        
        end
        
      end
      
      
      if ismember('null_tests',p_plot) && reg_motion %&& ts==316 
          rage = zeros(2+2*nC,nullp,3);
          page = zeros(2+2*nC,nullp,3);
          for j=1:nullp
            % correlations with age:
              % flexibility overall
              [rage(1,j,1),page(1,j,1)] = correlate(sAges(ib),mean(flexRand(~~nodeix,:,j),1),...
                  'type','Spearman','partial',sMotion(ib),'NoFigure');
              % recruitment overall
              [rage(2,j,1),page(2,j,1)] = correlate(sAges(ib),mean(mean(recruRand(~~nodeix,:,:,j),3),1),...
                  'type','Spearman','partial',sMotion(ib),'NoFigure');
            for cl=1:nC  
              % recruitment classwise
              [rage(2+cl,j,1),page(2+cl,j,1)] = correlate(sAges(ib),mean(clRecruRand(cl,:,:,j),3),...
                  'type','Spearman','partial',sMotion(ib),'NoFigure');
            end
              
            % correlations with dprime:
              % flexibility overall
              [rage(1,j,2),page(1,j,2)] = correlate(sAges(ib),mean(flexRand(~~nodeix,:,j),1),...
                  'type','Pearson','partial',sMotion(ib),'NoFigure');
              % recruitment overall
              [rage(2,j,2),page(2,j,2)] = correlate(sAges(ib),mean(mean(recruRand(~~nodeix,:,:,j),3),1),...
                  'type','Pearson','partial',sMotion(ib),'NoFigure');
            for cl=1:nC  
              % recruitment classwise
              [rage(2+cl,j,2),page(2+cl,j,2)] = correlate(sAges(ib),mean(clRecruRand(cl,:,:,j),3),...
                  'type','Pearson','partial',sMotion(ib),'NoFigure');
            end

            % correlations with cshift:
              % flexibility overall
              [rage(1,j,3),page(1,j,3)] = correlate(sAges(ib),mean(flexRand(~~nodeix,:,j),1),...
                  'type','Pearson','partial',sMotion(ib),'NoFigure');
              % recruitment overall
              [rage(2,j,3),page(2,j,3)] = correlate(sAges(ib),mean(mean(recruRand(~~nodeix,:,:,j),3),1),...
                  'type','Pearson','partial',sMotion(ib),'NoFigure');
            for cl=1:nC  
              % recruitment classwise
              [rage(2+cl,j,3),page(2+cl,j,3)] = correlate(sAges(ib),mean(clRecruRand(cl,:,:,j),3),...
                  'type','Pearson','partial',sMotion(ib),'NoFigure');
            end
          end
          
          % null checks and visualization
        figure;
          [rr,pp] = correlate(sAges(ib),mean(flexS(~~nodeix,:),1),...
                  'type','Spearman','partial',sMotion(ib),'ExtFigureCmd','RemoveLeg'); hold all;
          correlate(sAges(ib),mean(mean(flexRand(~~nodeix,:,:),1),3),...
                  'type','Spearman','partial',sMotion(ib),'ExtFigureCmd','RemoveLeg');
          title('age v. flexibility overall and null, MC');
          legend(['p = ' num2str(sum(rage(1,:,1)>rr)./nullp)]); hold off;
          sum(page(1,:,1)<0.05);
        figure; hist(rage(1,:,1));
          legend(['r = ' num2str(rr) ', p = ' num2str(sum(rage(1,:,1)>rr)./nullp)]);
          
        figure;
          [rr,pp] = correlate(sAges(ib),mean(mean(recruS(~~nodeix,:,:),3),1),...
                  'type','Spearman','partial',sMotion(ib),'ExtFigureCmd','RemoveLeg'); hold all;
          correlate(sAges(ib),mean(mean(mean(recruRand(~~nodeix,:,:,:),3),1),4),...
                  'type','Spearman','partial',sMotion(ib),'ExtFigureCmd','RemoveLeg');
          title('age v. recruitment overall and null, MC');
          legend(['p = ' num2str(sum(rage(2,:,1)<rr)./nullp)]); hold off;
        figure; hist(rage(2,:,1));
          legend(['r = ' num2str(rr) ', p = ' num2str(sum(rage(2,:,1)<rr)./nullp)]);
            title(['age v. recruitment overall and null, MC']);
        for cl=1:nC
        figure;
          [rr,pp] = correlate(sAges(ib),mean(clRecru(cl,:,:),3),...
                  'type','Spearman','partial',sMotion(ib),'ExtFigureCmd','RemoveLeg'); hold all;
          correlate(sAges(ib),mean(mean(clRecruRand(cl,:,:,:),3),4),...
                  'type','Spearman','partial',sMotion(ib),'ExtFigureCmd','RemoveLeg');
          title(['age v. recruitment ' legClass{cl} ' and null, MC']);
          legend(['p = ' num2str(sum(rage(2+cl,:,1)<rr)./nullp)]); hold off;
        figure; hist(rage(2+cl,:,1));
          legend(['r = ' num2str(rr) ', p = ' num2str(sum(rage(2+cl,:,1)<rr)./nullp)]);
          title(['age v. recruitment ' legClass{cl} ' and null, MC']);
        end
        
      end
      
      if ismember('save_stats',p_plot)
          filenameAge = ['rps_corr_age_all_' stats_tag '.csv'];
          filenameDP = ['rps_corr_DP_all_' stats_tag '.csv'];
          filenameCS = ['rps_corr_CS_all_' stats_tag '.csv'];
          dlmwrite(filenameAge,rpsAgeAll);
          dlmwrite(filenameDP,rpsDPAll);
          dlmwrite(filenameCS,rpsCSAll);
      end
      
      total_conn = zeros(numel(ib),1);
      for ii=1:numel(ib)
          adjs = A(ii).adj;
          total_conn(ii) = 0;
          for T=1:size(adjs,1)
              total_conn(ii) = total_conn(ii) + sum(sum(adjs{T}));
          end
      end
      correlate(sAges(ib),total_conn,'type','Spearman'); % they are correlated
      xlabel('age'); ylabel('total summed connectivity, all timesteps');
      
      correlate(total_conn,mean(flexS));
      xlabel('total connectivity'); ylabel('flexibility');
      
      correlate(total_conn,ncommS(end,:),'type','Spearman');
      xlabel('total connectivity'); ylabel('number of communities');
      
      correlate(sAges(ib),mean(flexS),'partial',total_conn,'type','Spearman');
      xlabel('age'); ylabel('flexibility'); 
      title('total connectivity partialed out');
      
      correlate(sAges(ib),ncommS(end,:),'partial',total_conn,'type','Spearman');
      xlabel('age'); ylabel('number of communities'); 
      title('total connectivity partialed out');
end
