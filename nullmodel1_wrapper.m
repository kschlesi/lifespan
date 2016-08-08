function output = nullmodel1_wrapper(myPwd,px,k,totalsubjs,missing_subjs,nruns,ts,ts_run,goRange,tag,p)
% for each shuffled subject, re-compute the community structure from the
% multislice shuffled connectivity matrix... compute stats... and save.

NNinpath = [myPwd '/NNs/'];
CDinpath = [myPwd '/CD/'];
extra_tag = '_nullspace';

  disp(['iteration ' num2str(P)]);
  tic;
  load([NNinpath 'Anull' tag extra_tag num2str(px) '.mat']);
  [t,~,ib] = load_CD_results(totalsubjs,missing_subjs,nruns,ts,ts_run,goRange,tag,CDinpath);
  output = age_communities_nulls(ib(k),Anull(k).adj,...
                            NNinpath,CDinpath,p,t,ts,ts_run,...
                            1,0,goRange,[],...
                            [tag extra_tag num2str(px)]);
  [~,~,~,~,Cplotall] = load_CD_results(totalsubjs,missing_subjs,nruns,...
                       ts,ts_run,goRange,[tag extra_tag num2str(px)],CDinpath,k);
  %if size(Cplotall,1)==1
  %    goix = 1;
  %end                      
  %partnRand = zeros(n,t); % find final partitions
  %flexRand = zeros(n);    % compute flexibilities
  ncommRand = zeros(t+1); % number of communities in each slice & overall
      origs = squeeze(Cplotall);
      partnRand = squeeze(mode(origs,1));
      flexRand = flexibility(partnRand',1);
      for T=1:t
          ncommRand(T) = numel(removeval(unique(partnRand(:,T)),0));
      end
      ncommRand(end) = max(max(partnRand(:,:)));
% save
  dlmwrite([CDinpath 'ageRandStats' num2str(ib(k)) tag extra_tag num2str(px) '_partn.txt'],partnRand);
  dlmwrite([CDinpath 'ageRandStats' num2str(ib(k)) tag extra_tag num2str(px) '_flex.txt'],flexRand);
  dlmwrite([CDinpath 'ageRandStats' num2str(ib(k)) tag extra_tag num2str(px) '_ncomm.txt'],ncommRand);
  toc;  
end