%%% removing vision -- plots comparing datasets.

clear all;
load('lifespan_rv.mat');
load('officer_rv.mat');
[n,t_OF,nsubs_OF] = size(partnS_rv_OF);
[~,t_LS,nsubs_LS] = size(partnS_rv_LS);

% load atlas classifications for n=194 hybrid atlas
fid = fopen('hybrid_labels_long.txt');
regName = textscan(fid,'%s',n,'Delimiter','\n');
regName = regName{:};
fclose(fid);
[regClass,regCN] = region_class_map(regName,[pwd '\']);
nodeix = ~ismemvar(regClass,'Visual');

% visualize changes in MEAN flexibility.
[sortedflex_LS,sortix_LS] = sort(mean(flex_rv_LS(~~nodeix,:),1));
[sortedflex_OF,sortix_OF] = sort(mean(flex_rv_OF(~~nodeix,:),1));
[sortedvflex_LS,sortvix_LS] = sort(mean(flex_v_LS(~~nodeix,:),1));
[sortedvflex_OF,sortvix_OF] = sort(mean(flex_v_OF(~~nodeix,:),1));
figure; plot(sortedflex_LS); hold all;
        plot(sortedflex_OF,'r'); hold all;
        plot(sortedvflex_LS,'--b');
        plot(sortedvflex_OF,'--r');
legend('lifespan, no vision','officer, no vision',...  
    'lifespan, vision','officer, vision',...
    'Location','NorthWest');
xlabel('subject ID (sorted)'); ylabel('mean flexibility');

% T-TESTS of this.
% (1) between vis and novis for officer
[t1,p1] = ttest(mean(flex_rv_OF(~~nodeix,:),1),mean(flex_v_OF(~~nodeix,:),1));
% (2) between vis and novis for lifespan
[t2,p2] = ttest(mean(flex_rv_LS(~~nodeix,:),1),mean(flex_v_LS(~~nodeix,:),1));
% (3) between vis-novis officer and vis-novis lifespan
[t3,p3] = ttest2( mean(flex_v_OF(~~nodeix,:),1)-mean(flex_rv_OF(~~nodeix,:),1),...
                 mean(flex_v_LS(~~nodeix,:),1)-mean(flex_rv_LS(~~nodeix,:),1));

% visualize changes in VAR flexibility.
[sortedflex_LS,sortix_LS] = sort(var(flex_rv_LS(~~nodeix,:)));
[sortedflex_OF,sortix_OF] = sort(var(flex_rv_OF(~~nodeix,:)));
[sortedvflex_LS,sortvix_LS] = sort(var(flex_v_LS(~~nodeix,:)));
[sortedvflex_OF,sortvix_OF] = sort(var(flex_v_OF(~~nodeix,:)));
figure; plot(sortedflex_LS); hold all;
        plot(sortedflex_OF,'r'); hold all;
        plot(sortedvflex_LS,'--b');
        plot(sortedvflex_OF,'--r');
legend('lifespan, no vision','officer, no vision',...  
    'lifespan, vision','officer, vision',...
    'Location','SouthEast');
xlabel('subject ID (sorted)'); ylabel('variance in flexibility');

% T-TESTS of this.
% (1) between vis and novis for officer
[t1v,p1v] = ttest(var(flex_rv_OF(~~nodeix,:)),var(flex_v_OF(~~nodeix,:)));
% (2) between vis and novis for lifespan
[t2v,p2v] = ttest(var(flex_rv_LS(~~nodeix,:)),var(flex_v_LS(~~nodeix,:)));
% (3) between vis-novis officer and vis-novis lifespan
[t3v,p3v] = ttest2(var(flex_v_OF(~~nodeix,:))-var(flex_rv_OF(~~nodeix,:)),...
                 var(flex_v_LS(~~nodeix,:))-var(flex_rv_LS(~~nodeix,:)));

        
% make some plots
figure; j=0;
for i=randi(nsubs_LS,[9,1])'; 
    j = j+1;
    subplot(3,3,j); bcolor(partnS_rv_LS(~~nodeix,:,i)); colorbar; 
end
        