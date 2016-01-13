% rearrangement: visualizing the measures.

%% first: static networks.
% make sure "A" adjacency matrices are loaded from age_analysis.
% set up network & parameters:
subj = 11;
run = 1;
gamma = 1; 
p = 100;
reg_remove = randi(n,[10,1]);

% grab data and create two matrices.
mat_orig = A(subj).adj{run};
[mat_rem, newix] = remove_missings(mat_orig, reg_remove);

% CD on both, p times, compute rearrangement
[C_orig,Q_orig] = genlouvainREPs(mat_orig,p,gamma,0);
[C_rem,Q_rem] = genlouvainREPs(mat_rem,p,gamma,0);
remat = rearrangement(C_orig,C_rem,newix);

% visualise partitions and partition changes.
figure; subplot(1,2,1); bcolor(C_orig'); colorbar; title('original partitions');
        subplot(1,2,2); bcolor(add_missings_dim(C_rem,newix,n,2,0)'); colorbar; title('new partitions');
figure; subplot(1,3,1); bcolor(C_orig(:,newix)'); colorbar; title('orig partitions');
        subplot(1,3,2); bcolor(C_rem'); colorbar; title('new partitions');
        subplot(1,3,3); bcolor(remat); colorbar; title('rearrangement');

% look at distributions over nodes.
figure; plot(sum(remat)./n); title('signed rearrangement per node');
    xlabel('node ID'); ylabel('rearrangement');
figure; plot(sum(abs(remat))./n); title('magnitude rearrangement per node');
    xlabel('node ID'); ylabel('rearrangement');
    
%% ok, now remove each node systematically & find resulting rearrangements.

[C_orig,Q_orig] = genlouvainREPs(mat_orig,p,gamma,0);
rem_nodes = zeros(p,n-1,n);
Q_nodes = zeros(p,n);
rearr_nodes = zeros(n-1,n-1,n);
for ii=1:n
    [mat_i,newix_i] = remove_missings(mat_orig, ii);
    [rem_nodes(:,:,ii),Q_nodes(:,ii)] = genlouvainREPs(mat_i,p,gamma,0);
    rearr_nodes(:,:,ii) = rearrangement(C_orig,rem_nodes(:,:,ii),newix_i);
end

%%
% now look for each node at which other node(s) it causes to rearrange.
% urgh, rearr_nodes is all wrong.
rn_fill = zeros(n,n,n);
for i=1:n
    intermediate = add_missings_dim(rearr_nodes(:,:,i),removeval(1:n,i),n,1,nan);
    rn_fill(:,:,i) = add_missings_dim(intermediate,removeval(1:n,i),n,2,nan);
end
figure; bcolor(squeeze(nanmean(abs(rn_fill),1))'); 
xlabel('ID of removed node'); ylabel('ID of affected node');
title('node-by-node rearrangement due to single node removal');
h = colorbar('southoutside'); h.Label.String = 'abs rearrangement of affected node';

retot = squeeze(nanmean(nanmean(rn_fill,1),2));
reabstot = squeeze(nanmean(nanmean(abs(rn_fill),1),2));

figure;
plot(retot);
hold all;
plot(reabstot);
title('mean rearrangement due to single node removal');
xlabel('ID of removed node');
ylabel('mean rearrangement of all other nodes');
legend('signed rearrangement','abs rearrangement');

figure;
plot(sort(retot));
hold all;
plot(sort(reabstot));
correlate(mean(mean(rearr_nodes,1),2),mean(mean(abs(rearr_nodes),1),2));

% look at node degree -- no apparent correlation
correlate(retot,sum(mat_orig));

% what about in- and out-degree?
[cSize,inND,outND] = comm_degrees(A,partnS,nodeix,n,t,nruns,nsubs);
correlate(retot,inND(:,run,subj));
correlate(retot,outND(:,run,subj));

% participation coefficient? VERY weak anti-correlation
correlate(partcoef(mode(C_orig)',{mat_orig}),retot);





