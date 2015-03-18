function [] = saveNNPPs(Nmatrix,Pmatrix,nregs,person,runnum,scheme,code,pval)

if pval < 1 && pval >=0
    thresh = 'T';
else
    thresh = [];
end

%dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
%    num2str(person+100) '/run' num2str(runnum) '/' scheme '/NN' ...
%    num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code thresh '.txt'],Nmatrix);

dlmwrite(['/Users/kimberly/Desktop/choking_NNs/NN' num2str(person+100) ...
    '_' num2str(nregs) scheme num2str(runnum) code thresh '.txt'],Nmatrix);


if numel(Pmatrix)
% dlmwrite(['/Users/kimberly/Google Drive/choking/' num2str(nregs) 'results/subj_' ...
%     num2str(person+100) '/run' num2str(runnum) '/' scheme '/PP' ...
%     num2str(person+100) '_' num2str(nregs) scheme num2str(runnum) code '.txt'],Pmatrix);

dlmwrite(['/Users/kimberly/Desktop/choking_NNs/NN' num2str(person+100) ...
    '_' num2str(nregs) scheme num2str(runnum) code thresh '.txt'],Pmatrix);

end

end