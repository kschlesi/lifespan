function C = class_corecruitment(MA,systemByNode)
% computes class corecruitment between each pair of systems in systemByNode
% and returns them in an nC x nC matrix C (where nC = number of systems)
% input: MA, module allegiance matrix
%        systemByNode, a vector or cell array of system names by node
% output: C

cList = unique(systemByNode);
nC = numel(cList);

C = zeros(nC);
for c=1:nC
    % compute co-recruitment of each node in system with cosystem c
    cCorec = corecruitment(MA,systemByNode,cList(c));
    for d=1:nC
        % average over nodes in d
        if iscell(systemByNode)
            C(c,d) = nanmean(cCorec(strcmp(systemByNode,cList{d})));
        else
            C(c,d) = nanmean(cCorec(systemByNode==cList(d)));
        end
    end
end

end