function R = corecruitment(MA,systemByNode,coSystem)
%CORECRUITMENT      corecruitment coefficient
%
%   R = CORECRUITMENT(MA,systemByNode,coSystem) calculates a corecruitment 
%   coefficient for each node of the network. The corecruitment coefficient 
%   of a node (with respect to a cosystem) corresponds to the average 
%   probability that this node is in the same network community as other 
%   nodes from that cosystem. (This quantity reduces to the recruitment
%   coefficient when coSystem==systemByNode; this is the default coSystem.)
%
%   Inputs:     MA,     Module Allegiance matrix, where element (i,j) 
%                       represents the probability that nodes i and j
%                       belong to the same community
%               systemByNode,	vector or cell array containing the system
%                       assignment for each node
%               coSystem,   vector or cell array containing the system for
%                       which to compute co-recruitment for each node.
%                       if scalar, will use the same value for each node;
%                       default value = systemByNode (recruitment coef)
%
%   Outputs:    R,      corecruitment coefficient for each node
%   _______________________________________________
%   Marcelo G Mattar (08/21/2014) 
%   extended by KJ Schlesinger (04/23/2015)

% Initialize output
R = zeros(length(systemByNode),1);

% Make sure the diagonal of the module allegiance is all nan
MA(logical(eye(size(MA)))) = nan;

% set default value of coSystem
if nargin<3
    coSystem = systemByNode;
end

% process char and/or scalar value of coSystem
if ischar(coSystem)
    coSystem = {coSystem};
end 
if isscalar(coSystem)
    coSystem = repmat(coSystem,length(systemByNode),1);
end
if ( iscell(coSystem) && ~iscell(systemByNode) ) || ...
   ( iscell(systemByNode) && ~iscell(coSystem) )
    error('systemByNode and coSystem types must match!');
end

% Calculate the recruitment for each node
if iscell(systemByNode)
    for i=1:length(systemByNode)
        thisSystem = coSystem{i};
        R(i) = nanmean(MA(i,strcmp(systemByNode,thisSystem)));
    end
else
    for i=1:length(systemByNode)
        thisSystem = coSystem(i);
        R(i) = nanmean(MA(i,systemByNode==thisSystem));
    end
end