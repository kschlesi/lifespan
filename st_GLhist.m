function [spatioC, tempC] = st_GLhist(C_hist,mnp,varargin)
% this function takes a community GL history and plots (and returns) the
% course of its staiotemporal consolidation. One can then observe these
% as a function of gamma and omega. yayayay.

% set parameters
assert(ndims(C_hist)==4);
[p,n,t,mu] = size(C_hist);
assert(p==length(mnp));

% compute spatial and temporal consolidations over history
spatioC = zeros(mu,p);
tempC = zeros(mu,p);
for m=1:mu
    for pp = 1:p
        comms = squeeze(C_hist(pp,:,:,m));
        spatioC(m,pp) = mean(vnique(comms,2));
        tempC(m,pp) = mean(vnique(comms,1));
    end
end

% plot the curves as a function of mu
if ~any(ismemvar(varargin,'ExtFigureCmd'));
figure;    
end
if any(ismemvar(varargin,'MeanTraj'));
spatioC = mean(spatioC,2);
tempC = mean(tempC,2);
end
plot((1+n-spatioC),(1+t-tempC),'.--');
%legend(num2str((1:pp)'),'Location','NorthWest');
ylabel('mean inter-slice similarity over nodes');
xlabel('mean spatial similarity over slices');

function [proj, u] = vnique(C,dim)
% takes a multidimensional array C and returns the number of
% unique elements projected along any dimension.

% currently supports only 2D arrays.
u = unique(C);
proj = zeros(size(C,dim),1);
for i=u'
    iproj = ~~sum((C==i),-1*(dim-3));
    proj = proj + iproj(:);
end

end

end