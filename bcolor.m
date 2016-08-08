function h = bcolor(inmat,varargin)
% provides a balanced color plot (no row/cols left out) with no edge lines
    if ~ismatrix(inmat)
        error('input matrix must be two-dimensional'); 
    end
    
    pad = mean(mean(inmat));
    if numel(varargin)==2
        varargin{1} = [varargin{1}(:);...
                       varargin{1}(end)+diff(varargin{1}(end-1:end))];
        varargin{2} = [varargin{2}(:);...
                       varargin{2}(end)+diff(varargin{2}(end-1:end))];
    end
    h = pcolor(varargin{:},padarray(inmat,[1 1],pad,'post'));
    set(h, 'EdgeColor', 'none');
    
    % reset axis labels to coincide with centers of appropriate
    % columns/rows
    if any(ismemvar(varargin,'AxisLabelSpace'))
        [axLs1,axLs2] = varargin{find(ismemvar(varargin,'AxisLabelSpace'),1,'first')+1};
    else
        [ylen,xlen] = size(inmat);
        axLs1 = floor(xlen/10);
        axLs2 = floor(ylen/10);
    end
    
    interp = @(vec_) vec_(1:end-1) + diff(vec_)./2;
    skip = @(vec_,skipno) vec_(1:skipno:end);
    if numel(varargin)==2
        set(gca, 'XTick', skip(interp(varargin{1}),axLs1),...
                 'XTickLabel', skip(varargin{1}(1:end-1),axLs1));
        set(gca, 'YTick', skip(interp(varargin{2}),axLs2),...
                 'YTickLabel', skip(varargin{2}(1:end-1),axLs2));
    else
        [ysz,xsz] = size(inmat);
        set(gca, 'XTick', skip(interp(1:1:xsz+1),axLs1),...
                 'XTickLabel', skip(1:1:xsz,axLs1));
        set(gca, 'YTick', skip(interp(1:1:ysz+1),axLs2),...
                 'YTickLabel', skip(1:1:ysz,axLs2));
    end
    
end