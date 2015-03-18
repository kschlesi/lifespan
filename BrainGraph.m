function[] = BrainGraph(data, person, labels, pos, nodedata, range)
%
% Input variables: data is the 1 by E matrix of edge weights(taken in
% order out of the upper triangular part of the node-node matrix), person
% is the (integer) number of the individual on which you are
% interested in plotting the data, nodedata is the 1 by N matrix of node
% weights (optional - put 0 if you just want black), and range is the five
% percentage values you want to see e.g. top 5%, 4%, 3%, 2%, 1% etc (also
% optional - will automatically threshold the top 5% as shown if there are
% no inputs)s
%
% Output: A plot of the brain where edges are the weights in data and nodes
% are 

if nargin < 6
    range = [5 4 3 2 1];
end

if nargin < 5
   nodedata = 0;
end

% labelFile = 'hybrid_labels_long.txt';
% posFile = 'hybrid_centroids.mat';
% labels = importdata(labelFile);
% B = load(posFile);
N = length(labels);
E = N*(N-1)/2;

% Getting the physical location data and initializing an N by 3 matrix of
% x,y,z positions for each of the N brain regions
%pos = squeeze(mean(positions.hybrid_c(person,:,:,:),2))';

% Iinitializing the matrix where we store edge data
%b = triu(ones(N),1);
%b(~~b)=data;

if numel(data) == N
    %sizes = 100*data;
    ecol = zeros(size(data,1),3);
    ecol(data==1,:) = repmat([0 0 1],size(ecol(data==1,:),1),1);
    ecol(data==0.5,:) = repmat([1 1 1],size(ecol(data==0.5,:),1),1);
    ecol(data==1.5,:) = repmat([0 1 0],size(ecol(data==1.5,:),1),1);
    ecol(data==2,:) = repmat([1 0 0],size(ecol(data==2,:),1),1);
    sizes = 100*ones(N,1);
else
    sizes = 100*ones(N,1);
    ecol = repmat([1 1 1],N,1);
end

figure;
% Plotting the nodes: you can also comment out the edge section if your
% just want to see nodes colored based on strength
if nodedata == 0;
    for i = 1:N;
        scatter3(pos(i,1),pos(i,2),pos(i,3),sizes(i),'k','fill',...
            'MarkerEdgeColor',ecol(i,:),'LineWidth', 1.25);
        hold all
    end
else
    for i = 1:N;
        scatter3(pos(i,1),pos(i,2),pos(i,3),sizes(i),nodedata(i),'fill',...
            'MarkerEdgeColor',ecol(i,:),'LineWidth', 1.25);
        hold all
    end
end

% % Plotting edges
% range = sort(range, 'descend');
% A=reshape(b,1,[]);
% bs = sort(A, 'descend');
% n=round(numel(bs)*range(1)/100)+1;
% T1=bs(n);
% n=round(numel(bs)*range(2)/100)+1;
% T2=bs(n);
% n=round(numel(bs)*range(3)/100)+1;
% T3=bs(n);
% n=round(numel(bs)*range(4)/100)+1;
% T4=bs(n);
% n=round(numel(bs)*range(5)/100)+1;
% T5=bs(n);
% 
% tic;
% for i=1:N;
%     for j=1:N;
%         if i>j;
%             if b(j,i) >= T1 && b(j,i) < T2;
%                 a = linspace(pos(i,1),pos(j,1));
%                 c = linspace(pos(i,2),pos(j,2));
%                 d = linspace(pos(i,3),pos(j,3));
%                 plot3(a,c,d,'LineWidth',2,'Color',[1 0.92549 0.545098]);
%                 hold all
%             end
%             if b(j,i) >= T2 && b(j,i) < T3;
%                 a = linspace(pos(i,1),pos(j,1));
%                 c = linspace(pos(i,2),pos(j,2));
%                 d = linspace(pos(i,3),pos(j,3));
%                 plot3(a,c,d,'LineWidth',2,'Color',[1 0.72549 0.058824]);
%                 hold all
%             end
%             if b(j,i) >= T3 && b(j,i) < T4;
%                 a = linspace(pos(i,1),pos(j,1));
%                 c = linspace(pos(i,2),pos(j,2));
%                 d = linspace(pos(i,3),pos(j,3));
%                 plot3(a,c,d,'LineWidth',2,'Color',[1 0.270588 0]);
%                 hold all
%             end
%             if b(j,i) >= T4 && b(j,i) < T5;
%                 a = linspace(pos(i,1),pos(j,1));
%                 c = linspace(pos(i,2),pos(j,2));
%                 d = linspace(pos(i,3),pos(j,3));
%                 plot3(a,c,d,'LineWidth',2,'Color',[0.545098 0 0]);
%                 hold all
%             end
%             if b(j,i) >= T5;
%                 a = linspace(pos(i,1),pos(j,1));
%                 c = linspace(pos(i,2),pos(j,2));
%                 d = linspace(pos(i,3),pos(j,3));
%                 plot3(a,c,d,'LineWidth',2,'Color',[0.2 0.09839 0]);
%                 hold all
%             end
%         end
%     end
% end
grid off
xlabel('x');
ylabel('y');
zlabel('z');

%%%% Then use these commands one at a time in the shell to get different
%%%% views.
% view(0,90); title('Axial View',  'FontWeight','bold');
% 
% view(0,0); title('Coronal View',  'FontWeight','bold');
% 
% view(90,0); title('Sagittal View',  'FontWeight','bold');

%toc;




end