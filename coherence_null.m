function [p_mcohr,p_ncomm,p_fskew,mcohr,ncomm,fskew] = coherence_null(idREGs,partn,reps,t_plot,alpha)

[n,t] = size(partn);
if nargin<5
    alpha = 0.05;
end
p_fskew = [];
fskew = [];

% create random partitions & find stats
s_partn = zeros(n,t,reps);  % store random partitions
mcohr = zeros(t,reps);      % store max cohr for each partition & slice
ncomm = zeros(t,reps);      % store #comms for idREGs per partn & slice
%fskew = zeros(t,reps);      % store distribution skew per partn & slice 
for i = 1:reps
  for T = 1:t
    s_partn(:,T,i) = partn(randperm(n),T);
    ncomm(T,i) = numel(unique([s_partn(idREGs,T,i);0]))-1;
  end
  mcohr(:,i) = coherence(idREGs,s_partn(:,:,i));
 % fskew(:,i) = partition_skew(idREGs,s_partn(:,:,i));
end

% find real values for each stat
r_mcohr = coherence(idREGs,partn);
r_ncomm = zeros(t,1);
for T=1:t; r_ncomm(T) = numel(unique([s_partn(idREGs,T,i);0]))-1; end;
%r_fskew = partition_skew(idREGs,partn);

% plots
if t_plot
  for T = t_plot
    figure;
        title('max coherence');
        hist(mcohr(T,:),reps/2);
        hold on;
        %plot(r_mcohr(T),reps);
    figure;
        title('max coherence');
        hist(ncomm(T,:),reps/2);
        hold on;
        %plot(r_ncomm(T),reps);
%     figure;
%         title('max coherence');
%         hist(fskew(T,:),reps/2);
%         hold on;
%         %plot(r_fskew(T),reps);
  end
end

% find p-values for each stat
%ix = -1*floor(-1*reps*(1-alpha));
hitail = sum(mcohr > repmat(r_mcohr,1,reps),2)/reps;
lotail = sum(mcohr < repmat(r_mcohr,1,reps),2)/reps;
p_mcohr = [(hitail < alpha), hitail, r_mcohr;...
           (lotail < alpha), lotail, r_mcohr;...
           ((hitail < alpha/2) + (lotail < alpha/2)), hitail+lotail, r_mcohr];
%disp(size(ncomm)); disp(size(r_ncomm)); disp(size(repmat(r_ncomm,1,reps)));
hitail = sum(ncomm > repmat(r_ncomm,1,reps),2)/reps;
lotail = sum(ncomm <= repmat(r_ncomm,1,reps),2)/reps;
p_ncomm = [(hitail < alpha), hitail, r_ncomm;...
           (lotail <= alpha), lotail, r_ncomm;...
           ((hitail < alpha/2) + (lotail <= alpha/2)), hitail+lotail, r_ncomm];
% hitail = sum(fskew > repmat(r_fskew,1,reps),2)/reps;
% lotail = sum(fskew < repmat(r_fskew,1,reps),2)/reps;
% p_fskew = [(hitail < alpha), hitail, r_fskew;...
%            (lotail < alpha), lotail, r_fskew;...
%            ((hitail < alpha/2) + (lotail < alpha/2)), hitail+lotail, r_fskew];      
       
end