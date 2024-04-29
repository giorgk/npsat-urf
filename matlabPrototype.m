%% Program Options
alpha = 0.32;
beta = 0.83;
minStep = 1;


%% First we read the streamlines.
S = readICHNOSgather('test_data/teststrmlinfit.traj');
%% Simulate one streamline
ii = 1;
pp = S(ii,1).p;
vv = S(ii,1).v;
%% 
% Calculate cumulative length
xp = zeros(size(pp,1),1);
vp = zeros(size(pp,1),1);
for k = 1:size(pp,1)-1
    stepLen = sqrt(sum((pp(k,:) - pp(k+1,:)).^2));
    if stepLen > minStep
        xp(k+1) = xp(k) + stepLen;
    end
end




