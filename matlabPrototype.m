%% Program Options
alpha = 0.32;
beta = 0.83;
Dm= 1.1578e-004;
minStep = 0.01;
minElemSize = 20;
lambda = 0; % or decay
K_d=0;
rho_b=1;%it doesnt matter since K_d=0

TimeStep = 1; % this is years
TotalTime = 800; %years
wmega=0.5;%crank Nickolson scheme
%% Should this be equal to porocity?
theta = 0.1;

%% First we read the streamlines.
S = readICHNOSgather('test_data/teststrmlinfit.traj');
%% Simulate one streamline
ii = 1;
pp = S(ii,1).p;
vv = S(ii,1).v;
%% Reverse the data
% The first point is near the well so we have to reverse the data
pp = flipud(pp);
vv = flipud(vv);
%% 
% Calculate cumulative length
xp = zeros(size(pp,1),1);
vp = zeros(size(pp,1),1);
vp(1,1) = vv(1,1);
for k = 1:size(pp,1)-1
    stepLen = sqrt(sum((pp(k,:) - pp(k+1,:)).^2));
    if stepLen > minStep
        xp(k+1) = xp(k) + stepLen;
        vp(k+1) = vv(k+1);
    end
end
%% Create the 1D mesh respecting the minimum discretization
p_1d = xp(1);
v_el = [];
for k = 1:length(xp)-1
    Dx=xp(k+1)-xp(k);
    if Dx > minElemSize
        Nseg=round(Dx/minElemSize);
        x_temp=linspace(xp(k),xp(k+1),Nseg+1)';
        Np_seg=length(x_temp);
        p_1d=[p_1d;x_temp(2:end)];
        v_el = [v_el; vp(k)*ones(length(x_temp)-1,1)];
    else
        p_1d=[p_1d;xp(k+1)];
        v_el = [v_el; vp(k)];
    end
end
%% Assemble matrix
Np = size(p_1d,1);
%Dglo = sparse(Np,Np);
%Mglo = sparse(Np,Np);
Dglo = zeros(Np,Np);
Mglo = zeros(Np,Np);

aL = alpha*p_1d(end)^beta;

for ii = 1:Np-1
    Del = aL*v_el(ii) + Dm;
    Lel = p_1d(ii+1) - p_1d(ii);
    %rho_term = theta + rho_b*K_d;

% %   Dx*theta   | 1  -1|    Vx  |-1  1|                          L |2  1|
% %   --------   |      | + ---- |     | +lambda(theta + pho*Kd)*---|    |
% %      L       |-1   1|     2  |-1  1|                          6 |1  2|
% 
%     D11 =  Del*theta/Lel  -  v_el(ii)/2 + 2*lambda*rho_term*Lel/6;
%     D12 = -Del*theta/Lel  +  v_el(ii)/2 + 1*lambda*rho_term*Lel/6;
%     D21 = -Del*theta/Lel  -  v_el(ii)/2 + 1*lambda*rho_term*Lel/6;
%     D22 =  Del*theta/Lel  +  v_el(ii)/2 + 2*lambda*rho_term*Lel/6;

% %                    L  |2  1|
% %  (pho*Kd + theta) --- |    |
% %                    6  |1  2|
%     A11 = 2*rho_term*Lel/6;
%     A12 = 1*rho_term*Lel/6;
%     A21 = 1*rho_term*Lel/6;
%     A22 = 2*rho_term*Lel/6;
    
    rho_term = 1 + rho_b*K_d/theta;

%      Dx   | 1  -1|    Vx  |-1  1|                         L |2  1|
%   --------|      | + ---- |     | + lambda(1 + pho*Kd/n)*---|    |
%      L    |-1   1|    2n  |-1  1|                         6 |1  2|
   
    D11 =  Del/Lel  -  v_el(ii)/2 + 2*lambda*rho_term*Lel/6;
    D12 = -Del/Lel  +  v_el(ii)/2 + 1*lambda*rho_term*Lel/6;
    D21 = -Del/Lel  -  v_el(ii)/2 + 1*lambda*rho_term*Lel/6;
    D22 =  Del/Lel  +  v_el(ii)/2 + 2*lambda*rho_term*Lel/6;

%                 L  |2  1|
%  (1+ pho*Kd/n) --- |    |
%                 6  |1  2|
    A11 = 2*rho_term*Lel/6;
    A12 = 1*rho_term*Lel/6;
    A21 = 1*rho_term*Lel/6;
    A22 = 2*rho_term*Lel/6;


    % This is prototype and we dotn care about efficiency. Otherwise use
    % sparse function
    Dglo(ii,ii) = Dglo(ii,ii) + D11;
    Dglo(ii,ii+1) = Dglo(ii,ii+1) + D12;
    Dglo(ii+1,ii) = Dglo(ii+1,ii) + D21;
    Dglo(ii+1,ii+1) = Dglo(ii+1,ii+1) + D22;

    Mglo(ii,ii) = Mglo(ii,ii) + A11;
    Mglo(ii,ii+1) = Mglo(ii,ii+1) + A12;
    Mglo(ii+1,ii) = Mglo(ii+1,ii) + A21;
    Mglo(ii+1,ii+1) = Mglo(ii+1,ii+1) + A22;
end
%% Apply boundary conditions and solve the system
T = (0:TimeStep:TotalTime)'*365;
Cinit=zeros(Np,1);
%F=sparse(Np,1);
F = zeros(Np,1);
% Concentration boundary conditions
% We use a constant concentration of 1  
CB = nan(Np,1);
CB(1,1) = 1;

Dt=diff(T);
C=zeros(length(Dt),length(CB));

ldnans=find(isnan(CB));
cnstHD=find(~isnan(CB));


Aglo=Mglo+wmega*Dt(1)*Dglo;
Bglo1 = (Mglo-(1-wmega)*Dt(1)*Dglo);
Bglo1red = Bglo1(2:end,:);

KK=Aglo(ldnans,ldnans);
GG=Aglo(ldnans,cnstHD);

for it=1:length(Dt)
    %Aglo=Mglo+wmega*Dt(it)*Dglo;
    %Bglo=(Mglo-(1-wmega)*Dt(it)*Dglo)*Cinit+Dt(it)*((1-wmega)*F+wmega*F);
    %Bglo=(Mglo-(1-wmega)*Dt(it)*Dglo)*Cinit;%+Dt(it)*((1-wmega)*F+wmega*F);
    RHS1 = Bglo1red*Cinit;
    
    % Apply boundary conditions
    
    %dd=CB(cnstHD);
    %CB(ldnans)=KK\(Bglo(ldnans)-GG*dd);
    CB(ldnans)=KK\(RHS1-GG);
    Cinit=CB;
    C(it,:)=CB';
end
val=C(:,end)';
negval=-1*val;
negval=[zeros(1,1) negval];
negval=negval(1,1:size(val,2));
URF=val+negval;
%% Lognormal Fit procedure in Matlab
% Fit function
ft = fittype( '(1/(x*b*sqrt(2*pi)))*exp((-(log(x)-a)^2)/(2*b^2))', 'independent', 'x', 'dependent', 'y' );
%ft = fittype( 'sqrt(b/(2*pi*x^3))*exp(-(b*(x-a)^2)/(2*a^2*x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.499112701571756 0.336174051321482];
%% Fitting
X = 1:length(T(:,1));
[xData, yData] = prepareCurveData( X, T(:,1)' );
[fitresult, gof] = fit( xData, yData, ft, opts );
ab = coeffvalues(fitresult);
%% lognormal as inline function
lgnrmfun = @(a,b,x)((1./(x.*b.*sqrt(2*pi))).*exp((-(log(x)-a).^2)./(2*b^2)));
%%
plot(T(:,1))
hold on
plot(lgnrmfun(ab(1),ab(2),1:267))












