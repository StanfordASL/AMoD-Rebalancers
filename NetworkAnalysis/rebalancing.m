% clear all;
% close all;

MODE = 2;
PLOT = 1;
DRIVEN = 0;
DATAFILE = 'timeVaryingDataSeattle121'
%DATAFILE = 'timeVaryingSimDataNY20'
% load data
% load('sysData');
if MODE == 1
    %load('data2')
    load('testData5Station')
    lambda = lambda_i;
elseif MODE == 2
    %load('NYData20Peak')
    load(DATAFILE)
    hour = 2;
    lambda = lambda_tv(:,hour);
    pij = pij_tv(:,:,hour);
    Tij = Tij_tv(:,:,hour);
elseif MODE == 3
    n = 20;
    rand('state',1);
    A = 5;    % length of environment
    B = 5;    % width of environment
    n = 20;     % number of stations

    lambdaMax = 0.05;
    lambdaMin = 0.025;

    % Generating station positions
    pos_i = rand(n, 2);
    pos_i(:,1) = pos_i(:,1).*A;
    pos_i(:,2) = pos_i(:,2).*B;
    % plot(pos_i(:,1), pos_i(:,2), 'o')

    % Generate Tij
    speed = 0.2;    % vehicles move at this speed per time step
    Tij = zeros(n,n);
    for i = 1:n
        for j = 1:n
            % distanceFunction(pos1, pos2, dist_type)
            % dist_type = 1 -> Euclidean distance
            % dist_type = 2 -> Manhattan distance
            Tij(i, j) = distanceFunction(pos_i(i,:), pos_i(j,:), 1)/speed;
        end
    end

    % Generate arrival rates
    % lambda = (lambdaMax-lambdaMin).*rand(n, 1) + lambdaMin;
    lambda = lambdaMax.*rand(n, 1);
    lambdaAvg = sum(lambda)/n;


    % departure Rate
    % mu = 2*max(lambda_i);
    % mu = 100;

    % Routing Probabilities
    % pij = 3 + rand(n);
    pij = rand(n);
    for i = 1:n
        pij(i,i) = 0;
        psum = sum(pij(i, :));
        pij(i,:) = pij(i,:)./psum;  % normalize to probability distribution
    end
end

if DRIVEN == 0
    fraction = 100;
    mTotal = 1064;
    m = round(fraction*(mTotal/(fraction+1)));
    %m2 = round(m/fraction);
    m2 = mTotal - m;
    m
end
n = length(lambda);
% plot unbalanced system
if PLOT
    [Run, LSun, LIun, Atun] = MVA_exact(lambda,pij,Tij,m);
    % plot unbalanced
    figure(1);
    hold on;
    for i = 1:n
        plot((1:m), Atun(i,:), 'Color',[i/n/2,0,i/n])
    end
    xlabel('Number of vehicles')
    ylabel('Probability at least 1 vehicle is available')
    title('Unbalanced System')
end

%%
% implement rebalancing

sum_ljpji = zeros(n,1);
for i = 1:n
    for j = 1:n
        if j ~= i
            sum_ljpji(i) = sum_ljpji(i) + lambda(j)*pij(j,i);
        end
    end
end
% D_i is the rate of depletion of vehicles from station i due to customers
D_i = -lambda + sum_ljpji;

lipij = zeros(n);
for i = 1:n
    for j = 1:n
        lipij(i,j) = lambda(i)*pij(i,j);
    end
end



% optimization for rebalancing people -- bb is for rebalancing people
cvx_begin
    variable bb(n,n)
    minimize (sum(sum(Tij.*bb)));
    subject to
        sum((bb - bb'),2) == -D_i;
        bb >= 0;
        bb <= lipij;
        for i = 1:n
            bb(i,i) == 0;
        end
cvx_end

Tb_star = cvx_optval;

eta = zeros(n,n);
lambdahat = sum(bb, 2);
tol = 1e-5;
for i = 1:n
    for j = 1:n
        if lambdahat(i) < tol
            eta(i,j) = 1/(n-1);
        else
            eta(i,j) = bb(i,j)./lambdahat(i);
        end
    end
    eta(i,i) = 0;
end

mu1 = lambda - lambdahat;
pijhat = zeros(n,n);
qi = mu1./lambda;
for i = 1:n
    for j = 1:n
        pijhat(i,j) = pij(i,j)/qi(i) - (1-qi(i))/qi(i)*eta(i,j);
    end
    pijhat(i,i) = 0;
end

%m = 1200;
[Rb, LSb, LIb, Atb] = MVA_exact(mu1, pijhat, Tij, m);


% optimization for rebalancing vehicles -- aa is for vehicles
% mu2 is for taxis
cvx_begin
    variable aa(n,n)
    minimize(sum(sum(Tij.*aa)));
    subject to
        sum((aa - aa'),2) == D_i;
        aa >= 0;
        for i = 1:n
            aa(i,i) == 0;
        end
cvx_end

Ta_star = cvx_optval;
si = zeros(n,n);
psi = sum(aa, 2);
for i = 1:n
    for j = 1:n
        if psi(i) < tol
            si(i,j) = 1/(n-1);
        else
            si(i,j) = aa(i,j)./psi(i);
        end
    end
    si(i,i) = 0;
end
mu2 = lambdahat + psi;
pij2 = zeros(n,n);
for i = 1:n
    for j = 1:n
        pij2(i,j) = eta(i,j)*lambdahat(i)/mu2(i) + psi(i)/mu2(i)*si(i,j);
    end
    pij2(i,i) = 0;
end




[Ra, LSa, LIa, Ata] = MVA_exact(mu2, pij2, Tij, m2);

% now to determine the real availabilities for all REAL passengers
through1 = LSb./Rb;
through2 = LSa./Ra;
%fraction = 4; % let's have 3 times as many empty vehicles as taxis (1/4 number of rebalancers as vehicles)

realA = zeros(n,mTotal);
%realA = zeros(n,m+floor(m/fraction));
j = 1; % index for cars (number of cars)
k = 1; % index for rebalancers (number of rebalancers)
firstTime = 1;
count = 1;
for i = 2:size(realA,2)
    numSys1 = round(fraction*(i/(fraction+1)));
    numSys2 = i - numSys1;
    if numSys2 == 0
        realA(:,i) = (through1(:,numSys1))./lambda;
    else
        realA(:,i) = (through1(:,numSys1) + (lambdahat./(lambdahat+psi)).*through2(:,numSys2))./lambda;
    end
    
%     realA(:,i) = (through1(:,j) + (lambdahat./(lambdahat+psi)).*through2(:,k))./lambda;
%     if count > fraction
%         count = 1;
%         k = k + 1;
%     else
%         j = j + 1;
%         count = count + 1;
%     end
end

%%
% Now do a aMOD and compare
alpha = zeros(n,n);
psi = sum(aa, 2);
for i = 1:n
    for j = 1:n
        if psi(i) < tol
            alpha(i,j) = 1/(n-1);
        else
            alpha(i,j) = aa(i,j)./psi(i);
        end
    end
end
lambda_new = lambda + psi;
pijauto = zeros(n,n);
pii = psi./lambda_new;
for i = 1:n
    for j = 1:n
        pijauto(i,j) = alpha(i,j)*pii(i) + pij(i,j)*(1-pii(i));
    end
    pijauto(i,i) = 0;
end
m3 = size(realA, 2);
sum(pijauto,2);
[Rauto, LSauto, LIauto, Atauto] = MVA_exact(lambda_new, pijauto, Tij, m3);

%% plot things
if PLOT
    figure(2);
    hold on;
    for i = 1:n
        plot((1:m), Atb(i,:), 'Color',[i/n/2,0,i/n])
    end
    xlabel('Number of empty vehicles')
    ylabel('Probability at least 1 vehicle is available')
    title('Balanced System - empty vehicles')

    figure(3);
    hold on;
    for i = 1:n
        plot((1:m2), Ata(i,:), 'Color',[i/n/2,0,i/n])
    end
    xlabel('Number of rebalancers/taxis')
    ylabel('Probability at least 1 vehicle is available')
    title('Balanced System - taxis')

    figure(4);
    hold on;
    plot((1:m3), Atauto(1,:), 'r', 'LineWidth',2)

    %plot((1:m3)*(fraction+1), Atauto(1,:), 'r', 'LineWidth',3)

    for i = 1:n
        plot((1:m3), realA(i,:), 'Color',[i/n/2,0,i/n])
    end
    xlabel('Number of total vehicles')
    ylabel('Probability at least 1 vehicle is available')
    title('Availabilities of Real Passengers')
    axis([0,m3,0,1])
end

pp = zeros(n,n);
for i = 1:n
    for j = 1:n
        pp(i,j) = lipij(i,j) - bb(i,j);
    end
    pp(i,i) = 0;
end

save('RebParameters','aa','bb','pp');
%save('RebParameters5Station','qi','pijhat','eta','psi','si');