% Veronica Perez
% December 5th 2023
% Homework 3 David Lagakos

close all; clear all; clc;

%% 1. define parameters
delta = 0.10;
theta = 0.4;
beta = 0.96;
q = 0.7;
sigma = 0.1;

%% 2. productivity process
% Z in {-sigma,+sigma}, P = [q, 1-q; 1-q, q];

Zgrid = [-sigma; sigma];
Znum = length(Zgrid); 
P = [q , 1-q; 1-q, q];

%% 3. Value function iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vfmaxit = 1000; %max number of VF iterations
vftol = 1e-7;       %tolerance on Bellman equation, in percentages

T=100;                    %final period for transition
transitspan=0:T-1; %time span for transition

% %exact SS quantities
% Abar=(1/beta-1+delta)/alpha;
% Abar;
% A=Abar+0.03;

%exact SS quantities
kss_z_zero=((1/beta-1+delta)/(theta))^(1/(theta-1));
k_z_low = ((theta/exp(-sigma))/(1/beta-1+delta))^(1/(1-theta));
k_z_high = ((theta/exp(sigma))/(1/beta-1+delta))^(1/(1-theta));

%define grid
kmin=kss_z_zero*0.5;
kmax=kss_z_zero*2;
knum=1000;
K0=linspace(kmin,kmax,knum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute return function (utility function). We use meshgrid
[kp,k]=meshgrid(K0);

Glow=zeros(knum, knum);
% RHS of budget constraint
Gmax_z_low = exp(-sigma) .* (k.^(theta)) + ((1-delta).*k);
Gmax_z_high = exp(sigma) .* (k.^(theta)) + ((1-delta).*k);

%restrict control variable to feasible set
kp1 = kp;
kp1(kp1>=Gmax_z_low | kp1<Glow) =NaN;
utility_low = log(exp(-sigma).*(k.^(theta)) + (1-delta).*k - kp1);

kp2 = kp;
kp2(kp2>=Gmax_z_high| kp2<Glow) =NaN;
utility_high = log(exp(sigma).*(k.^(theta)) + (1-delta).*k - kp2);
% In utility matrix, rows represent the state (capital
%today) columns are the control (capital tomorrow).

%initialize the policy and VF arrays
Vold = zeros(Znum*knum,1);   %this will store old VF - max possible utility of the state today 
V = zeros(Znum*knum,1);      %this will store new VF
kprimeind=zeros(2*knum,1); %the best place in the grid of the best one

disp('%%%%')
disp('Starting VFI now')
disp(' ')

RHSMAT = zeros(Znum*knum,knum);
% VFI -- took from sabhya, much better than what i did
for vfit = 1:vfmaxit
    %form Bellman equation
    exp_Vold(1:knum,1) = P(1,1).*Vold(1:knum,1) + P(1,2).*Vold((knum+1):2*knum, 1);
    exp_Vold((knum+1):2*knum,1) = P(2,1).*Vold(1:knum,1) + P(2,2).*Vold((knum+1):2*knum, 1);

    RHSMAT(1:knum, 1:knum) = utility_low + beta*repmat((exp_Vold(1:knum,1))', knum,1);
    RHSMAT((knum+1):2*knum, 1:knum) = utility_high + beta*repmat((exp_Vold((knum+1):2*knum,1))', knum,1);
    [V,kprimeind] = max(RHSMAT,[],2); % a column vector containing max in each row
    absdiff=abs(V-Vold);
    vferr = max(absdiff); %maximum absolute error - one metric of distance
    if (mod(vfit,50)==1)
        disp(['VF error = ' num2str(vferr) ' on VF iteration ' num2str(vfit) '.'])
    end
       %exit if converged
    if (vferr<vftol)
        disp(['VFI complete at iteration ' num2str(vfit) '.'])
        break; 
    end
    %if not converged, then Vold <- V, and redo
    Vold = V;
end

Kpol_zlow = K0(kprimeind(1:knum)); %capital policy function-  capital at the best value
Kpol_zhigh = K0(kprimeind(knum+1:2*knum));

Cpol_zlow = (exp(Zgrid(1)) .* (K0.^theta)) + (1-delta).*K0 - Kpol_zlow;
Cpol_zhigh = (exp(Zgrid(2)) .* (K0.^theta)) + (1-delta).*K0 - Kpol_zhigh;

%% 4. Graph for policy functions

%% Plotting Capital Policy Functions
figure;
set(gcf,'color','w'); % Set background color to white
plot(K0, Kpol_zlow, 'b-', 'LineWidth', 2, 'DisplayName', 'Low Productivity');
hold on;
plot(K0, Kpol_zhigh, 'r-', 'LineWidth', 2, 'DisplayName', 'High Productivity');

% Plotting the 45-degree line
plot(K0, K0, 'k--', 'LineWidth', 1.5, 'DisplayName', '45-degree Line');

title('Capital Policy Functions');
xlabel('Current Capital (k)');
ylabel('Next Period Capital (k'')');
legend('show');
grid on;
hold off;
saveas(gcf,"/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS_11_q_2_2_1.png")


%% Plotting Value Functions
figure;
set(gcf,'color','w'); % Set background color to white

% Plotting the value functions
plot(K0, V(1:knum), 'b-', 'LineWidth', 2, 'DisplayName', 'Low Productivity');
hold on;
plot(K0, V(knum+1:2*knum), 'r-', 'LineWidth', 2, 'DisplayName', 'High Productivity');

title('Value Functions');
xlabel('Current Capital (k)');
ylabel('Value Function (V)');
legend('show');
grid on;
hold off;
saveas(gcf,"/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS_11_q_2_2_2.png")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM THE STEADY STATE BEFORE THE SHOCK %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5. Simulating the shocks
rng(70223); %random number generator seed, for reproducibility
T = 50; 
Zinit = 1;

Zsimind = zeros(T,1);

% draw random numbers on [0,1]
z_shocks = rand(T,1);

% Simulate the z process by comparing uniform shocks to intervals on the transition matrix
Zsimind(1) = Zinit;
for t=2:T
    %extract yesterday's Z
    Zct = Zsimind(t-1);
    
    %compare today's uniform shock to transition matrix intervals
    if (z_shocks(t)<P(Zct,1))
        Zprimect = 1;
    else
        Zprimect = 2;
    end
    
    %store simulated Z
    Zsimind(t) = Zprimect;
end
Zsim = Zgrid(Zsimind);

%%%%%%%%%
% Plot simulated path for Z
%%%%%%%%%
figure; 
scatter(1:T-1,Zsim(1:T-1),'LineWidth',3);
xlabel('t'); ylabel('Z'); title('Simulation of Two-State Markov Chain')
set(gca,'FontSize',14); ylim([-0.2 0.2])


%% 6. Transit for i and c
kstart = 0.1 * kss_z_zero;  % Initial condition, 10 percent of the non-stochastic steady-state capital

[~, kindex] = min(abs(K0 - kstart)); % Find the grid point closer to the initial condition
kinit = K0(kindex);  % Point in the grid closer to kstart as the initial condition

ktransit = zeros(T, 1);  % Preallocate the capital transition vector
itransit = zeros(T, 1);  % Preallocate the investment transition vector
ctransit = zeros(T, 1);  % Preallocate the consumption transition vector

ktransit(1) = kinit;  % First entry of the capital transition is the initial condition

for it = 2:T
    % Update capital using policy functions
    if Zsim(it) == sigma  % High shock
        ktransit(it) = interp1(K0, Kpol_zhigh, kinit, 'linear', 'extrap');
        ctransit(it) = (exp(sigma) .* (kinit.^theta)) + (1-delta).* kinit - ktransit(it);
    else  % Low shock
        ktransit(it) = interp1(K0, Kpol_zlow, kinit, 'linear', 'extrap');
        ctransit(it) = (exp(-sigma) .* (kinit.^theta)) + (1-delta).* kinit - ktransit(it);
    end

    % Update investment
    
    itransit(it) = ktransit(it) - (1 - delta) * kinit;

    % Update kinit
    kinit = ktransit(it);  % Update initial capital value
end

%% Plotting the Transitions (excluding period zero)
figure;
set(gcf, 'color', 'w');  % Set background color to white

subplot(3, 1, 1);
plot(2:T, ktransit(2:end), 'b-', 'LineWidth', 2);  % Exclude period zero
title('Capital Transition');
xlabel('Time Period');
ylabel('Capital (k)');
grid on;

subplot(3, 1, 2);
plot(2:T, ctransit(2:end), 'r-', 'LineWidth', 2);  % Exclude period zero
title('Consumption Transition');
xlabel('Time Period');
ylabel('Consumption (c)');
grid on;

subplot(3, 1, 3);
plot(2:T, itransit(2:end), 'g-', 'LineWidth', 2);  % Exclude period zero
title('Investment Transition');
xlabel('Time Period');
ylabel('Investment (i)');
grid on;

sgtitle('Transitions of Capital, Consumption, and Investment');

saveas(gcf,"/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS_11_q_2_3.png")
