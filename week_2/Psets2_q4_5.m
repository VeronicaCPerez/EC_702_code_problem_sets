%%%%%%%%%%%%%
% Simulates Solow model in discrete time
% Original code by Stefano Pica Fall 2019 and Hannah Rhodenhiser 2021
% Edited by Jiayin Liu Fall 2023 for Macro Pset2
%%%%%%%%%%%%%

clear; close all; clc

% Model Parameters - annual calibration
s = 0.2;        %average share of real investment in real GDP (around 20%)
delta = 0.05;   %average ratio of depreciation to GDP (around 5%)
n = 0.04;       %population growth (around 2%)
g = 0.04;       %technology growth
alpha = 1/3;    %current level of capital share in the economy (around 33%)

%time objects
T = 100;        %end of time %V review: why 101
time = 1:T;     %grid of time periods

%preallocate vectors     
x = NaN(100,T);       %capital per capita; small x
y = NaN(100,T);       %output per capita; small y
mpk = NaN(100,T);     %marginal product of capital
sharek = NaN(100,T);  %capital share

% generate random variables eps from normal distribution
rng('default') % For reproducibility
eps =normrnd(0,1,[100,1]) % generate epsilon for 100 countries
e_eps = exp(eps)

% SS quantities: analytical solution
% st-st equation
%I=ones(100,1)
x_ss = (s/((g+n+delta))^(1/(1-alpha))); %SS capital per capita
% use a function
y_ss = prodfCD(x_ss, alpha);                       %SS output per capita
mpk_ss =  mpkCD(x_ss, alpha);                   %SS marginal product of capital
sharek_ss = mpk_ss.*x_ss/y_ss;                           %SS capital share   
%first point in time
% initial station, show the relative location to st-st
x(:,1) = e_eps.*x_ss;                 %initial condition on state variable: half of st-st, twice of stst
y(:,1) = prodfCD(x(:,1), alpha);    %output
mpk(:,1) = mpkCD(x(:,1), alpha); %marginal product of capital
sharek(:,1) = mpk(:,1).*x(:,1)./y(:,1);          %capital share

% The actual simulation part
for t = 2:T %solve recursively with the law of motion
    for i = 1:100
       x(i,t) = (1/((1+n+g))).*((1-delta).*x(i,t-1) + s.*x(i,t-1).^alpha);
       y(i,t) = prodfCD(x(i,t), alpha);
       mpk(i,t) = mpkCD(x(i,t), alpha); 
       sharek(i,t) = mpk(i,t).*x(i,t)/y(i,t);
    end
end



% Computing Aggregates per capita
% look at notes for sep15
A0 = 1 ;                              % normalization 
y_pc = y.*(1+g).^time.*A0 ; %gdp per capita; "." means elementwise operation of matrix
k_pc = x.*(1+g).^time.*A0 ; %capital per capita
w = (y-mpk.*x).*(1+g).^time.*A0 ; % wages
R = mpk ;                            %return to capital


% log all the variables that are growing in time
ln_y_pc = log(y_pc) ;
ln_k_pc = log(k_pc) ;
ln_w = log(w) ;

% Potential
ln_k_pot = log(x_ss.*(1+g).^time.*A0) ;
ln_y_pot = log(y_ss.*(1+g).^time.*A0) ;
%potential gdp not in logs
y_pot = exp(ln_y_pot)

ln_R_pot = log(mpk_ss) ;
ln_w_pot = log((y_ss-mpk_ss.*x_ss).*(1+g).^time.*A0) ;



% FIGURE 4 % 
% plot 10, 20 countries explain the patterns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select 5 random countries
num_random_countries = 10; % Change this as needed

% Generate random row indices
random_c_index = randperm(size(y_pc, 1), num_random_countries);

% Select the random rows from the matrix
gdp_subset = y_pc(random_c_index, :);

% make the figure
figure('Color', 'white'); % Create a new figure
hold on; % Allow multiple lines on the same plot
for i = 1:num_random_countries
    plot(time, gdp_subset(i, :),'LineWidth', 2);
end

% Customize the plot as needed (e.g., labels, title, legend)
xlabel('Time', 'FontSize', 14);
ylabel('Value', 'FontSize', 14);
legend('Country 1', 'Country 2', 'Country 3', 'Country 4', 'Country 5','Country 6', 'Country 7', 'Country 8', 'Country 9', 'Country 10', 'Location', 'southeast'); % Customize legend
% Set x-axis limits to restrict the range to 1 to 100
xlim([1, 100]);

% plot potential
% plot(time, y_pot, 'b','LineWidth', lwidnum);

% Hold off to prevent further additions to the current plot
hold off;

% saveas(gcf, '/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS2_graph4.png');
% saveas(gcf, '/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS2_graph4.png');

% POINT 5 % 
% exporting the data for the regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Export the data to a CSV file
writematrix(y_pc, "ps2_q5_simulated_data_2.csv");









% Function definitions 

function y = prodfCD(x, alpha)
% This function recalls the Cobb Douglas production function

y = x.^alpha;
    
end

function mpk = mpkCD(x, alpha)
% This function recalls the expression for the marginal product of capital
% per efficency labor
% for a Cobb Douglas production function

mpk = alpha*x.^(alpha-1);

end
