%%%%%%%%%%%%%
% Simulates Solow model in discrete time with 100 countries
% Original code by Stefano Pica Fall 2019
% Edits by Hannah Rhodenhiser Fall 2021
%%%%%%%%%%%%%

clear; close all; clc

% Model Parameters - annual calibration
s     = 0.2;    % average share of real investment in real GDP (around 20%)
delta = 0.05;   % average ratio of depreciation to GDP (around 5%)
alpha = 1/3;    % current level of capital share in the economy (around 33%)
n = 0.01;       % population growth (around 2%)
g = 0.02;       % technology growth

% Number of countries
num_countries = 100;

% Time objects
T = 100;        % end of time
time = 1:T;     % grid of time periods

% Standard deviation for ε and η
sigma_epsilon = 1;
sigma_eta = 0.5;

% Preallocate arrays for each country
x = NaN(num_countries, T);     % capital per effective unit of labor
y = NaN(num_countries, T);     % output per effective unit of labor
mpk = NaN(num_countries, T);   % marginal product of capital
sharek = NaN(num_countries, T);% capital share

% Preallocate A0 vector for each country
A0 = NaN(num_countries, 1);

% Loop over countries
for country = 1:num_countries
    % SS quantities: analytical solution
    x_ss = (s/((1+g)*(1+n) - (1-delta)))^(1/(1-alpha)); % SS capital per effective unit of labor
    y_ss = prodfCD(x_ss, 1, 1, alpha);                  % SS output per effective unit of labor
    mpk_ss =  mpkCD(x_ss, 1, 1, alpha);                 % SS marginal product of capital
    sharek_ss = mpk_ss*x_ss/y_ss;                      % SS capital share
    
    rng(country+1) % For reproducibility
    % Generate random ε (drawn from normal distribution with σ = 1)
    epsilon = sigma_epsilon * randn(1, 1);
    
    % Generate random η (drawn from normal distribution with σ = 0.5)
    eta = sigma_eta * randn(1, 1);
    
    % Calculate A0 for each country
    A0(country) = exp(eta);
    
    % Initial condition for each country
    x(country, 1) = x_ss * exp(epsilon);
    y(country, 1) = prodfCD(x(country, 1), 1, 1, alpha);
    mpk(country, 1) = mpkCD(x(country, 1), 1, 1, alpha);
    sharek(country, 1) = mpk(country, 1) * x(country, 1) / y(country, 1);
    
    % Simulate growth for each country
    for t = 2:T
        x(country, t) = (1/((1+n)*(1+g)))*((1-delta)*x(country, t-1) + s*y(country, t-1));
        y(country, t) = prodfCD(x(country, t), 1, 1, alpha);
        mpk(country, t) = mpkCD(x(country, t), 1, 1, alpha);
        sharek(country, t) = mpk(country, t) * x(country, t) / y(country, t);
    end
end

% Compute aggregates per capita for each country
y_pc = y .* (1+g).^time .* A0; % gdp per capita
k_pc = x .* (1+g).^time .* A0; % capital per capita
w = (y - mpk .* x) .* (1+g).^time .* A0; % wages
R = mpk; % return to capital

% Log all the variables that are growing in time for each country
ln_y_pc = log(y_pc);
ln_k_pc = log(k_pc);
ln_w = log(w);

% Potential for each country
ln_k_pot = log(x_ss .* (1+g).^time .* A0);
ln_y_pot = log(y_ss .* (1+g).^time .* A0);
ln_R_pot = log(mpk_ss);
ln_w_pot = log((y_ss - mpk_ss .* x_ss) .* (1+g).^time .* A0);

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
legend('Country 1', 'Country 2', 'Country 3', 'Country 4', 'Country 5','Country 6', 'Country 7', 'Country 8', 'Country 9', 'Country 10', 'Location', 'northwest'); % Customize legend
% Set x-axis limits to restrict the range to 1 to 100
xlim([1, 100]);

% plot potential
% plot(time, y_pot, 'b','LineWidth', lwidnum);

% Hold off to prevent further additions to the current plot
hold off;

saveas(gcf, '/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS2_graph6.png');
saveas(gcf, '/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez/figures/PS2_graph6.png');

% POINT 5 % 
% exporting the data for the regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Export the data to a CSV file
writematrix(y_pc, "ps2_q5_simulated_data_2.csv");


% Function definitions 

function y = prodfCD(K, L, A, alpha)
% This function recalls the Cobb Douglas production function

y = K.^alpha*(A*L)^(1-alpha);
    
end

function mpk = mpkCD(K, L, A, alpha)
% This function recalls the expression for the marginal product of capital
% for a Cobb Douglas production function

mpk = alpha*K.^(alpha-1)*(A*L)^(1-alpha);

end
