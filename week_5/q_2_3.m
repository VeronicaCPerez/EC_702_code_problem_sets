%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the Planner's problem using VFI %
% Matt Hong and Pascual Restrepo, EC 702 Fall 2020
% Edits Hannah Rhodenhiser Fall 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;

%reference: Stokey_Lucas_1989 5.1
% model to be solved
% V(k) = max_{kp<=f(k)+(1-delta)k} [ u(c) + beta V(kprime) ]
% BC: c + kprime = f(k) + (1-delta)*k
% prod function: f(k)=k^alpha*ell^(1-alpha)
% utility function: u(c) = c^(1-gamma) / (1-gamma)

%start tracking timing
tic; 
disp('%%%%%%%%%%%%%%%%%%%')
disp('Solving the planners problem')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS BLOCK SETS UP THE MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vfmaxit = 1000; %max number of VF iterations
vftol = 1e-7;       %tolerance on Bellman equation, in percentages

T=40;                    %final period for transition
transitspan=0:T-1; %time span for transition


%model parameters (note, this is parametrized in a way such that one unit
%of time is a year)
gamma=1.5;  %this means that there is a low willingness to subs c over time
delta=0.06;
alpha=0.3;
beta=0.95;
ell=1;      %normalization

%exact SS quantities
Abar=(1/beta-1+delta)/alpha;
Abar;
A=0.5;

%exact SS quantities
kss=A*(alpha/(1/beta-1+delta))^(1/(1-alpha))*ell;
Rss=A*alpha*(kss/ell)^(alpha-1);
Wss=A*(1-alpha)*(kss/ell)^(alpha);
rss=1/beta-1;

%define grid
kmin=kss*0.5;
kmax=kss*2;
knum=1000;
K0=linspace(kmin,kmax,knum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute return function (utility function). We use meshgrid
[kp,k]=meshgrid(K0);
Glow=zeros(knum, knum);
Gmax=k.^alpha*ell^(1-alpha)+(1-delta)*k;
kp(kp>=Gmax | kp<Glow)=NaN; %restrict control variable to feasible set
utility=(k.^alpha*ell^(1-alpha)+(1-delta)*k-kp).^(1-gamma)/(1-gamma);
% In utility matrix, rows represent the state (capital
%today) columns are the control (capital tomorrow).

%initialize the policy and VF arrays
Vold = zeros(knum,1);   %this will store old VF
V = zeros(knum,1);      %this will store new VF
kprimeind=zeros(knum,1);

disp('%%%%')
disp('Starting VFI now')
disp(' ')

for vfit = 1:vfmaxit
    %form Bellman equation
    RHSMAT = utility + beta*repmat(Vold',knum,1);
   	[V,kprimeind] = max(RHSMAT,[],2); % a column vector containing max in each row
    absdiff=abs(V-Vold);
    vferr = max(absdiff); %maximum absolute error
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
toc;
disp('%%%%')
disp(' ')

Kpol = K0(kprimeind); %capital policy function
Cpol = K0.^alpha*ell^(1-alpha) + (1-delta)*K0 - Kpol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT POLICY FUNCTIONS AND VERIFY THEY MATCH THE SLIDES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lwidnum = 2;    %line width on graphs
fsizenum = 14;  %font size on graphs

policy=figure; 
subplot(1,3,1)
plot(K0,V,'b','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Value');
title('Value Function')
set(gca,'FontSize',fsizenum)

subplot(1,3,2)
plot(K0,Kpol,'r',K0,K0,'k--', 'LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Capital Tomorrow')
title(' Capital Policy')
legend('Capital Policy','45 degree line','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(1,3,3)
plot(K0,Cpol,'k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Consumption');
title('Consumption Policy')
set(gca,'FontSize',fsizenum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM AN ARBITRARY INITIAL CONDITION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kstart=0.5*kss;                 %initial condition---half of steady state
[~,kindex]=min(abs(K0-kstart)); %grid point closer to initial condition
kinit=K0(kindex);               %point in grid closer to kstart as initial condition

ktransit=zeros(T,1); %preallocate the capital transition vector
ktransit(1)=kinit;   %first entry of capital transition is initial condition
for it=2:T
    ktransit(it)=Kpol(K0==kinit); %get from the grid the initial value and apply the PF
    kinit=ktransit(it);           %update initial capital value
end

%transition for other interesting equilibrium objects
Wtransit=(1-alpha)*(ktransit/ell).^(alpha); %transition for wages
Rtransit=alpha*(ktransit/ell).^(alpha-1);   %transition for Rental rate
rtransit=Rtransit-delta;                    %transition for Rental rate

transition=figure;
subplot(2,2,1)
plot(transitspan,ktransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[kss,kss],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Capital');
title('Transition for Capital')
legend('Capital','KSS','Location','SouthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,2)
plot(transitspan,Rtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[Rss,Rss],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Rental Rate');
title('Transition for Rental Rate')
legend('Rental Rate','RSS','Location','NorthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,3)
plot(transitspan,rtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[rss,rss],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Interest Rate');
title('Transition for Interest Rate')
legend('Interest Rate','rSS','Location','NorthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,4)
plot(transitspan,Wtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[Wss,Wss],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Real Wage');
title('Transition for Real Wage')
legend('Real Wage','Wss','Location','SouthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% TRANSITIONAL DYNAMICS FOLLOWING A SHOCK %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% THIS IS NEW MATERIAL                    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note:  in this example, we consider an increase in population by 10%
%new model parameters (this is the place to include your shocks)    %
%note also that after shocks, VF changes, so we need to recompute it%
ell=1.1;

%new exact SS quantities
kss_new=(alpha/(1/beta-1+delta))^(1/(1-alpha))*ell;
Rss_new=alpha*(kss_new/ell)^(alpha-1);
Wss_new=(1-alpha)*(kss_new/ell)^(alpha);
rss_new=1/beta-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% THIS BLOCK SETS UP AND SOLVES THE NEW VALUE FUNCTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute return function (utility function). We use meshgrid
[kp,k]=meshgrid(K0);
Glow=zeros(knum, knum);
Gmax=k.^alpha*ell^(1-alpha)+(1-delta)*k;
kp(kp>=Gmax | kp<Glow)=NaN; %restrict control variable to feasible set
utility=(k.^alpha*ell^(1-alpha)+(1-delta)*k-kp).^(1-gamma)/(1-gamma);

%initialize the policy and VF arrays (note: we can initialize at V---the
%current value function before the shock)
Vold = V; %this will store real value of old VF
kprimeind=zeros(knum,1);

disp('%%%%')
disp('Second VFI')
disp(' ')

for vfit = 1:vfmaxit
    %form Bellman equation
    RHSMAT = utility + beta*repmat(Vold',knum,1);
   	[V,kprimeind] = max(RHSMAT,[],2);
    absdiff=abs(V-Vold);
    vferr = max(absdiff); %maximum absolute error
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

toc;
disp('%%%%')
disp(' ')

% obtain the new policy functions %
Kpol_new = K0(kprimeind); %real values of capital policy function
Cpol_new = K0.^alpha*ell^(1-alpha) + (1-delta)*K0 - Kpol_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARISSON OF NEW AND OLD POLICY FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

policy_comparisson=figure; 
subplot(1,2,1)
plot(K0,Kpol_new,'r', K0,Kpol,'--r', K0,K0,'k--', 'LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Capital Tomorrow')
title(' Capital Policy')
legend('New Policy', 'Old Policy', '45 degree line','Location','NorthWest');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(1,2,2)
plot(K0,Cpol_new,'k',K0,Cpol,'--k','LineWidth',lwidnum)
xlabel('Capital Today')
ylabel('Consumption');
legend('New Policy', 'Old Policy');
title(' Consumption Policy')
set(gca,'FontSize',fsizenum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS BLOCK COMPUTES THE TRANSITION FROM THE STEADY STATE BEFORE THE SHOCK %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kstart=kss;                     %initial condition---past steady state
[~,kindex]=min(abs(K0-kstart)); %grid point closer to initial condition
kinit=K0(kindex);               %point in grid closer to kstart as initial condition

ktransit=zeros(T,1);    %preallocate the capital transition vector
ktransit(1)=kinit;      %first entry of capital transition is initial condition
for it=2:T
    ktransit(it)=Kpol_new(K0==kinit); %get from the grid the initial value and apply the PF
    kinit=ktransit(it); %update initial capital value
end

%transition for other interesting equilibrium objects
Wtransit=(1-alpha)*(ktransit/ell).^(alpha); %transition for wages
Rtransit=alpha*(ktransit/ell).^(alpha-1); %transition for Rental rate
rtransit=Rtransit-delta; %transition for Rental rate

%includes some points before the shock to capture the pre-shock st st %
transitspan=[-10 0 transitspan];
ktransit=[kss kss ktransit'];
Wtransit=[Wss Wss Wtransit'];
Rtransit=[Rss Rss Rtransit'];
rtransit=[rss rss rtransit'];

transition_new=figure;
subplot(2,2,1)
plot(transitspan,ktransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[kss_new,kss_new],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Capital');
title('Transition for Capital')
legend('Capital','New st st','Location','SouthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,2)
plot(transitspan,Rtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[Rss_new,Rss_new],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Rental Rate');
title('Transition for Rental Rate')
legend('Rental Rate','New st st','Location','NorthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,3)
plot(transitspan,rtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[rss_new,rss_new],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Interest Rate');
title('Transition for Interest Rate')
legend('Interest Rate','New st st','Location','NorthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)

subplot(2,2,4)
plot(transitspan,Wtransit,'b--','LineWidth',lwidnum)
line([transitspan(1) transitspan(end)],[Wss_new,Wss_new],'Color',[1 0 0],'LineStyle','--')
xlabel('Time')
ylabel('Real Wage');
title('Transition for Real Wage')
legend('Real Wage','New st st','Location','SouthEast');
legend boxoff;
set(gca,'FontSize',fsizenum)


