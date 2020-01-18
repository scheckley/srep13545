function [t,y] = ATR_simulation()

%in vitro model parameters
parameters(25,1)=0;
parameters(1) = 2.60217; %h
parameters(2) = 0.0811185; %k1
parameters(3) = 0.757926; %k2
parameters(4) = 1766.78; %k3
parameters(5) = 0.130153; %k4
parameters(6) = 0.110204; %k5
parameters(7) = 0.0903928; %k6
parameters(8) = 0.153196; %ki
parameters(9) = 2.95487; % z_in_vitro
parameters(10) = 0.906667; %z_in_vivo
parameters(11) = 0.117798; %kon
parameters(12) = 0.0290027; %koff
parameters(13) = 0.32909; %kscale
parameters(14) = 0.288711; %fu
parameters(15) = 0.018; %k_IR_repair
parameters(16) = 0.039; %k_IR_apoptosis
parameters(17) = 0.101; %v1
parameters(18) = 9.17; %v2
parameters(19) = 0.880; %q
parameters(20) = 0.603; %cl
parameters(21) = 0.592; %ka
parameters(22) = 0.54; %pk_factor
parameters(23) = 1; %in vitro toggle (in vitro=1, in vivo=0)
parameters(24) = 1; %washout toggle (1=yes, 0=no)
parameters(25) = 0; % IR dose at t=0hr (1=yes, 0=no)

timespan = [0 72];

%%%%%%%%%%%%%%%%%%%%%%
% initial conditions %
%%%%%%%%%%%%%%%%%%%%%%
%set up the steady state cell fractions
cell_cycle_frac_G1=0.47;
cell_cycle_frac_nDD=0.0184;
cell_cycle_frac_Stotal=0.32;
cell_cycle_frac_G2=0.21;
cell_cycle_frac_S=(cell_cycle_frac_Stotal-cell_cycle_frac_nDD);

%set up initial cell count
G1 = cell_cycle_frac_G1*1100;
G2 = cell_cycle_frac_G2*1100;
S = cell_cycle_frac_S*1100;
apoptosis = 0.0;
nDD = cell_cycle_frac_nDD*1100;

G1_IR = parameters(25)*cell_cycle_frac_G1*1100;
G2_IR = parameters(25)*cell_cycle_frac_G2*1100;
S_IR = parameters(25)*cell_cycle_frac_S*1100;
nDD_IR = parameters(25)*cell_cycle_frac_nDD*1100;

ATR = 0.01;
ATR_ATRi = 0;
apop_IR = 0.0;

%[drug dose in vitro] = 30, 10, 3, 1, 0.3, 0
%[drug dose in vivo] = 75, 37.5, 25, 10, 1.65
if parameters(23) == 0
    GUT = 75*(1000/429);
    ATRi = GUT;
else
    ATRi = 30.0;
    GUT = ATRi;
end

CEN=0;
PER=0;

initial_conditions=[G1; G2; S; apoptosis; nDD; G1_IR; G2_IR; S_IR; nDD_IR; ATR; ATRi; ATR_ATRi; GUT; CEN; PER; apop_IR];

%% single trajectory
[t,y] = ode23s(@(timespan, initial_conditions) ATR_model(timespan, initial_conditions, parameters), timespan, initial_conditions);
%yH2AX calculation depending on vitro or in vivo
if parameters(23)==1 %in vitro scaling
    yH2AXpan = (((y(:,5)+y(:,9)) ./ (y(:,5)+y(:,3)+y(:,1)+y(:,2)+y(:,9)+y(:,8)+y(:,6)+y(:,7)))*parameters(9))*100;
    yH2AXfoci = (((y(:,6)+y(:,9)+y(:,8)+y(:,7)) ./ (y(:,5)+y(:,3)+y(:,1)+y(:,2)+y(:,9)+y(:,8)+y(:,6)+y(:,7))))*100;
else %in vivo scaling
    yH2AXpan = (((y(:,5)+y(:,9)) ./ (y(:,5)+y(:,3)+y(:,1)+y(:,2)+y(:,9)+y(:,8)+y(:,6)+y(:,7)))*parameters(10))*100;
    yH2AXfoci = (((y(:,6)+y(:,9)+y(:,8)+y(:,7)) ./ (y(:,5)+y(:,3)+y(:,1)+y(:,2)+y(:,9)+y(:,8)+y(:,6)+y(:,7))))*100;
end

%variables for plotting
totalcells = (y(:,5)+y(:,3)+y(:,1)+y(:,2)+y(:,9)+y(:,8)+y(:,6)+y(:,7));
yH2AXtotal = (yH2AXpan(:,1)+yH2AXfoci(:,1));

fractionG1=100*(y(:,1)+y(:,6))./totalcells;
fractionnDD=100*(y(:,5)+y(:,9))./totalcells;
fractionS=100*(y(:,3)+y(:,8))./totalcells;
fractionG2=100*(y(:,2)+y(:,7))./totalcells;

Stotal=(fractionnDD+fractionS);
totalapoptosis=(y(:,4)+y(:,16));
endtotal = (totalcells(end,1));
endpoints=[y(end,1)/endtotal,y(end,2)/endtotal,y(end,3)/endtotal,y(end,5)/endtotal];
names={'G1','S','G2','nDD'};
G1track=(y(:,1)./totalcells(:,1));
G2track=(y(:,2)./totalcells(:,1));
Strack=(y(:,3)./totalcells(:,1));
nDDtrack=(y(:,5)./totalcells(:,1));

%model output
figure(1)
subplot(2,2,1)
hold on
plot(t,yH2AXtotal,'LineWidth',1.5)
box off; grid on;
title('yH2AX','fontsize',12)
xlabel('time (h)','fontsize',10)
ylabel('yH2AX+ve (%)','fontsize',10)

subplot(2,2,2)
hold on
plot(t,totalcells,'LineWidth',1.5)
box off; grid on;
title('total cells','fontsize',12)
xlabel('time (h)','fontsize',10)
ylabel('cell number (#)','fontsize',10)

subplot(2,2,3)
bar(endpoints)
box off;
title('% cell cycle phase distribution at t=end','fontsize',12)
ylabel('percentage of total population (%)','fontsize',10)
set(gca,'xticklabel',names)

subplot(2,2,4)
plot(t,G1track,t,Strack,t,G2track,t,nDDtrack,'LineWidth',1.5)
title('cell cycle phase distribution during time course','fontsize',12)
ylabel('percentage of total population (%)','fontsize',10)
xlabel('time (hrs)','fontsize',10)
legend('G1','S','G2/M','nDD')
box off; grid on

end
