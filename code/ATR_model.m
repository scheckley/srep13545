
function dydt = ATR_model(t,y,parameters)
%%%%%%%%%%%%%
% variables %
%%%%%%%%%%%%%
G1 = y(1);
G2 = y(2);
S = y(3);
apoptosis = y(4);
nDD = y(5);
G1_IR = y(6);
G2_IR = y(7);
S_IR = y(8);
nDD_IR = y(9);

ATR = y(10);
ATRi = y(11);
ATR_ATRi = y(12);

GUT = y(13);
CEN = y(14);
PER = y(15);

apop_IR = y(16);

%%%%%%%%%%%%%%
% parameters %
%%%%%%%%%%%%%%
%in vitro
h = parameters(1);
k1 = parameters(2);
k3 = parameters(3);
k4 = parameters(4);
k5 = parameters(5);
k6 = parameters(6);
k7 = parameters(7);
ki = parameters(8);
z_in_vitro = parameters(9);
z_in_vivo = parameters(10);
kon = parameters(11);
koff = parameters(12);
kscale = parameters(13);
fu = parameters(14);
k_IR_repair = parameters(15);
k_IR_apop = parameters(16);
%in vivo (PK)
v1 = parameters(17);
v2 = parameters(18);
q = parameters(19);
cl = parameters(20);
ka = parameters(21);
%model switches
pk_factor = parameters(22);
in_vitro = parameters(23);
washout = parameters(24);
IR_dose = parameters(25);

if in_vitro < 1
    Cp = (CEN/v1);
    y(11) = Cp;
end

%washout event at t=16hrs. Toggle with parameters(24)
if washout == 1
    if t < 16.0
        y(11)=y(11);
    else
        y(11)=0;
    end
end

%%%%%%%%%%%%%%%%%
%   Equations   %
%%%%%%%%%%%%%%%%%
Reaction_flux(16,1)=0;
Reaction_flux(1)=(k1*k3*y(1));
Reaction_flux(2)=k1*(1-k3)*y(1);
Reaction_flux(3)=k5*y(5);
Reaction_flux(4)=(k4*y(5))*(1-(1/(1+(y(10)/ki)^h)));
Reaction_flux(5)=k7*y(3);
Reaction_flux(6)=k6*y(2);
Reaction_flux(7)=k_IR_repair*y(6);
Reaction_flux(8)=k_IR_repair*y(9);
Reaction_flux(9)=k_IR_repair*y(8);
Reaction_flux(10)=k_IR_repair*y(7);
Reaction_flux(11)=k_IR_apop*y(6);
Reaction_flux(12)=k_IR_apop*y(9);
Reaction_flux(13)=k_IR_apop*y(8);
Reaction_flux(14)=k_IR_apop*y(7);
%switch for in vitro or in vivo scaling on ATRi binding kinetics
if in_vitro == 1
    Reaction_flux(15)=kon*y(10)*y(11); %Rxn_ATR_binding
else
    Reaction_flux(15)=kon*y(10)*kscale*y(11); %Rxn_ATR_binding
end
Reaction_flux(16)=koff*y(12); %Rxn_ATR_unbinding

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      ODE - in vitro     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1 = (2*Reaction_flux(6))-Reaction_flux(1)-Reaction_flux(2)+(2*Reaction_flux(10));
S = Reaction_flux(2)+Reaction_flux(4)-Reaction_flux(5)+(1-k3)*Reaction_flux(7);
nDD = Reaction_flux(1)-Reaction_flux(3)-Reaction_flux(4)+k3*Reaction_flux(7);
G2 = Reaction_flux(5)-Reaction_flux(6)+Reaction_flux(9)+Reaction_flux(8);
apoptosis = Reaction_flux(3);
ATR_ATRi = Reaction_flux(15)-Reaction_flux(16);
ATR = Reaction_flux(16)-Reaction_flux(15);
apop_IR = Reaction_flux(11)+Reaction_flux(12)+Reaction_flux(13)+Reaction_flux(14);
    
% IR damage species
if IR_dose == 1
    G1_IR = -Reaction_flux(7)-Reaction_flux(11);
    nDD_IR = -Reaction_flux(8)-Reaction_flux(12);
    S_IR = -Reaction_flux(9)-Reaction_flux(13);
    G2_IR = -Reaction_flux(10)-Reaction_flux(14);
else
    G1_IR = 0;
    nDD_IR = 0;
    S_IR = 0;
    G2_IR = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      ODE - in vivo      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if in_vitro == 1 % in vitro: no PD d[ATRi]/dt
    GUT=0;
    CEN=0;
    PER=0;
    Cp=0;
    ATRi=0;
else % in vivo d[ATRi]/dt calculated by PK model
    GUT = -ka*y(13);
    CEN = ka*y(13)-(q+cl)*(y(14)/v1)+q*(y(15)/v2);
    PER = q*(y(14)/v1-y(15)/v2);
end

dydt = [G1; G2; S; apoptosis; nDD; G1_IR; G2_IR; S_IR; nDD_IR; ATR; ATRi; ATR_ATRi; GUT; CEN; PER; apop_IR];

end
