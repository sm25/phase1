phase1
======

ece4784_project_phase1

%%Hodgkin-Huxley
t=0:0.0001:0.1;
gK=0.036;
gNa=0.12;
gL=0.003;
E_K=(-1*0.012);
E_Na=0.115;
E_L=0.0106;
V_rest=(-1*0.07);
C_m=1*(10^(-6));
plot(t,V_rest,'r-');subplot(3,1,1);plot(t,gK,'r-');subplot(3,1,2);plot(t,gNa,'b-');subplot(3,1,3);
% Gating variables
V=0;
alpha_m= 0.1*((25-V)/exp((25-V)/10)-1);
beta_m=4*exp(((-1)*V)/18);
alpha_n=0.01*((10-V)/(exp((10-V)/10)-1));
beta_n=0.0125*exp(-V/80);
alpha_h=0.07*exp(-V/20);
beta_h=1/(exp((30-V)/10)+1);

% m0 = alpha_m/(alpha_m+beta_m);
% n0 = alpha_n/(alpha_n+beta_n);
% h0 = alpha_h/(alpha_h+beta_h);

% Currents

I_Na=(m^3)*gNa*h*(V-E_Na);
I_K=(n^4)*gK*(V-E_K);
I_L=gL*(V-E_L);
I_ion=I-(I_Na)-(I_K)-(I_L);

%Derivatives

frac=0.00005;
V_new=V+(frac*(I_ion/C_m));
m = m0 + frac*(alpha_m*(1-m) - (beta_m*m));
n = n0 + frac*(alpha_n*(1-n) - (beta_n*n));
h = h0 + frac*(alpha_h*(1-h) - (beta_h*h));
