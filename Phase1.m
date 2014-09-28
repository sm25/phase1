%%Hodgkin-Huxley
t=0:0.00025:0.1;    %Define time
gK=0.036;           %Initial conductance for K
gNa=0.12;           %Initial conductance for Na
gL=0.003;           %Initial conductance for the leakage current
E_K=(-1*0.012);     %Initial potential for K channels
E_Na=0.115;         %Initial potential for Na channels
E_L=0.0106;         %Initial potential associated with leakage
V_rest=(-1*0.07);   %Resting potential
C_m=1*(10^(-6));    %Membrane capacitance
frac=0.00025;       %factor by which time is incrementing

%%Initializing vectors!!

%Initialize V with the first value as V_rest
V_new=zeros([1 length(t)]);     
V_i=zeros([1 (length(t)-1)]);
V=[V_rest V_i];
%Initialize the gating variables
alpha_m=zeros([1 length(t)]);
alpha_n=zeros([1 length(t)]);
alpha_h=zeros([1 length(t)]);
beta_m=zeros([1 length(t)]);
beta_n=zeros([1 length(t)]);
beta_h=zeros([1 length(t)]);
%Initialize the conductance related variables
m0=zeros([1 length(t)]);
n0=zeros([1 length(t)]);
h0=zeros([1 length(t)]);
m=zeros([1 length(t)]);
n=zeros([1 length(t)]);
h=zeros([1 length(t)]);
%Initialize the various current vectors
I_Na=zeros([1 length(t)]);
I_K=zeros([1 length(t)]);
I_L=zeros([1 length(t)]);
I_ion=zeros([1 length(t)]);


%%%%%Trying to provide stumulus
I_inj=zeros([1 length(t)]);
% % I_inj_1=ones([1 20]); % Specify the duration of the stimulus in terms of the no. of samples 
% % a=length(I_inj_1);
% % I_inj(1, 10:(9+a))=((5*10^(-6))*I_inj_1);


for i=1:(length(t)-1)
    
%Calculate values for the gating variables
alpha_m(i)= 0.1*((25-V(i))/exp((25-V(i))/10)-1);
beta_m(i)=4.*exp(((-1)*V(i))/18);
alpha_n(i)=0.01*((10-V(i))/(exp((10-V(i))/10)-1));
beta_n(i)=0.0125*exp(-V(i)/80);
alpha_h(i)=0.07*exp(-V(i)/20);
beta_h(i)=1./(exp((30-V(i))/10)+1);


%Calculate values for the conductance related variables

m(1) = alpha_m(1)./(alpha_m(1)+beta_m(1));
n(1) = alpha_n(1)./(alpha_n(1)+beta_n(1));
h(1) = alpha_h(1)./(alpha_h(1)+beta_h(1));


%Update the variables to the next value using Euler's method

m(i) = m(i) + frac.*(alpha_m(i).*(1-m(i)) - (beta_m(i).*m(i)));
n(i) = n(i) + frac.*(alpha_n(i).*(1-n(i)) - (beta_n(i).*n(i)));
h(i) = h(i) + frac.*(alpha_h(i).*(1-h(i)) - (beta_h(i).*h(i)));

% Currents, calculate the current values

I_Na(i)=(m(i).^3)*gNa.*h(i).*(V(i)-E_Na);
I_K(i)=(n(i).^4)*gK.*(V(i)-E_K);
I_L(i)=gL.*(V(i)-E_L);
I_ion(i)=I_inj(i)-(I_Na(i))-(I_K(i))-(I_L(i));

%Update the Voltage to the next value using Euler's method

V_new(i)=V(i)+(frac.*(I_ion(i)/C_m));
V(i+1)=V_new(i);
end

% plot(t,V_rest,'r-');subplot(3,1,1);plot(t,gK,'r-');subplot(3,1,2);plot(t,gNa,'b-');subplot(3,1,3);
% plot(t,V(1:length(t)),'r-');
% plot(t,n,'b-');
% plot(t,m,'g-');

