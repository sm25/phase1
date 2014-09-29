%%Hodgkin-Huxley

%Initialize time
frac=0.01;
t=0:frac:100;  


gK=36*10^(-3);              %Initial conductance for K
gNa=120*10^(-3);            %Initial conductance for Na
gL=0.3*10^(-3);             %Initial conductance for the leakage current
E_K=(-12*10^(-3));          %Initial potential for K channels
E_Na=115*10^(-3);           %Initial potential for Na channels
E_L=10.6*10^(-3);           %Initial potential associated with leakage
V_rest=(-70)*10^(-3);       %Resting potential
C_m=1*(10^(-6));            %Membrane capacitance
            
%%Initializing vectors!!

V=zeros([1 length(t)]);
V_m=zeros([1 length(t)]);

%Initialize the gating variables

alpha_m=zeros([1 length(t)]);
alpha_n=zeros([1 length(t)]);
alpha_h=zeros([1 length(t)]);
beta_m=zeros([1 length(t)]);
beta_n=zeros([1 length(t)]);
beta_h=zeros([1 length(t)]);

%Initialize the conductance related variables

m=zeros([1 length(t)]);
n=zeros([1 length(t)]);
h=zeros([1 length(t)]);

%Initialize the various current vectors
I_ion=zeros([1 length(t)]);


%%%%%Trying to provide stumulus
I_inj=zeros([1 length(t)]);
I_inj_1=ones([1 500]); % Specify the duration of the stimulus in terms of the no. of samples 
a=length(I_inj_1);
I_inj(10: 9+a)=((5*10^(-6))*I_inj_1);

%Provide intial conditions!

I_ion(1)=0;

V(1)=V_rest/(10^(-3));
V_m(1)=V_rest;

alpha_m(1)= 0.1*((25-V(1))/(exp((25-V(1))/10)-1));
beta_m(1)=4*(exp((-V(1))/18));
alpha_n(1)=0.01*((10-V(1))/(exp((10-V(1))/10)-1));
beta_n(1)=0.125*(exp(-V(1)/80));
alpha_h(1)=0.07*(exp(-V(1)/20));
beta_h(1)=1/(exp((30-V(1))/10)+1);

%Calculate values for the conductance related variables

m(1) = alpha_m(1)/(alpha_m(1)+beta_m(1));
n(1) = alpha_n(1)/(alpha_n(1)+beta_n(1));
h(1) = alpha_h(1)/(alpha_h(1)+beta_h(1));


    for i=2:(length(t)-1)
 
        V(i)=V_m(i-1)/(10^(-3));    %Scale the voltage to fit the equations, use the previous value of V_m
        
        %Calculate values for the gating variables
        
        alpha_m(i) = 0.1*((25-V(i))/(exp((25-V(i))/10)-1));
        beta_m(i) = 4*(exp(-V(i)/18));
        alpha_n(i) = 0.01*((10-V(i))/(exp((10-V(i))/10)-1));
        beta_n(i) = 0.125*(exp(-V(i)/(80)));
        alpha_h(i) = 0.07*(exp(-V(i)/(20)));
        beta_h(i) = 1/(exp((30-V(i))/10)+1);

        %Update the variables to the next value using Euler's method

        m(i) = m(i-1) + frac*(alpha_m(i)*(1-m(i-1))-beta_m(i)*m(i-1));
        n(i) = n(i-1) + frac*(alpha_n(i)*(1-n(i-1))-beta_n(i)*n(i-1));
        h(i) = h(i-1) + frac*(alpha_h(i)*(1-h(i-1))-beta_h(i)*h(i-1));
        
        % Currents, calculate the current values
        
        I_Na = (m(i)^3)*gNa*h(i)*(V_m(i-1)-E_Na);
        I_K = (n(i)^4)*gK*(V_m(i-1)-E_K);
        I_L = gL*(V_m(i-1)-E_L);
        
        I_ion(i)=I_inj(i)-(I_Na)-(I_K)-(I_L);

        %Update the Voltage to the next value using Euler's method

        V_m(i)=V_m(i-1)+(frac* 10^(-3)*(I_ion(i)/C_m));
        
    end
    
    V_m(i)=V_m(i)-(V_rest);         %Calculate the HH voltage
    plot(t, V_m)                    %The voltage plot
    
%     plot(t,m);                    %the Potassium conductance
%     hold on
%     plot(t,n,'r-')                %the Sodium conductance
    
