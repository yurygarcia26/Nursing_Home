%This material corresponds to the result presented in the article: 
%Effect of Non-Pharmaceutical Interventions to Contain the Spread of SARS-CoV-2 in Nursing Homes
%In which a study is carried out on the impact that the implementation of 
%different non-pharmaceutical measures could have to prevent the spread of COVID-19 within a nursing home.

%A complete description of the model implemented in this code can be found
%in: Protecting residential care facilities from pandemic influenza
%
%Description about the Categories, parameters, etc can be found in:
%Effect of Non-Pharmaceutical Interventions to Contain the Spread of SARS-CoV-2 in Nursing Homes


%TO GET RESULTS
%--------------------------------------------------------------------------
%1. Change the number of simulation in line 167, it can take a while 
%   depending of the number of simulations.
%2. Run this script. It is important to change the name of
%   the file to save at the end of this script according with the Category you
%   are running.
%3. After get the files, save them into a folder named Data to run the scripts
%   to plot results


%IMPORTANT:!!!!
%--------------------------------------------------------------------------
%You have to remove the comment to the scenario you want to study
%eta = parameter to be changed according with the % of asymptomatic
%isolation
% eta = (1-%of isolation)
%--------------------------------------------------------------------------

%%BASELINE
%%Probability of escape monitoring efforts
p_E=1;  pI=1;   pA=1;  
    
%%Time that visitors and staff expend into the facilities
tau_VF = 12;     %Visit of 2h a week
tau_VC =.1428;   %Time in community
tau_SF = 3;      %Staff work 8h for week
tau_SC = 1.5;    %Staff 16h in community

%%Transmission reduction parameter (Since same for all groups-no-index here)
rhoi  = 1;      %Transmission reduction for the use of face mask, etc
rho_R = 1;      %Residents don't use face mask
eta   = 1;      %proportion of asymtomatic isolation
pi_R  = 1;      %Transmission reduction for infected residen isolation
pi_SF = 1;      %Transmission reduction for infected staff isolation  
pi_VF = 1;      %Transmission reduction for infected staff isolation
pi    = 1;      %Transmission reduction in the rest of community 
pi_C  = pi; pi_SC = pi; pi_VC = pi; 

%--------------------------------------------------------------------------
%%CATEGORY 1-2
%--------------------------------------------------------------------------
%%Probability of escape monitoring efforts
% p_E=1;  pI=1;   pA=1;
%     
% %%Time that visitors and staff expend into the facilities
% tau_VF = 12;     %Visit of 2h a week
% tau_VC =.1428;   %Time in community
% tau_SF = 3;      %Staff work 8h for week
% tau_SC = 1.5;    %Staff 16h in community
% 
% %%Transmission reduction parameter (Since same for all groups-no-index here)
% rhoi  = 0.86;   %Transmission reduction for the use of face mask, etc
% rho_R = 1;      %Residents don't use face mask, etc
% eta   = (1-0.5);    %Proportion of asymtomatic isolation
% eta_E = 1;      %Infectiousness of pepople in latency period  
% eta_A = 1;      %Infectiousness of asymptomatic relative to symptomatic 
% pi_R  = 0;      %Transmission reduction for infected residen isolation
% pi_SF = 0;      %Transmission reduction for infected staff isolation  
% pi_VF = 0;      %visitor with symptoms are not allowed to enter
% pi    = 1;      %Transmission reduction in the rest of community (It is not possible to control this outside the facility)
% pi_C  = pi; pi_SC = pi; pi_VC = pi; 

%--------------------------------------------------------------------------
%%CATEGORY 3-4
%%Probability of escape monitoring efforts
%--------------------------------------------------------------------------
% p_E=0.5;  pI=0.1;   pA=1;
%     
% %%Time that visitors and staff expend into the facilities
% tau_VF = 24;     %Visit of 2h a week
% tau_VC =.1428;   %Time in community
% tau_SF = 2;      %Staff work 12h for week
% tau_SC = 2;      %Staff 16h in community
% 
% %%Transmission reduction parameter (Since same for all groups-no-index here)
% rhoi  = 0.86;   %Transmission reduction for the use of face mask, etc
% rho_R = 1;      %Residents don't use face mask, etc
% eta   = (1-0.5);    %Infectiousness of people in latency period
% eta_E = 1;
% eta_A = 1;
% pi_R  = 0;      %Transmission reduction for infected residen isolation
% pi_SF = 0;      %Transmission reduction for infected staff isolation  
% pi_VF = 0;      %visitor with symptoms are not allowed to enter
% pi    = 1;      %Transmission reduction in the rest of community (It is not possible to control this outside the facility)
% pi_C  = pi; pi_SC = pi; pi_VC = pi; 


% %--------------------------------------------------------------------------
% %%CATEGORY 5
% %%Probability of escape monitoring efforts
% %--------------------------------------------------------------------------
% p_E=0.14;  pI=0.1;   pA=1;   
% %%Time that visitors and staff expend into the facilities
% tau_VF = 0;     %No visit
% tau_VC = 0;     %No visit
% tau_SF = 1/4;   %Staff work 4-day-on
% tau_SC = 1/4;   %Staff in community 4 days
% 
% %%Transmission reduction parameter (Since same for all groups-no-index here)
% rhoi  = 0.86;   %Transmission reduction for the use of face mask, etc
% rho_R = 1;      %Residents don't use face mask, etc
% eta   = (1-0.5);    %Infectiousness of people in latency period
% eta_E = 1;
% eta_A = 1;
% pi_R  = 0;      %Transmission reduction for infected residen isolation
% pi_SF = 0;      %Transmission reduction for infected staff isolation  
% pi_VF = 0;      %visitor with symptoms are not allowed to enter
% pi    = 1;      %Transmission reduction in the rest of community (It is not possible to control this outside the facility)
% pi_C  = pi; pi_SC = pi; pi_VC = pi; 
% 


%PARAMETERS
%Beta values

for k=1:11

values = 0.61:0.06:1.21;
beta   = values(k);

%beta in 0.61-1.21 gives R0=2-4 (increments of 0.06 produce increments of 0.2 in R0)
beta_R = beta;
beta_SF= beta; 
beta_VF= beta; 
beta_SC= beta; 
beta_VC= beta; 
beta_C = beta;
    
  
%Recovery rate
gamma_A =0.07142;  %1/14 
gamma_I =0.07142;  %1/14;
   
  
phi     = 0.2;     %1/phi average latency period 
m       = 0.667;   %Proportion of people that moves from E to I 
theta_h = 0.16/4;  %Hospitalization rate
theta_u = 0.32/4;  %UCI rate
mu_h    = 0.15/11; %Hospitalized resident mortality rate
mu_u    = 0.33/11; %UCI resident mortality rate
gamma_h = 1/12;    %Hospitalized resident recovery rate
gamma_u = 1/13;    %UCI resident recovery rate

vectorResI    =[];
vectorResD    =[];
vectorResE    =[];
vectorStaffI  =[];
vectorVisI    =[];
vectorResR    =[];
vectorVisR    =[];
vectorStaffR  =[];
vectorResH    =[];
vectorResU    =[];
vectorResHT   =[];
vectorResUT   =[];
vectorResIT   =[];
vectorResAT   =[];  

%time step
dt=0.001;

%epidemic duration (days)
T=200;

%time vector
ts=0:dt:T;

tic
parfor realizations=1:100
    %Poisson simulation of the ADL and community model

    %initial conditions(resident popn)
    S_R0=200;
    E_R0=0;
    I_R0=0;
    A_R0=0;
    D_R0=0;
    R_R0=0;
    H_R0=0;
    U_R0=0;
    HT_0=0;
    UT_0=0;
    IT_0=0;
    AT_0=0;
    
    
    %initial conditions(visitors in community and facility pops)
    S_VC0=27;
    E_VC0=1;
    A_VC0=0;
    I_VC0=1;
    R_VC0=0;
    S_VF0=13;
    E_VF0=0;
    A_VF0=0;
    I_VF0=0;
    R_VF0=0;
    %initial conditions(staff in community and facility)
    S_SC0=8;
    E_SC0=1;
    A_SC0=0;
    I_SC0=1;
    R_SC0=0;
    S_SF0=75;
    E_SF0=0;
    A_SF0=0;
    I_SF0=0;
    R_SF0=0;
    
    %initial conditions(general community)
    S_C0=50000;
    E_C0=1;
    A_C0=0;
    I_C0=1;
    R_C0=0;
       
    set(0,'RecursionLimit',500)

    %create vectors with variables (Residents)
    S_R=zeros(length(ts),1);
    E_R=zeros(length(ts),1);
    I_R=zeros(length(ts),1);
    A_R=zeros(length(ts),1);
    D_R=zeros(length(ts),1);
    R_R=zeros(length(ts),1);
    H_R=zeros(length(ts),1);
    U_R=zeros(length(ts),1);
    HT =zeros(length(ts),1);
    UT =zeros(length(ts),1);
    IT =zeros(length(ts),1);
    AT =zeros(length(ts),1);
     
  
    %create vectors with variables (Visitors in community and Facility)
    S_VC=zeros(length(ts),1);
    E_VC=zeros(length(ts),1);
    A_VC=zeros(length(ts),1);
    I_VC=zeros(length(ts),1);
    R_VC=zeros(length(ts),1);
    S_VF=zeros(length(ts),1);
    E_VF=zeros(length(ts),1);
    A_VF=zeros(length(ts),1);
    I_VF=zeros(length(ts),1);
    R_VF=zeros(length(ts),1);

    %create vectors with variables (Staff in community and Facility)
    S_SC=zeros(length(ts),1);
    E_SC=zeros(length(ts),1);
    A_SC=zeros(length(ts),1);
    I_SC=zeros(length(ts),1);
    R_SC=zeros(length(ts),1);
    S_SF=zeros(length(ts),1);
    E_SF=zeros(length(ts),1);
    A_SF=zeros(length(ts),1);
    I_SF=zeros(length(ts),1);
    R_SF=zeros(length(ts),1);

    %create vectors with variables (General Community)
    S_C=zeros(length(ts),1);
    E_C=zeros(length(ts),1);
    A_C=zeros(length(ts),1);
    I_C=zeros(length(ts),1);
    R_C=zeros(length(ts),1);
    
    %initialize variables(Residents)
    S_R(1)= S_R0;
    E_R(1)= E_R0;
    I_R(1)= I_R0;
    A_R(1)= A_R0;
    D_R(1)= D_R0;
    R_R(1)= R_R0;
    H_R(1)= H_R0;
    U_R(1)= U_R0;
    HT(1) = HT_0;
    UT(1) = UT_0;
    IT(1) = IT_0;
    AT(1) = AT_0;

    %initialize variables(Visitors:Community and Facility)
    S_VC(1)= S_VC0;
    E_VC(1)= E_VC0;
    A_VC(1)= A_VC0;
    I_VC(1)= I_VC0;
    R_VC(1)= R_VC0;
    
    S_VF(1)= S_VF0;
    E_VF(1)= E_VF0;
    A_VF(1)= A_VF0;
    I_VF(1)= I_VF0;
    R_VF(1)= R_VF0;
    
    %initialize variables(Staff:Community and Facility)
    S_SC(1)=S_SC0;
    E_SC(1)=E_SC0;
    A_SC(1)=A_SC0;
    I_SC(1)=I_SC0;
    R_SC(1)=R_SC0;
   
    S_SF(1)=S_SF0;
    E_SF(1)=E_SF0;
    A_SF(1)=A_SF0;
    I_SF(1)=I_SF0;
    R_SF(1)=R_SF0;
    
    %initialize variables(General Community)
    S_C(1)=S_C0;
    E_C(1)=E_C0;
    A_C(1)=A_C0;
    I_C(1)=I_C0;
    R_C(1)=R_C0;
   
    
    t=0;

    for i=1:1:length(ts)-1
        %ts(i)
          
        N_R=S_R(i)+E_R(i)+A_R(i)+I_R(i)+R_R(i)+H_R(i)+U_R(i);              %tot resident pop
        N_C=S_C(i)+E_C(i)+I_C(i)+A_C(i)+R_C(i);                            %tot general community

        N_VC=S_VC(i)+E_VC(i)+I_VC(i)+A_VC(i)+R_VC(i);                      %tot visitors in comm
        N_VF=S_VF(i)+E_VF(i)+I_VF(i)+A_VF(i)+R_VF(i);                      %tot visitors in facility
        N_SF=S_SF(i)+E_SF(i)+I_SF(i)+A_SF(i)+R_SF(i);                      %tot staff in facility
        N_SC=S_SC(i)+E_SC(i)+I_SC(i)+A_SC(i)+R_SC(i);                      %tot staff in comm

        N_out=N_VC+N_SC +N_C;                                              %tot outside ADL
        N_in =N_R +N_SF +N_VF;                                             %tot in ADL
        N_tot=N_in+N_out;                                                  %tot popn modeled here

         %Forces of infection: In facility, vacc, prophylaxis, in community
        L_in =(beta_R*rho_R*(pi_R*I_R(i)+eta*(A_R(i)+E_R(i)))+beta_VF*rhoi*(pi_VF*I_VF(i)+eta*(A_VF(i)+E_VF(i)))+beta_SF*rhoi*(pi_SF*I_SF(i)+eta*(A_SF(i)+E_SF(i))))/N_in;

        L_out=(beta_C*rhoi*(pi_C*I_C(i)+A_C(i) + E_C(i))+beta_SC*rhoi*(pi_SC*I_SC(i)+A_SC(i)+ E_SC(i))+beta_VC*rhoi*(pi_VC*I_VC(i)+A_VC(i)+E_VC(i)))/N_out;
        
        %Events Pertaining to activity of ADL RESIDENTS
        event1  =poissrnd(dt*L_in*S_R(i))/dt;           %event 1
        event2  =poissrnd(dt*m*phi*E_R(i))/dt;          %event 2
        event3  =poissrnd(dt*(1-m)*phi*E_R(i))/dt;      %event 3
        event4  =poissrnd(dt*gamma_A*A_R(i))/dt;        %event 4
        event5  =poissrnd(dt*gamma_I*I_R(i))/dt;        %event 5
        event6  =poissrnd(dt*theta_h*I_R(i))/dt;        %event 6
        event52 =poissrnd(dt*theta_u*H_R(i))/dt;        %event 52
        event53 =poissrnd(dt*gamma_h*H_R(i))/dt;        %event 53
        event54 =poissrnd(dt*gamma_u*U_R(i))/dt;        %event 54
        event55 =poissrnd(dt*mu_h*H_R(i))/dt;           %event 55
        event56 =poissrnd(dt*mu_u*U_R(i))/dt;           %event 56
       
        
        %Events Pertaining to VISITORS in the Community and Facility
        event7 =poissrnd(dt*tau_VF*S_VF(i))/dt;                      %event 1
        event8 =poissrnd(dt*L_out*S_VC(i))/dt;                       %event 2
        event9 =poissrnd(dt*tau_VC*S_VC(i))/dt;                      %event 3
        event10=poissrnd(dt*tau_VF*E_VF(i))/dt;                      %event 4
        event11=poissrnd(dt*p_E*tau_VC*E_VC(i))/dt;                  %event 5
        event12=poissrnd(dt*m*(1-p_E)*phi*E_VC(i))/dt;               %event 6
        event13=poissrnd(dt*(1-m)*(1-p_E)*phi*E_VC(i))/dt;           %event 7
        event14=poissrnd(dt*tau_VF*A_VF(i))/dt;                      %event 8
        event15=poissrnd(dt*gamma_A*(1-pA)*A_VC(i))/dt;              %event 9
        event16=poissrnd(dt*pA*tau_VC*A_VC(i))/dt;                   %event 10
        event17=poissrnd(dt*tau_VF*I_VF(i))/dt;                      %event 11
        event18=poissrnd(dt*gamma_I*I_VC(i))/dt;                     %event 12
        event19=poissrnd(dt*pI*tau_VC*I_VC(i))/dt;                   %event 13
        event20=poissrnd(dt*tau_VF*R_VF(i))/dt;                      %event 14
        event21=poissrnd(dt*tau_VC*R_VC(i))/dt;                      %event 15
        event22=poissrnd(dt*L_in*S_VF(i))/dt;                        %event 16
        event23=poissrnd(dt*m*phi*E_VF(i))/dt;                       %event 17
        event24=poissrnd(dt*(1-m)*phi*E_VF(i))/dt;                   %event 18
        event25=poissrnd(dt*gamma_A*A_VF(i))/dt;                     %event 19
        event26=poissrnd(dt*gamma_I*I_VF(i))/dt;                     %event 20

        %Events Pertaining to STAFF in the Community and Facility
        event27=poissrnd(dt*tau_SF*S_SF(i))/dt;                       %event 1
        event28=poissrnd(dt*L_out*S_SC(i))/dt;                        %event 2
        event29=poissrnd(dt*tau_SC*S_SC(i))/dt;                       %event 3
        event30=poissrnd(dt*tau_SF*E_SF(i))/dt;                       %event 4
        event31=poissrnd(dt*p_E*tau_SC*E_SC(i))/dt;                   %event 5
        event32=poissrnd(dt*m*(1-p_E)*phi*E_SC(i))/dt;                %event 6
        event33=poissrnd(dt*(1-m)*(1-p_E)*phi*E_SC(i))/dt;            %event 7
        event34=poissrnd(dt*tau_SF*A_SF(i))/dt;                       %event 8
        event35=poissrnd(dt*gamma_A*(1-pA)*A_SC(i))/dt;               %event 9
        event36=poissrnd(dt*pA*tau_SC*A_SC(i))/dt;                    %event 10
        event37=poissrnd(dt*tau_SF*I_SF(i))/dt;                       %event 11
        event38=poissrnd(dt*gamma_I*I_SC(i))/dt;                      %event 12
        event39=poissrnd(dt*pI*tau_SC*I_SC(i))/dt;                    %event 13
        event40=poissrnd(dt*tau_SF*R_SF(i))/dt;                       %event 14
        event41=poissrnd(dt*tau_SC*R_SC(i))/dt;                       %event 15
        event42=poissrnd(dt*L_in*S_SF(i))/dt;                         %event 16
        event43=poissrnd(dt*m*phi*E_SF(i))/dt;                        %event 17
        event44=poissrnd(dt*(1-m)*phi*E_SF(i))/dt;                    %event 18
        event45=poissrnd(dt*gamma_A*A_SF(i))/dt;                      %event 19
        event46=poissrnd(dt*gamma_I*I_SF(i))/dt;                      %event 20

        %Events Pertaining General COMMUNITY
        event47=poissrnd(dt*L_out*S_C(i))/dt;                    %event 1
        event48=poissrnd(dt*m*phi*E_C(i))/dt;                    %event 2
        event49=poissrnd(dt*(1-m)*phi*E_C(i))/dt;                %event 3
        event50=poissrnd(dt*gamma_A*A_C(i))/dt;                  %event 4
        event51=poissrnd(dt*gamma_I*I_C(i))/dt;                  %event 5
      
        %Equations for facility RESIDENTS
        S_R(i+1)=S_R(i)-dt*event1;
        S_R(i+1)=S_R(i+1)*(S_R(i+1)>0);
        E_R(i+1)=E_R(i)+dt*event1-dt*(event2+event3);
        E_R(i+1)=E_R(i+1)*(E_R(i+1)>0);
        I_R(i+1)=I_R(i)+dt*event2-dt*(event5+event6);
        I_R(i+1)=I_R(i+1)*(I_R(i+1)>0);
        A_R(i+1)=A_R(i)+dt*event3-dt*event4;
        A_R(i+1)=A_R(i+1)*(A_R(i+1)>0);
        R_R(i+1)=R_R(i)+dt*(event4+event5+event53+event54);
        R_R(i+1)=R_R(i+1)*(R_R(i+1)>0);
        D_R(i+1)=D_R(i)+dt*event55+dt*event56;
        D_R(i+1)=D_R(i+1)*(D_R(i+1)>0);
        H_R(i+1)=H_R(i)+dt*event6-dt*(event52+event53+event55);
        H_R(i+1)=H_R(i+1)*(H_R(i+1)>0);
        U_R(i+1)=U_R(i)+dt*event52-dt*(event54+event56);
        U_R(i+1)=U_R(i+1)*(U_R(i+1)>0);
        HT(i+1) =HT(i)+ dt*event6;
        HT(i+1) =HT(i+1)*(HT(i+1)>0);
        UT(i+1) =UT(i)+ dt*event52; 
        UT(i+1) =UT(i+1)*(UT(i+1)>0);
        IT(i+1) =IT(i)+ dt*event2;
        IT(i+1) =IT(i+1)*(IT(i+1)>0);
        AT(i+1) =AT(i)+ dt*event3; 
        AT(i+1) =AT(i+1)*(AT(i+1)>0);

        %Equations for VISITORS in Community and Facility
        S_VC(i+1)=S_VC(i)+dt*event7-dt*(event8+event9);
        S_VC(i+1)=S_VC(i+1)*(S_VC(i+1)>0);
        E_VC(i+1)=E_VC(i)+dt*(event10+event8)-dt*(event11+event12+event13);
        E_VC(i+1)=E_VC(i+1)*(E_VC(i+1)>0);
        A_VC(i+1)=A_VC(i)+dt*(event14+event13)-dt*(event15+event16);
        A_VC(i+1)=A_VC(i+1)*(A_VC(i+1)>0);
        I_VC(i+1)=I_VC(i)+dt*(event17+event12)-dt*(event18+event19);
        I_VC(i+1)=I_VC(i+1)*(I_VC(i+1)>0);
        R_VC(i+1)=R_VC(i)+dt*(event20+event18+event15)-dt*event21;
        R_VC(i+1)=R_VC(i+1)*(R_VC(i+1)>0);
        
        S_VF(i+1)=S_VF(i)+dt*event9-dt*(event22+event7);
        S_VF(i+1)=S_VF(i+1)*(S_VF(i+1)>0);
        E_VF(i+1)=E_VF(i)+dt*(event11+event22)-dt*(event23+event24+event10);
        E_VF(i+1)=E_VF(i+1)*(E_VF(i+1)>0);
        A_VF(i+1)=A_VF(i)+dt*(event24+event16)-dt*(event25+event14);
        A_VF(i+1)=A_VF(i+1)*(A_VF(i+1)>0);
        I_VF(i+1)=I_VF(i)+dt*(event23+event19)-dt*(event26+event17);
        I_VF(i+1)=I_VF(i+1)*(I_VF(i+1)>0);
        R_VF(i+1)=R_VF(i)+dt*(event21+event25+event26)-dt*event20;
        R_VF(i+1)=R_VF(i+1)*(R_VF(i+1)>0);
        
        %Equations for STAFF in Community and Facility
        S_SC(i+1)=S_SC(i)+dt*event27-dt*(event28+event29);
        S_SC(i+1)=S_SC(i+1)*(S_SC(i+1)>0);
        E_SC(i+1)=E_SC(i)+dt*(event30+event28)-dt*(event31+event32+event33);
        E_SC(i+1)=E_SC(i+1)*(E_SC(i+1)>0);
        A_SC(i+1)=A_SC(i)+dt*(event34+event33)-dt*(event35+event36);
        A_SC(i+1)=A_SC(i+1)*(A_SC(i+1)>0);
        I_SC(i+1)=I_SC(i)+dt*(event37+event32)-dt*(event38+event39);
        I_SC(i+1)=I_SC(i+1)*(I_SC(i+1)>0);
        R_SC(i+1)=R_SC(i)+dt*(event40+event35+event38)-dt*event41;
        R_SC(i+1)=R_SC(i+1)*(R_SC(i+1)>0);
        
        S_SF(i+1)=S_SF(i)+dt*event29-dt*(event42+event27);
        S_SF(i+1)=S_SF(i+1)*(S_SF(i+1)>0);
        E_SF(i+1)=E_SF(i)+dt*(event31+event42)-dt*(event43+event44+event30);
        E_SF(i+1)=E_SF(i+1)*(E_SF(i+1)>0);
        A_SF(i+1)=A_SF(i)+dt*(event44+event36)-dt*(event45+event34);
        A_SF(i+1)=A_SF(i+1)*(A_SF(i+1)>0);
        I_SF(i+1)=I_SF(i)+dt*(event43+event39)-dt*(event46+event37);
        I_SF(i+1)=I_SF(i+1)*(I_SF(i+1)>0);
        R_SF(i+1)=R_SF(i)+dt*(event41+event45+event46)-dt*event40;
        R_SF(i+1)=R_SF(i+1)*(R_SF(i+1)>0);

        %Equations for the general community
        S_C(i+1)=S_C(i)-dt*event47;
        S_C(i+1)=S_C(i+1)*(S_C(i+1)>0);
        E_C(i+1)=E_C(i)+dt*event47-dt*(event48+event49);
        E_C(i+1)=E_C(i+1)*(E_C(i+1)>0);
        A_C(i+1)=A_C(i)+dt*event49-dt*event50;
        A_C(i+1)=A_C(i+1)*(A_C(i+1)>0);
        I_C(i+1)=I_C(i)+dt*event48-dt*event51;
        I_C(i+1)=I_C(i+1)*(I_C(i+1)>0);
        R_C(i+1)=R_C(i)+dt*(event51+event50);
        R_C(i+1)=R_C(i+1)*(R_C(i+1)>0);
      
    end
    vectorResI   =[vectorResI I_R];       % save I curves (one per column)
    vectorStaffI =[vectorStaffI I_SF];    % save I curves (one per column)
    vectorVisI   =[vectorVisI I_VF];      % save I curves (one per column)
    vectorResD   =[vectorResD D_R];       % save D curves (one per column)
    vectorResR   =[vectorResR R_R];       % save H curves (one per column)
    vectorVisR   =[vectorVisR R_VF];      % save U curves (one per column)
    vectorStaffR =[vectorStaffR R_SF];    % save I curves (one per column)
    vectorResH   =[vectorResH H_R];       % save I curves (one per column)
    vectorResU   =[vectorResU U_R];       % save I curves (one per column)
    vectorResHT  =[vectorResHT HT];       % save I curves (one per column)
    vectorResUT  =[vectorResUT UT];       % save I curves (one per column)
    vectorResIT  =[vectorResIT IT];       % save I curves (one per column)
    vectorResAT  =[vectorResAT AT];       % save I curves (one per column)



end
%Change name according with Category you are running
save(['Category0_' num2str(k) '.mat'], 'vectorResI','vectorResR','vectorResH','vectorResU','vectorResD','vectorResIT','vectorResHT','vectorResUT','vectorStaffR','vectorStaffI','vectorVisR')

  
end
toc



