% Element testing based on a stress-dependent elastoplastic damage model
% author: Lang Liu (lliu@ualberta)
% (Sign convention: COMPRESSION = -; TENSION = +)

clear;
addpath('./Fcns');
DataFolder='.\Data\'; 
%% MODEL PARAMETERS

% AniSD Model
fcn_sig=@fcn_DPIsodwAniDM;
name_para={'E1','\nu','\kappa','c_0','c_1','c_2','b','B'};

% calibrated model parameters for three multi-staged triaxial tests
% reported by Graesle and Plischke (2011):

E1=1200;nu=0.2;kappa=16;c0=1E-3;dmax=0.5;c1=1.5e-3;b=0.5;B=1.4; % test 09026
% E1=800;nu=0.2;kappa=16;c0=1E-3;dmax=0.5;c1=1.5e-3;b=0.5;B=1.2; % test 09027
% E1=900;nu=0.2;kappa=16;c0=1E-3;dmax=0.5;c1=1.5e-3;b=0.5;B=1; % test 09054

val_para={E1,nu,kappa,c0,dmax,c1,b,B};

%% TEST DATA and BOUNDARY VALUES

% _______ TN2010-86 (Grasle & Plischke, 2011)____________
DATA='GP2011_Triax_Full'; 
DATA=[DataFolder,DATA];
load(DATA);

% Multi-stage Test Data 09026 
eps_A=Data_GP_1MPa(:,1);
sig_A=Data_GP_1MPa(:,2);
deviatoricS=Data_GP_1MPa(:,3);
epsA_offset=0.5604e-3;% initial offset to zero start of the data
data_end=25954;eps_A(data_end:end)=[];deviatoricS(data_end:end)=[];
% unload starting and ending strains:
epsA_max(1)=-0.03008;epsA_min(1)=-0.028;SigC(1)=-1; 
epsA_max(2)=-0.05185;epsA_min(2)=-0.0486;SigC(2)=-3; 
epsA_max(3)=-0.0646;epsA_min(3)=-0.0606;SigC(3)=-6; 
epsA_max(4)=-0.07674;epsA_min(4)=-0.0718;SigC(4)=-10; 
epsA_max(5)=-0.0881;epsA_min(5)=-0.0826;SigC(5)=-15; 
epsA_max(6)=-0.09907;epsA_min(6)=-0.0931;SigC(6)=-20; 

% Multi-stage Test Data 09027
% eps_A=Data_GP_3MPa(:,1);
% sig_A=Data_GP_3MPa(:,2);
% deviatoricS=Data_GP_3MPa(:,3);
% epsA_offset=0.4095e-3+eps_A(1,1);
% data_end=18297;eps_A(data_end:end)=[];deviatoricS(data_end:end)=[];
% SigC(1)=-3; epsA_max(1)=-2.71764e-2;epsA_min(1)=-2.27347e-2;
% SigC(2)=-1; epsA_max(2)=-3.26949e-2;epsA_min(2)=-3.03201e-2;
% SigC(3)=-6; epsA_max(3)=-5.72776e-2;epsA_min(3)=-5.2535e-2;
% SigC(4)=-10; epsA_max(4)=-7.51047e-2;epsA_min(4)=-6.9646e-2;
% SigC(5)=-15; epsA_max(5)=-10.0808e-2;epsA_min(5)=-9.47718e-2;
% SigC(6)=-20; epsA_max(6)=-12.9351e-2;epsA_min(6)=-12.084e-2;

% Multi-stage Test Data 09054
% eps_A=Data_GP_6MPa(:,1);
% sig_A=Data_GP_6MPa(:,2);
% deviatoricS=Data_GP_6MPa(:,3);
% epsA_offset=2.341e-3+eps_A(1,1);
% data_end=18325;eps_A(data_end:end)=[];deviatoricS(data_end:end)=[];
% SigC(1)=-6; epsA_max(1)=-4.8375e-2;epsA_min(1)=-4.4216e-2;
% SigC(2)=-1; epsA_max(2)=-6.955e-2;epsA_min(2)=-6.7054e-2;
% SigC(3)=-3; epsA_max(3)=-8.2013e-2;epsA_min(3)=-7.956e-2;
% SigC(4)=-10; epsA_max(4)=-10.009e-2;epsA_min(4)=-9.6432e-2;
% SigC(5)=-15; epsA_max(5)=-12.247e-2;epsA_min(5)=-11.796e-2;
% SigC(6)=-20; epsA_max(6)=-16.264e-2;epsA_min(6)=-15.72e-2;

%% INITIALIZATION
  
% strain/stress step size
depsA=-1e-6;% axial strain increment 
depsA_un=1e-5;% axial strain decrement
dsigC=-0.05; % confining stress increment 

d=0; 
alpha(1)=0; 
alpha(2)=0;
alpha(3)=0;
epsA=0;
sigC=0;
q=0;
deps=zeros(3,3);
eps=zeros(3,3);
eps_e=zeros(3,3);
sigeff=zeros(3,3);
sigeff(1,1)=sigC;
sigeff(2,2)=sigC;
sigeff(3,3)=sigC;

% initialize other state variables
i=1; % computation step
variables={deps,sigeff,eps,eps_e,d,alpha};
Output=fcn_sig(variables{:},val_para);
sig=Output{1};ddsdde=Output{2};C66=Output{3};variables=Output(4:end);

% matrices to store results
SIG(:,i)=[sig(1,1);sig(2,2);sig(3,3)];
SIGEFF(:,i)=[sigeff(1,1);sigeff(2,2);sigeff(3,3)];
EPS_E(:,i)=[eps_e(1,1);eps_e(2,2);eps_e(3,3)];
EPS(:,i)=[eps(1,1);eps(2,2);eps(3,3)];
D(i)=d;

ALPHA(:,i)=alpha;
Epara=Cal_OT(C66);
% calculate engineering elastic parameters
E11(i)=Epara(1,1);E33(i)=Epara(1,3);NU12(i)=Epara(2,1);NU21(i)=Epara(2,2);
NU13(i)=Epara(2,3);NU31(i)=Epara(2,4);
G12(i)=Epara(1,4);G13(i)=Epara(1,5);G23(i)=Epara(1,6);

%% Triaxial Loading

% 1. Isotropic Compression
step0=SigC/dsigC;
for i=length(D)+1:length(D)+step0
    
    sigC=dsigC+sigC;
    prescribed(:,1)=[1;1;1;0;0;0]; % prescribed boundary flags (strain: 0; stress: 1);
    prescribed(:,2)=[0;0;0;0;0;0]; % prescribed boundary strain 
    prescribed(:,3)=[sigC;sigC;sigC;0;0;0]; % prescribed boundary stress
   
    % solve for the strain increment that satisfies prescribed boundaries using
    % Newton Raphson method
    deps=NR(fcn_sig,prescribed,val_para,variables);
    variables{1}=deps;
    % update stress, strain and state variables
    Output=fcn_sig(variables{:},val_para);
    sig=Output{1};ddsdde=Output{2};C66=Output{3};variables=Output(4:end);
    eps=variables{3};sigeff=variables{2};d=variables{5};alpha=variables{6};eps_e=variables{4}; 
    q=-sig(3,3)+sig(1,1);
    
    % store results
    SIG(:,i)=[sig(1,1);sig(2,2);sig(3,3)];
    SIGEFF(:,i)=[sigeff(1,1);sigeff(2,2);sigeff(3,3)];
    EPS_E(:,i)=[eps_e(1,1);eps_e(2,2);eps_e(3,3)];
    EPS(:,i)=[eps(1,1);eps(2,2);eps(3,3)];
    D(i)=d;
    ALPHA(:,i)=alpha;
    Q(i)=q;
    
    % calculate engineering elastic parameters
    Epara=Cal_OT(C66);
    E11(i)=Epara(1,1);E33(i)=Epara(1,3);NU12(i)=Epara(2,1);NU21(i)=Epara(2,2);
    NU13(i)=Epara(2,3);NU31(i)=Epara(2,4);
    G12(i)=Epara(1,4);G13(i)=Epara(1,5);G23(i)=Epara(1,6);

    
end
eps_ini=EPS(:,end); % initial strain after 1st stage compression
sig_ini=SIG(:,end); % initial stress after 1st stage compression


% 2. Axial Loading/Unloading
for istep=1:length(epsA_max)

    if istep==1 % first loading to failure
        epsA_inc=epsA_max(istep);
        sigC_inc=0;
        dsigC_inc=dsigC;
    else % for testing stages after failure
        epsA_inc=epsA_max(istep)-epsA_min(istep-1);
        sigC_inc=SigC(istep)-SigC(istep-1);
        dsigC_inc=-sign(sigC_inc-1e-5)*dsigC;
    end
    
    % Increase the confining pressure to a new level
    step=sigC_inc/dsigC_inc;
    deps=zeros(3); % set the inital trial strain increment to zero 
    
    for i=length(D)+1:length(D)+step
        sigC=sigC+dsigC_inc;
        prescribed(:,1)=[1;1;0;0;0;0]; % prescribed boundary flags (strain: 0; stress: 1);
        prescribed(:,2)=[0;0;eps(3,3);0;0;0]; % prescribed boundary strain
        prescribed(:,3)=[sigC;sigC;0;0;0;0]; 

        % solve for the strain increment that satisfies prescribed boundaries using
        % Newton Raphson method
        deps=NR(fcn_sig,prescribed,val_para,variables);
        variables{1}=deps;
        % update stress, strain and state variables
        Output=fcn_sig(variables{:},val_para);
        sig=Output{1};ddsdde=Output{2};C66=Output{3};variables=Output(4:end);
        eps=variables{3};sigeff=variables{2};d=variables{5};alpha=variables{6};eps_e=variables{4}; 
        q=-sig(3,3)+sig(1,1);
        
        % store results
        SIG(:,i)=[sig(1,1);sig(2,2);sig(3,3)];
        SIGEFF(:,i)=[sigeff(1,1);sigeff(2,2);sigeff(3,3)];
        EPS_E(:,i)=[eps_e(1,1);eps_e(2,2);eps_e(3,3)];
        EPS(:,i)=[eps(1,1);eps(2,2);eps(3,3)];
        D(i)=d;
        ALPHA(:,i)=alpha;
        Q(i)=q;
        
        % calculate engineering elastic parameters
        Epara=Cal_OT(C66);
        E11(i)=Epara(1,1);E33(i)=Epara(1,3);NU12(i)=Epara(2,1);NU21(i)=Epara(2,2);
        NU13(i)=Epara(2,3);NU31(i)=Epara(2,4);
        G12(i)=Epara(1,4);G13(i)=Epara(1,5);G23(i)=Epara(1,6);

    end
    
    % Axial loading

    step=epsA_inc/depsA; % strain controlled
    for i=length(D)+1:length(D)+step
         
        %_____________TRAXIAL COMPRESSION_______________
        epsA=eps(3,3)+depsA;
        prescribed(:,1)=[1;1;0;0;0;0]; % prescribed boundary flags (strain: 0; stress: 1);
        prescribed(:,2)=[0;0;epsA;0;0;0]; % prescribed boundary strain
        prescribed(:,3)=[sigC;sigC;0;0;0;0]; % prescribed boundary stress

        % solve for the strain increment that satisfies prescribed boundaries using
        % Newton Raphson method
        deps=NR(fcn_sig,prescribed,val_para,variables);
        variables{1}=deps;
        % update stress, strain and state variables
        Output=fcn_sig(variables{:},val_para);
        sig=Output{1};ddsdde=Output{2};C66=Output{3};variables=Output(4:end);
        eps=variables{3};sigeff=variables{2};d=variables{5};alpha=variables{6};eps_e=variables{4}; 
        q=-sig(3,3)+sig(1,1);
        
        % store results
        SIG(:,i)=[sig(1,1);sig(2,2);sig(3,3)];
        SIGEFF(:,i)=[sigeff(1,1);sigeff(2,2);sigeff(3,3)];
        EPS_E(:,i)=[eps_e(1,1);eps_e(2,2);eps_e(3,3)];
        EPS(:,i)=[eps(1,1);eps(2,2);eps(3,3)];
        D(i)=d;
        ALPHA(:,i)=alpha;
        Q(i)=q;
        
        % calculate engineering elastic parameters
        Epara=Cal_OT(C66);
        E11(i)=Epara(1,1);E33(i)=Epara(1,3);NU12(i)=Epara(2,1);NU21(i)=Epara(2,2);
        NU13(i)=Epara(2,3);NU31(i)=Epara(2,4);
        G12(i)=Epara(1,4);G13(i)=Epara(1,5);G23(i)=Epara(1,6);
    end
    
    % Axial unloading
    epsA_inc=epsA_min(istep)-epsA_max(istep);
    step=epsA_inc/depsA_un; 
    for i=length(D)+1:length(D)+step
 
        epsA=eps(3,3)+depsA_un;
        prescribed(:,1)=[1;1;0;0;0;0]; % prescribed boundary flags (strain: 0; stress: 1);
        prescribed(:,2)=[0;0;epsA;0;0;0]; % prescribed boundary strain
        prescribed(:,3)=[sigC;sigC;0;0;0;0]; % prescribed boundary stress

        % solve for the strain increment that satisfies prescribed boundaries using
        % Newton Raphson method
        deps=NR(fcn_sig,prescribed,val_para,variables);
        variables{1}=deps;
        % update stress, strain and state variables
        Output=fcn_sig(variables{:},val_para);
        sig=Output{1};ddsdde=Output{2};C66=Output{3};variables=Output(4:end);
        eps=variables{3};sigeff=variables{2};d=variables{5};alpha=variables{6};eps_e=variables{4}; 
        q=-sig(3,3)+sig(1,1);
        
        % store results
        SIG(:,i)=[sig(1,1);sig(2,2);sig(3,3)];
        SIGEFF(:,i)=[sigeff(1,1);sigeff(2,2);sigeff(3,3)];
        EPS_E(:,i)=[eps_e(1,1);eps_e(2,2);eps_e(3,3)];
        EPS(:,i)=[eps(1,1);eps(2,2);eps(3,3)];
        D(i)=d;
        ALPHA(:,i)=alpha;
        Q(i)=q;
        
        % calculate engineering elastic parameters
        Epara=Cal_OT(C66);
        E11(i)=Epara(1,1);E33(i)=Epara(1,3);NU12(i)=Epara(2,1);NU21(i)=Epara(2,2);
        NU13(i)=Epara(2,3);NU31(i)=Epara(2,4);
        G12(i)=Epara(1,4);G13(i)=Epara(1,5);G23(i)=Epara(1,6);

    end
end

% Determine net stress/strain change after the stage of initial isotropic compression
EPS=EPS-eps_ini;
SIG=SIG-sig_ini;
SIGEFF=SIGEFF-sig_ini;

%% Plotting
figure; hold on

% ____Stress-Strain Curve_____ 
% *********Plot Test Data******************
plot((eps_A-epsA_offset)*100,deviatoricS,'ro','displayn','Data'...
    ,'markersize',4); 
% *************Model Prediction********************
plot(-EPS(3,:)*100,Q,'r-','displayn','Deviator Stress'); 
set(gca,'ylim',[0,max(ylim)]);
yyaxis right
plot(-EPS(3,:)*100,D,'--','displayn','Damage Variable','linewidth',1);
legend('show')
xlabel('Strain (%)');
ylabel('Damage/Recovery Variable');
set(gca,'xlim',[0,max(xlim)]);


% ____Elastic Parameters_________
% figure; hold on;
% subplot(2,1,2); 
% plot(-EPS(3,:)*100,E11,'k--','displayn','E_{1}');
% hold on
% plot(-EPS(3,:)*100,E33,'k-','displayn','E_{3}');
% ylabel('Elastic Modulus (MPa)');
% xlabel('Axial Strain (%)');
% yyaxis right
% plot(-EPS(3,:)*100,NU12,'b-','displayn','\nu_{12}');
% plot(-EPS(3,:)*100,NU12,'m-','displayn','\nu_{21}');
% plot(-EPS(3,:)*100,NU13,'g-','displayn','\nu_{13}');
% plot(-EPS(3,:)*100,NU31,'r-','displayn','\nu_{31}');
% ylabel('Poissons Ratio');
% xlabel('Axial Strain (%)');
% legend('show');

rmpath('./Fcns');