% Elastoplastic damage model with elastic stiffness depedent on stress
% tensor (AniDM)

function Output=fcn_DPIsodwAniDM(deps,sigeff_old,eps_old,eps_e_old,D,alpha,para)

% elastic parameters
E1=para{1};
nu=para{2};
% damage evolution parameters
c0d=para{4};
dmax=para{5};
c1d=para{6};
% stress dependency parameters
b=para{7};
% drucker prager parameters
kappa=para{3};
B=para{8};

elambda1=E1*nu/((1+nu)*(1-2*nu));emu1=E1/(2*(1+nu));
eK1=elambda1+2/3*emu1;
C0=zeros(3,3,3,3); 
I2=eye(3,3);

% INITIAL STIFFNESS
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                I2I2(i,j,k,l)=I2(i,j)*I2(k,l);
                I4(i,j,k,l)=0.5*(I2(i,k)*I2(j,l)+I2(i,l)*I2(j,k));
                C0(i,j,k,l)=elambda1*I2I2(i,j,k,l)+2*emu1*I4(i,j,k,l);
            end
        end
    end
end
C099=reshape(C0,[9,9]);

% 1. EVALUATE DAMAGE FUNCTION and UPDATE DAMAGE VARIABLE
eps=eps_old+deps;
[eps_plus,~]=PosProj(eps); % as a function of total strain
norm_eps_plus=norm(eps_plus,'fro');

% exponential function with L2-norm of tensile strain and limit of
% maximum damage variable
fd=dmax*(1-exp((c0d-norm_eps_plus)/c1d))-D;
if fd>0 % damage activated
    D=dmax*(1-exp((c0d-norm_eps_plus)/c1d)); 
end
kappa_d=(1-D)*kappa;


% 2. UPDATE STRESS-DEPENDENT STIFFNES

% I) Iterate to find alpha at the current step
del_R=1;
tol_R=1e-5;
while abs(del_R)>tol_R

% II) Use the alpha from the last increment to compute the C_sd

    alpha_tr=alpha;
    [v_P,~]=eig(sigeff_old);
    
    Q=zeros(3,3);
    Qk=zeros(3,3);
    for K=1:3
        Q=Q+kron(v_P(:,K),v_P(:,K)');
        Qk=Qk+(alpha_tr(K))^(1/4)*kron(v_P(:,K),v_P(:,K)');
    end
    
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    P(i,j,k,l)=0;
                    for m=1:3
                        for n=1:3
                            P(i,j,k,l)=P(i,j,k,l)+Qk(i,m)*Qk(j,n)*Q(k,m)*Q(l,n);
                        end
                    end
                end
            end
        end
    end

    P99=reshape(P,[9,9]);
    Csd99=P99*(C099)*P99; 
    C_sd=reshape(Csd99,[3,3,3,3]);

    eps_e=eps_e_old+deps;
    sigeff_tr=Aijkl_Bkl(C_sd,eps_e);
    peff_tr=trace(sigeff_tr)/3;
    Seff_tr=sigeff_tr-peff_tr*I2;

% 3. EVALUATE YIELD FUNCTION AND UPDATE PLASTIC STRAIN
    norm_Seff_tr=norm(Seff_tr,'fro');
    fp=sqrt(3/2)*norm_Seff_tr+B*peff_tr-kappa_d; 
    eps_p_old=eps_old-eps_e_old;
    
    if fp<=0 % in elastic domain
        sigeff=sigeff_tr;
        eps_p=eps_p_old;

    else % yield

        % computation constants
        lambda=0;
        P_dev=I4-1/3*I2I2; %deviatoric projection tensor
        tol_fp=1e-8;
        tol_Resid=1e-8;
        C_sd66=Voigt(C_sd,1);
        S_sd66=C_sd66^-1;
        S_sd=Inv_Voigt(S_sd66,2);


        % initialize all the variables used in the iteration:
        eps_p=eps_p_old;
        peff=peff_tr;
        Seff=Seff_tr;
        sigeff=Seff+peff*I2;

        SeffSeff=Aij_Bkl(Seff,Seff);
        norm_Seff=norm(Seff,'fro');
        n_hat=Seff/norm_Seff;
        df_dsig=sqrt(3/2)*n_hat+B/3*I2;
        df2_dsig2=sqrt(3/2)*(P_dev/norm_Seff-SeffSeff/norm_Seff^3);
        Resid=-eps_p+eps_p_old+lambda*df_dsig;
        norm_Resid=norm(Resid,'fro');

        % General return-mapping algorithm (Simo & Hughes, 1997)
        while fp>tol_fp||(norm_Resid>tol_Resid)||(lambda<0)

            
            %step1
            Xi=Inv_Voigt((S_sd66+lambda*Voigt(df2_dsig2,2))^-1,1);
            lambda2=(fp-Aij_Bij(Resid,Aijkl_Bkl(Xi,df_dsig)))/Aij_Bij(df_dsig,Aijkl_Bkl(Xi,df_dsig));
            
            %step2
            del_sigeff=Aijkl_Bkl(Xi,-Resid-lambda2*df_dsig);
            del_eps_p=Aijkl_Bkl(-S_sd,del_sigeff);
            
            %step3
            eps_p=eps_p+del_eps_p;
            lambda=lambda+lambda2;
            sigeff=sigeff+del_sigeff;

            %evaluate the termination criteria
            peff=trace(sigeff)/3;
            Seff=sigeff-peff*I2;
            SeffSeff=Aij_Bkl(Seff,Seff);
            norm_Seff=norm(Seff,'fro');
            n_hat=Seff/norm_Seff;
            df_dsig=sqrt(3/2)*n_hat+B/3*I2;
            df2_dsig2=sqrt(3/2)*(P_dev/norm_Seff-SeffSeff/norm_Seff^3);
            fp=sqrt(3/2)*norm_Seff+B*peff-kappa_d;  
            Resid=-eps_p+eps_p_old+lambda*df_dsig;
            norm_Resid=norm(Resid,'fro');

        end
    end
    
    eps_e=eps-eps_p;
    
    [v_P,~]=eig(sigeff);
    for K=1:3
        sig_K=v_P(:,K)'*sigeff*v_P(:,K);
        if sig_K>-0.01
            sig_K=-0.01;
        end

        alpha(K)=(-sig_K)^b;
    end
    
% for "while abs(del_R)>tol_R"    
    del_R=norm(alpha_tr-alpha,'fro');
end % 


% 7. UPDATE STIFFNESS AND STRESS
C=C_sd*(1-D);
sig=Aijkl_Bkl(C,eps_e);
C66=Voigt(C,1);

% For undrained loading (Ref: lecture note from CIV E 799)
% Kf=2.1e3; % bulk modulus for the fluid, MPa
% n_poro=0.1; % porosity of the rock
% m_vec=[1;1;1;0;0;0]; % vector to calculate trace from a vector
% C66=C66+m_vec*(Kf/n_poro)*m_vec';
% C=Inv_Voigt(C66,1);
% sig=Aijkl_Bkl(C,eps_e);

ddsdde=C66/(1-b);
Output={sig,ddsdde,C66,deps,sigeff,eps,eps_e,D,alpha};


%%-------------------BELOW ARE FUNCTIONS-----------------------------------
% calculate the positive projection of strain
function [eps_plus,P_plus]=PosProj(eps)
[v_eps,eps_k]=eig(eps);

% calculate the generalized tensile strain (Ortiz, 1985)
eps_plus=zeros(3);
for I=1:3
    eps_plus=eps_plus+eps_k(I,I)*heaviside(eps_k(I,I))*(v_eps(:,I)*v_eps(:,I)');
end

% compute with the positive projection tensor (Hansen & Schreyer, 1995)
Q_plus=zeros(3);
for I=1:3
    Q_plus=Q_plus+heaviside(eps_k(I,I))*(v_eps(:,I)*v_eps(:,I)');
end
P_plus=Aik_Bjl(Q_plus,Q_plus);
% eps_plus=Aijkl_Bkl(P_plus,eps);
end

function C=Aij_Bij(A,B) % dyadic composition product
C=0;
for i=1:3
    for j=1:3
        C=C+A(i,j)*B(i,j);
    end
end
end

function C=Aij_Bkl(A,B) % dyadic composition product
C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l)=A(i,j)*B(k,l);
            end
        end
    end
end
end

function C=Aik_Bjl(A,B) % cross composition product
C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l)=A(i,k)*B(j,l);
            end
        end
    end
end
end

function C=Ail_Bjk(A,B) % cross composition product
C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l)=A(i,l)*B(j,k);
            end
        end
    end
end
end


function C=Aijkl_Bkl(A,B)
C=zeros(3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j)=C(i,j)+A(i,j,k,l)*B(k,l);
            end
        end
    end
end

end

function C=Aijmn_Bmnkl(A,B)
C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,m,n)*B(k,l,m,n);
                    end
                end
            end
        end
    end
end
end


end


