%% Objective function for search of strain increment
% Solve the prescribed boundaries condition using Newton-Raphson method
% reference: IncrementalDriver2017 by Niemunis
% prescribed BC - STRAIN: prescribed(:,1),STRESS: prescribed(:,3)
function deps=NR(fcn_sig,prescribed,para,variables)
    deps_old=variables{1};
    eps_old=variables{3};
    stran_old=[eps_old(1,1),eps_old(2,2),eps_old(3,3),2*eps_old(1,2),2*eps_old(1,3),2*eps_old(2,3)]';
    isstran=find(prescribed(:,1)==0);
    isstress=find(prescribed(:,1)==1);
    dstran=[deps_old(1,1),deps_old(2,2),deps_old(3,3),2*deps_old(1,2),2*deps_old(1,3),2*deps_old(2,3)]';
    dstran(isstran)=prescribed(isstran,2)-stran_old(isstran);
    if length(isstran)<6 % mixed boundaries or pure stress boundaries
        tol=1e-5; 
        c_dstran=zeros(6,1);
        ddsdde_pr=zeros(length(isstress));
        err_norm_kp1=1;  %norm of error in the current step
        err_norm_k=2; %norm of error in the last step
        while (err_norm_kp1>tol)&&(abs((err_norm_k-err_norm_kp1))>tol)
%         while (err_norm_kp1>tol)
            err_norm_k=err_norm_kp1;
            dstran=dstran+c_dstran;
            deps=[dstran(1),dstran(4)/2,dstran(5)/2;...
                  dstran(4)/2,dstran(2),dstran(6)/2;...
                  dstran(5)/2,dstran(6)/2,dstran(3)];
            variables{1}=deps;
            % try with updated dstran
            Output=fcn_sig(variables{:},para);
            sig=Output{1};ddsdde=Output{2};
            stress=[sig(1,1),sig(2,2),sig(3,3),sig(1,2),sig(1,3),sig(2,3)]';
            % error between computed stress and prescribed stress
            err_dstress=prescribed(isstress,3)-stress(isstress);
            err_norm_kp1=norm(err_dstress,'fro');
            % construct the stiffness matrix to derive the dstran corrector
            for i=1:length(isstress)
                for j=1:length(isstress)
                    ddsdde_pr(i,j)=ddsdde(isstress(i),isstress(j));
                end
            end
            c_dstran(isstress)=ddsdde_pr^-1*err_dstress; % dstran corrector
        end
    end
    deps=[dstran(1),dstran(4)/2,dstran(5)/2;...
          dstran(4)/2,dstran(2),dstran(6)/2;...
          dstran(5)/2,dstran(6)/2,dstran(3)];
end
