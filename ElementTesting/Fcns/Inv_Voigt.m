% Inverse Voigt Transformation
function T3333=Inv_Voigt(M66,ICOE)
T3333=zeros(3,3,3,3);
switch ICOE
    case 1 % Stiffness matrix
        COE1=1;
        COE2=2;
    case 2 % Compliance matrix
        COE1=2;
        COE2=4;
end

I=[1,4,6;4,2,5;6,5,3];
J=[1,4,6;4,2,5;6,5,3];

COE=[1,1,1,COE1,COE1,COE1;1,1,1,COE1,COE1,COE1;1,1,1,COE1,COE1,COE1;...
    COE1,COE1,COE1,COE2,COE2,COE2;COE1,COE1,COE1,COE2,COE2,COE2;COE1,COE1,COE1,COE2,COE2,COE2];

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                T3333(i,j,k,l)=M66(I(i,j),J(k,l))/COE(I(i,j),J(k,l));
            end
        end
    end
end
end
        
    
