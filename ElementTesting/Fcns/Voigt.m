% Voigt transformation
function M66=Voigt(T3333,ICOE)

% M66=zeros(6,6);
switch ICOE
    case 1 % Stiffness matrix
        COE1=1;
        COE2=1;
    case 2 % Compliance matrix
        COE1=2;
        COE2=4;
end

i=[1,2,3,1,2,1];
j=[1,2,3,2,3,3];
k=[1,2,3,1,2,1];
l=[1,2,3,2,3,3];

COE=[1,1,1,COE1,COE1,COE1;1,1,1,COE1,COE1,COE1;1,1,1,COE1,COE1,COE1;...
    COE1,COE1,COE1,COE2,COE2,COE2;COE1,COE1,COE1,COE2,COE2,COE2;COE1,COE1,COE1,COE2,COE2,COE2];

for I=1:6
    for J=1:6
        M66(I,J)=T3333(i(I),j(I),k(J),l(J))*COE(I,J);
    end
end

end