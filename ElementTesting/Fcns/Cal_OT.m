% Convert the orthotropic stiffness matix into engineering parameters
function Epara=Cal_OT(C)
S=C^-1;
E11=1/S(1,1);
E22=1/S(2,2);
E33=1/S(3,3);
nu12=-S(2,1)*E11;
nu21=-S(1,2)*E22;
nu13=-S(3,1)*E11;
nu31=-S(1,3)*E33;
nu23=-S(3,2)*E22;
nu32=-S(2,3)*E33;
G12=1/S(6,6);
G13=1/S(5,5);
G23=1/S(4,4);
Epara=[E11,E22,E33,G12,G13,G23;nu12,nu21,nu13,nu31,nu23,nu32];
end