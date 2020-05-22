function beta1 = Initbeta1(Ps,Xs,Ys,d)
class = length(unique(Ys));
Y = Construct_Y(Ys,length(Ys),class); Y=Y';

%update beta1
A=Ps'*Xs;
rho=1;%[0,0.01,0.1,1,10];
beta1 = (A*A'+rho*eye(d))\(A*Y);
end