function [beta1,beta2]=updateD(Ps,Pt,Xs,Xt,Ys,beta1,options)   
lambda  = options.lambda;
Td      = options.Td;
d       = options.d;
class   = length(unique(Ys));
ns      = size(Xs,2);
nt      = size(Xt,2);

rho     = 1;
A       = Ps'*Xs;
B       = Pt'*Xt;
Y       = Construct_Y(Ys,length(Ys),class); Y=Y';
beta2   = beta1;
Spar    = 1;
for iterd = 1:Td  
    % checking convergence
    conv(iterd) = Spar*norm(A'*beta1-Y,'fro')+Spar*norm(A'*beta2-Y,'fro')-lambda*norm(B'*beta1-B'*beta2,'fro');
    
    % update beta2
    beta2 = (Spar*A*A'-lambda*B*B'+rho*eye(d))\(Spar*A*Y-lambda*B*B'*beta1);

    % update beta1
    beta1 = (Spar*A*A'-lambda*B*B'+rho*eye(d))\(Spar*A*Y-lambda*B*B'*beta2);

    if iterd > 2
         if  abs((conv(iterd)-conv(iterd-1)))< 10^-6
              break
         end
    end
end
end

