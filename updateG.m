function [Ps,Pt]=updateG(Ps0,Pt0,Xs,Xt,Ys,Yt_f,beta1,beta2,options) 
%用梯度下降的方式求
d       = options.d;
alpha   = options.alpha;
lambda  = options.lambda;
lambdaMMD=options.lambdaMMD;


Tg= options.Tg;
class = length(unique(Ys));

[m,ns] = size(Xs); nt = size(Xt,2);
convergence = 10^-6;
rhoPt=0.1;
% ----------------------------------------------
%               Initialization
% ----------------------------------------------
[Ts, Tt, Tst, Tts] = constructMMD1(ns,nt,Ys,Yt_f,class);
Ms = Xs*Ts*Xs';
Mt = Xt*Tt*Xt';
Mst = Xs*Tst*Xt';
Mts = Xt*Tts*Xs';
% ------------------------------------------------
%                   Main Loop
% ------------------------------------------------
conv(1)=100;
for iter = 1:Tg      
    % updating Pt 
    if (iter == 1)
        Pt = Pt0;
    else
        A  = 2*alpha*eye(m)+2*lambdaMMD*Mt;
        B  = 2*lambda*Xt*Xt';
        C  = beta1*beta1'-beta2*beta1'-beta1*beta2'+beta2*beta2';
        D  = 2*alpha*Ps-lambdaMMD*Mts*Ps-lambdaMMD*Mst'*Ps;
        deltaPt = A*Pt+B*Pt*C-D;
        Pt = Pt-rhoPt*deltaPt/norm(deltaPt,'fro');
    end
    
%     % updating Ps
%     if (iter == 1)
%         Ps = Ps0;
%     else
%          Ps = (2*alpha*eye(m)+lambdaMMD*Ms)\(2*alpha*Pt-lambdaMMD*Mts'*Pt-lambdaMMD*Mst*Pt);
         Ps = (2*alpha*eye(m)+2*lambdaMMD*Ms)\(2*alpha*Pt-lambdaMMD*Mts'*Pt-lambdaMMD*Mst*Pt);
%     end       

    % checking convergence
    conv(iter)=alpha*norm(Ps-Pt,'fro')+lambda*norm((Pt'*Xt)'*beta1-(Pt'*Xt)'*beta2,'fro')+ lambdaMMD*trace([Ps',Pt']*[Ms,Mst;Mts,Mt]*[Ps;Pt]);
    cov1(iter)=alpha*norm(Ps-Pt,'fro');
    cov2(iter)=lambda*norm((Pt'*Xt)'*beta1-(Pt'*Xt)'*beta2,'fro');
    
    cov3(iter)=lambdaMMD*trace([Ps',Pt']*[Ms,Mst;Mts,Mt]*[Ps;Pt]);
    if iter >= 2
         cov2_diff(iter-1)= abs((cov2(iter)-cov2(iter-1)));
         if abs((conv(iter)-conv(iter-1)))<convergence
              break
         end
    end    
end
end

function [Ms, Mt, Mst, Mts] = constructMMD1(ns,nt,Ys,Yt0,C)
e = [1/ns*ones(ns,1);-1/nt*ones(nt,1)];
es = 1/ns*ones(ns,1);
et = -1/nt*ones(nt,1);

M = e*e'*C;
% Ms = es*es'*C;
% Mt = et*et'*C;
% Mst = es*et'*C;
% Mts = et*es'*C;
Ms=zeros(ns);Mt=zeros(nt);Mst=zeros(ns,nt);Mts=zeros(nt,ns);
if ~isempty(Yt0) && length(Yt0)==nt
    for c = reshape(unique(Ys),1,C)
        es = zeros(ns,1);
        et = zeros(nt,1);
        es(Ys==c) = 1/length(find(Ys==c));
        et(Yt0==c) = -1/length(find(Yt0==c));
        es(isinf(es)) = 0;
        et(isinf(et)) = 0;
        Ms  = Ms + es*es';
        Mt  = Mt + et*et';
        Mst = Mst + es*et';
        Mts = Mts + et*es';
    end
end

Ms = Ms/norm(M,'fro');
Mt = Mt/norm(M,'fro');
Mst = Mst/norm(M,'fro');
Mts = Mts/norm(M,'fro');
end