function  [Acc,Ps,Pt,beta1,beta2] = RDA(Xs,Ys,Xt,Yt,options)
d        = options.d;
max_iter = options.T;
S        = options.S;
[m,nt]   = size(Xt);

%initializing Ps,Pt
Pss = pca(Xs, nt);
Ptt = pca(Xt, nt);
Pt0 = Ptt(:, 1:d);
Ps0 = Pss(:, 1:d);
Ps0 = Ps0*Ps0'*Pt0;% subspace alignment

if strcmp(options.kernel,'primal')
    for outloop = 1:max_iter
        if outloop == 1
            Ps=Ps0;Pt=Pt0;
            [Acc(outloop),Yt_f] = labeling_target(Xs, Xt, Ys, Yt, eye(m),eye(m), S);
        else 
            [Acc(outloop),Yt_f] = labeling_target(Xs, Xt, Ys, Yt, Ps, Pt, S);  
        end

        % init beta1 
        % beta1 = eye(d,20);
        beta1 = Initbeta1(Ps, Xs, Ys, d);

        % update D:beta1 and beta2
        [beta1,beta2] = updateD(Ps, Pt, Xs, Xt, Ys, beta1, options) ; 

        % update G:Ps and Pt
        [Ps,Pt] = updateG(Ps, Pt, Xs, Xt, Ys, Yt_f, beta1, beta2, options) ; 
    end
else
    X     = [Xs,Xt];
    X     = X*diag(1./sqrt(sum(X.^2)));
    [~,n] = size(X);
    ker   = options.kernel;
    gamma = 1;
    K     = kernel(ker, X, [], gamma);
    Ks    = kernel(ker, X, Xs, gamma);
    Kt    = kernel(ker, X, Xt, gamma);

    for outloop = 1:max_iter
        if outloop == 1 
            [Acc(outloop),Yt_f] = labeling_target(Xs, Xt, Ys, Yt, eye(m),eye(m), S);   
            Ps0 = X\Ps0;
            Pt0 = X\Pt0;
            Ps  = Ps0;
            Pt  = Pt0;
        else
            [Acc(outloop),Yt_f] = labeling_target(Ks, Kt, Ys, Yt, Ps, Pt, S);
        end
        
        %init beta1 
        beta1 = Initbeta1(Ps, Ks, Ys, d);

        %update D:beta1 and beta2
        [beta1,beta2] = updateD(Ps, Pt, Ks, Kt, Ys, beta1, options) ;
 
        %update G:Ps and Pt
        [Ps,Pt] = updateG_ker(Ps, Pt, Ks, Kt, K, Ys, Yt_f, beta1, beta2, options) ; 
    end
end
end

function [acc,Yt] = labeling_target(Xs, Xt, Xs_label, Xt_label, Ps, Pt,S)
C = [0.001 0.01 0.1 1.0 10 100 1000 10000];  
for chsvm = 1 :length(C)
    tmd = ['-s ' num2str(S) ' -c ' num2str(C(chsvm)) ' -B 1 -q'];
    model(chsvm) = train(Xs_label, sparse(double(Xs'*Ps)),tmd);
    [~,acc, ~] = predict(Xt_label, sparse(double(Xt'*Pt)), model(chsvm), '-q');
    [~,accs, ~] = predict(Xs_label, sparse(double(Xs'*Ps)), model(chsvm), '-q');
    acc1(chsvm)=acc(1);
end	
[acc,bestsvm_id]=max(acc1);
fprintf(' svm acc=%2.2f %%\n',acc);
model=model(bestsvm_id);
c=C(bestsvm_id);
score = model.w * [Pt'*Xt; ones(1, size(Xt, 2))];

[~, C] = max(score, [], 1);
Yt = C';
end