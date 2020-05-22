% DEMO for testing RDA on C1-C2
clc;clear
addpath('liblinear-2.1/matlab');

% Preprocess data using L2-norm
src   = 'COIL_1';
tgt   = 'COIL_2';
fprintf(' %s vs %s ', src, tgt);
load(['data/' src '.mat']);
% Normalization
X_src = X_src./repmat(sum(X_src,2),1,size(X_src,2));
X_tar = X_tar./repmat(sum(X_tar,2),1,size(X_tar,2));
% Z-score
X_src = zscore(X_src',1)';
X_tar = zscore(X_tar',1)';
% L2-norm
X_src = X_src*diag(sparse(1./sqrt(sum(X_src.^2))));
X_tar = X_tar*diag(sparse(1./sqrt(sum(X_tar.^2))));
Xs    = X_src;
Ys    = Y_src;
Xt    = X_tar;
Yt    = Y_tar;

% Set hyderparameters
options.kernel    = 'primal'; % 'primal' or 'linear'
options.d         = 20;

options.alpha     = 10;
options.lambda    = 0.1;
options.lambdaMMD = 1;   
options.Ipara     = 1; 

if strcmp(options.kernel,'primal')
    options.T = 10; 
else
    options.T = 13; 
end
options.Td        = 10;
options.Tg        = 30;
options.S         = 4;
fprintf('alpha = %d , lambda = %d , lambdaMMD = %d ,Ipara_all = %d \n', options.alpha, options.lambda, options.lambdaMMD, options.Ipara);

% do RDA
[Acc,~,~,~,~] = RDA(Xs, Ys, Xt, Yt, options);    


