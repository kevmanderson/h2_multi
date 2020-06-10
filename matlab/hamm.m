function [h2, p_perm] = h2_mat(P, K, X, n_perm)
%%
%
% Heritability estimation from phenotypic and genetic similary mattrices
%
% Input:
% P: an Nsubj x Nsubj phenotypic similarity matrix
% K: an Nsubj x Nsubj genetic similarity matrix
% X: an Nsubj x Ncov matrix of covariates
% n_perm: number of permutations; set Nperm = 0 if permutation inference is not needed
%
% Output:
% h2: heritability estimate
% p_perm: nonparametric permutation p-value; if Nperm = 0, PermPval = NaN
%
%%

[n_subj, n_cov] = size(X);

% -----
P0 = eye(n_subj) - X/(X'*X)*X';
[U,~,~] = svd(P0); U = U(:,1:n_subj-n_cov);
P = U'*P*U; K = U'*K*U;
n_subj = n_subj-n_cov;
%%
kappa = trace(K^2)/n_subj;
tau = trace(K)/n_subj;
vK = n_subj*(kappa-tau^2);

Qg = K-tau*eye(n_subj);
Qe = kappa*eye(n_subj)-tau*K;

tg = trace(Qg*P)/vK;
te = trace(Qe*P)/vK;
tp = tg+te;

h2 = max(min(tg/tp,1),0);
%%
if n_perm == 0
    p_perm = NaN;
    return
end

h2_perm = zeros(1,n_perm);
for s = 1:n_perm
    disp(['----- Permutation-', num2str(s), ' -----'])

    if s == 1
        K_perm = K;
    else
        subj_perm = randperm(n_subj);
        K_perm = K(subj_perm,subj_perm);
    end

    Qg_perm = K_perm-tau*eye(n_subj);
    Qe_perm = kappa*eye(n_subj)-tau*K_perm;

    tg_perm = trace(Qg_perm*P)/vK;
    te_perm = trace(Qe_perm*P)/vK;

    h2_perm(s) = max(min(tg_perm/(tg_perm+te_perm),1),0);
end

p_perm = sum(h2_perm>=h2)/n_perm;
%%

