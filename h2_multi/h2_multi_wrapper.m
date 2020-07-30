

% ------------------
% set up paths/dirs
% ------------------
source('/gpfs/milgram/project/holmes/kma52/h2_multi/R/h2_multi.R')
base_dir = '/gpfs/milgram/project/holmes/kma52/h2_multi'


K = csvread(fullfile(base_dir, 'example_data', 'K.csv'), 1);
X = csvread(fullfile(base_dir, 'example_data', 'covar.csv'), 1);
F = csvread(fullfile(base_dir, 'example_data', 'F.csv'), 1);
P = csvread(fullfile(base_dir, 'example_data', 'lh_overall_dice_P.csv'), 0);


% ------------------
% kinship matrix
% ------------------
kin   = csvread(fullfile(base_dir, 'example_data', 'K.csv'), 1);

% ------------------
% covariates
% ------------------
covar = csvread(fullfile(base_dir, 'example_data', 'covar.csv'), 1);


% ------------------
% family ids
% ------------------
F = csvread(fullfile(base_dir, 'example_data', 'F.csv'), 1);


% --------------------------------------------
% heritability of overall network topography
% --------------------------------------------
pheno = csvread(fullfile(base_dir, 'example_data', 'lh_overall_dice_P.csv'), 0);


% estimate heritability, permuted p-values
[h2, p_perm] = h2_multi(pheno, kin, covar, 1000);


% same as above, but get jack-knife SE
[h2, p_perm, jack_se] = h2_multi(pheno, kin, covar, 10, F);


