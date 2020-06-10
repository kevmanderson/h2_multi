

% ------------------
% set up paths/dirs
% ------------------
addpath('/gpfs/milgram/project/holmes/kma52/h2_multi/matlab')
base_dir = '/gpfs/milgram/project/holmes/kma52/h2_multi'


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
covar = csvread(fullfile(base_dir, 'example_data', 'F.csv'), 1);


% --------------------------------------------
% heritability of overall network topography
% --------------------------------------------
pheno = csvread(fullfile(base_dir, 'example_data', 'lh_overall_dice_P.csv'), 0);


% estimate heritability, permuted p-values
[h2, p_perm] = h2_multi(pheno, kin, covar, 1000);


    for (i = 1:17)
        disp(i)
        pheno = csvread(fullfile(base_dir, [hemi '_dice_net_', num2str(i),'_P.csv']), 0);
        [h2, p_perm]   = h2_mat(pheno, kin, covar, 1000);
        out_mat(i+1,:) = {h2, p_perm, num2str(i), hemi};
    end

    out_path = fullfile(base_dir, [hemi '_dice_network_topology_h2.csv']);
    writetable(out_mat, out_path)
end


pheno = csvread(fullfile(base_dir, 'height_P.csv'), 0);
[h2, p_perm] = h2_mat(pheno, kin, covar, 1000)

out_mat = {};
pheno = csvread(fullfile(base_dir, 'dice_net_12_P.csv'), 0);
[h2, p_perm] = h2_mat(pheno, kin, covar, 1000);
cur_row = [h2, p_perm];