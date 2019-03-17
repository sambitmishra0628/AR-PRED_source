function [] = PredictActiveSiteResidues(outputdir, pdb_ca_file)
    % Make a consolidated file containing all the features for active site
    % prediction in their required order and then predict active site
    % residues using a set of trained models.
    
    % Weigh the probability scores for each residue from each model by the
    % individual model's MCC. The net score for a residue 'i' from all the 10 models is
    % then: sum(An*Pn)/sum(An), where An is the MCC of model 'n' and
    % Pn is the probability score assigned to residue 'i' by model 'n'.
    % This function is for a single protein file and we will write the
    % prediction output to a file.
    
    model_files = {'model_1.mat','model_2.mat','model_3.mat','model_4.mat','model_5.mat','model_6.mat','model_7.mat','model_8.mat','model_9.mat','model_10.mat'};
    model_wts = [0.6935;0.6296;0.6930;0.6329;0.7156;0.6697;0.6885;0.6769;0.7300;0.6752]; % consider weights as MCCs of individual models
    arpred_path = getenv('ARPRED_HOME');
    if isempty(arpred_path)
        error ('ERROR!!! ARPRED_HOME environment variable not set! Cannot continue with prediction!!!\n');
		exit;
    end
    model_path = strcat(arpred_path, '/', 'rf_models_active_site');
    feature_matrix = GetConsolidatedFeatureMatrix(outputdir); 
    size(feature_matrix)
    Y_score_all = []; % Store the scores for all residues
    Y_pred_all = []; % Store the predicted class for all residues
    model_count = 0;
    METRICS = [];
    for i=1:numel(model_files)
        model_i = strcat(model_path,'/',model_files{i});
        load (model_i);
        model_count = model_count+1;
        fprintf ('Running predictions with model %d...\n', model_count);
        %X_test = feature_matrix(:,1:37);
        X_test = feature_matrix;
        [Y_pred,Y_score] = predict(rf_model, X_test); % run predictions on the test set
        Y_pred = str2num([Y_pred{:}]');
        Y_score_all = [Y_score_all;[Y_score(:,2)]'];
        Y_pred_all = [Y_pred_all;Y_pred'];
        fprintf('Done!\n');
        clear rf_model;
    end
    [~, nres] = size(Y_score_all);
    
    % Weigh the residue probability scores by weights assigned to each model.
    model_wts_rep = repmat(model_wts,1,nres); 
    Y_score_weighted = Y_score_all.*model_wts_rep;
    Y_score_weighted_residue_level = sum(Y_score_weighted,1)./sum(model_wts); % sum elements from each column and normalize by the sum of model weights to get a single score for each residue
    [Y_score_weighted_residue_level_sorted, sorted_ind] = sort(Y_score_weighted_residue_level,'descend'); % Sort in descending order
    
    % Read the PDB file and get the PDB indices for residues
    P = readPDB(strcat(outputdir, '/', pdb_ca_file),1);
    resind = P.resSeq;
    resind_sorted = resind(sorted_ind); % Get the corresponding PDB indices for the residues
    all_predictions = [resind_sorted(:) Y_score_weighted_residue_level_sorted(:)];
    
    % Write the prediction output to file
    pred_output_file = strcat(outputdir, '/', 'activesite_predictions.csv');
    fh = fopen(pred_output_file, 'w');
    fprintf(fh, 'Residue, Score\n');
    fclose(fh);
    dlmwrite(pred_output_file, all_predictions, 'delimiter', ',', '-append');
    fprintf('Wrote weighted predictions and scores to file %s\n', pred_output_file);
end

%% 
% Combine all the features for active site residue prediction in the 
% expected order into a single matrix and return that matrix.
function [feature_matrix] = GetConsolidatedFeatureMatrix(outputdir)
    feature_file_names = {'aa_type.csv', 'aa_identity.csv', 'sasa_features.csv', 'sse.csv', 'msf.csv', 'aa_hpathy.csv', 'dfi_activesite.csv', 'conservation.csv', 'ntw_cen_adj.csv', ...
        'ntw_cen_dist_int.csv', 'ntw_cen_bt_activesite.csv', 'ntw_cen_corr_activesite.csv', 'pocket_features.csv'};
    feature_matrix = [];
    for i=1:length(feature_file_names)
        feature_file_i = strcat(outputdir, '/', cell2mat(feature_file_names(i)));
        rec = csvread(feature_file_i);
       feature_matrix = [feature_matrix;rec];
    end
    feature_matrix = feature_matrix'; % take the transpose so that features are columns
end



%%% LGPL 
%    This file is part of AR-PRED: Active and Regulatory site Prediction
%
%    AR-PRED is a free software. You can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    AR-PRED is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    and the GNU Lesser General Public License along with the AR-PRED source code.
%    If not, see <http://www.gnu.org/licenses/>.
    
