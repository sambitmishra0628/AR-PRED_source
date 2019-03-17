function [errorflag] = matlab_feature_calc_wrapper(pdb_ca_file, outputdir, pred_type, active_site_file)
    % Wrapper to call each feature calculation function and write the calculated
    % features into the respective feature files. 
    % Before running this module, make sure that the paths to all the matlab modules
    % provided in this tool are included in the search path.
    %
    % Arguments-
    %   pdb_ca_file - The C-alpha PDB file of the protein
    %	outputdir - The name of the output directory
    %   pred_type - 1 for active site prediction, 
    %               any other value for allosteric site residue
    %               prediction
    %   active_site_file - CSV file having the list of residues forming the
    %   active site. Only required if pred_type is not 1 i.e., for prediction
    %   of allosteric site residues.

    numargs = nargin;
    if numargs < 3 || isempty(pdb_ca_file) || isempty(outputdir) || isempty(pred_type) % Ensure that all the required arguments are passed
        disp('Error!!! Missing required arguments!\n');
        return;
    end
	logfile = strcat(outputdir, '/', 'log.txt');
	pdb_ca_file = strcat(outputdir, '/', pdb_ca_file);
    errorflag = 0;
    % Open the log file to write the progress
    [LH, LH_msg] = fopen(logfile, 'a');
    if ~isempty(LH_msg) % Some problem opening the log file
		fprintf('ERROR opening log file: %s', LH_msg);
		errorflag = 1;
		return;
    end
    if exist('active_site_file', 'var') && ~isempty(active_site_file)
        active_site_ind = ParseActiveSiteFile(active_site_file);
    else
        active_site_ind = [];
    end
    
    % We will perform calculations for the individual features and if in
    % any calculation we encounter error, we will write it into the log
    % file and exit the calculations.
    
	% Perform calculations for amino acid type
    feature_file = strcat(outputdir, '/', 'aa_type.csv');
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for amino acid type as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for amino acid type as feature file already exists!\n');
    else    
        fprintf ('Mapping amino acid to respective type ...\n');
        fprintf (LH, 'Mapping amino acid to respective type ...\n');
        try
            GetAAType(pdb_ca_file, feature_file);
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function GetAAType. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function GetAAType: %s\n', ERROR.message);
            return
        end    
        fprintf ('Done!\n');    
        fprintf (LH, 'Done!\n');    
    end
    
	% Perform calculations for amino acid identity
	feature_file = strcat(outputdir, '/', 'aa_identity.csv');
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for amino acid identity as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for amino acid identity as feature file already exists!\n');
    else    
        fprintf ('Mapping amino acid to its respective identity ...\n');
        fprintf (LH, 'Mapping amino acid to its respective identity ...\n');
        try
            AminoAcidIdentity(pdb_ca_file, feature_file);
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function AminoAcidIdentity. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function AminoAcidIdentity: %s\n', ERROR.message);
            return
        end
        fprintf ('Done!\n');    
        fprintf (LH, 'Done!\n');    
    end
    
	% Perform calculations for amino acid hydropathy using Kyte-Doolittle index
	feature_file = strcat(outputdir, '/', 'aa_hpathy.csv'); 
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for amino acid hydropathy index as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for amino acid hydropathy index as feature file already exists!\n');
    else    
        fprintf ('Mapping amino acid to its respective hydropathy index ...\n');
        fprintf (LH, 'Mapping amino acid to its respective hydropathy index ...\n');
        try
            GetSequenceHydropathy(pdb_ca_file, feature_file);
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function GetSequenceHydropathy. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function GetSequenceHydropathy: %s\n', ERROR.message);
            return
        end
        fprintf ('Done!\n');    
        fprintf (LH, 'Done!\n'); 
    end
    
	% Perform calculations for residue mean-squared fluctuations
    feature_file = strcat(outputdir, '/', 'msf.csv'); 
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for residue mean square fluctuations as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for residue mean square fluctuations as feature file already exists!\n');
    else    
        fprintf('Running calculations for residue mean square fluctuations...\n');
        fprintf(LH, 'Running calculations for residue mean square fluctuations...\n');
        try
            GetResidueMSF(pdb_ca_file, feature_file); % Run calculations for residue mean-square fluctuations
        catch ERROR    
            errorflag = 1;
            fprintf ('ERROR encountered in function GetResidueMSF. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function GetResidueMSF: %s\n', ERROR.message);
            return
        end    
        fprintf ('Done!\n');    
        fprintf (LH, 'Done!\n'); 
    end
    
	% Perform calculations for residue dynamic flexibility index
    if pred_type == 1
        feature_file = strcat(outputdir, '/', 'dfi_activesite.csv');
    else
        feature_file = strcat(outputdir, '/', 'dfi_allostericsite.csv');
    end
    
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for perturbation response as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for perturbation response as feature file already exists!\n');
    else    
        fprintf('Running calculations for perturbation response...\n');
        fprintf(LH, 'Running calculations for perturbation response...\n');
      
        try
            DFIandAlloResponse(pdb_ca_file, feature_file, active_site_ind); % Run calculations for DFI and allosteric response
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function DFIandAlloResponse. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function DFIandAlloResponse: %s\n', ERROR.message);
            return
        end
        fprintf ('Done!\n');    
        fprintf (LH, 'Done!\n'); 
    end
    
    
	% Perform calculations for the different network centrality features
    fprintf('Running calculations for network centralities...\n');
    fprintf(LH, 'Running calculations for network centralities...\n');
    
    feature_file = strcat(outputdir, '/', 'ntw_cen_adj.csv');
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for network centrality with adjacency matrix as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for network centrality with adjacency matrix as feature file already exists!\n');
    else    
        try
            NetworkCentralities_AdjMat(pdb_ca_file, feature_file); % Run calculations for network centralities with distance based adjacency matrix
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function NetworkCentralities_AdjMat. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function NetworkCentralities_AdjMat: %s\n', ERROR.message);
            return
        end
    end
    
    feature_file = strcat(outputdir, '/', 'ntw_cen_dist_int.csv');
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for network centrality with inverse distance matrix as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for network centrality with inverse distance matrix as feature file already exists!\n');
    else    
        try
            NetworkCentralities_DistBasedInteractionMatrix(pdb_ca_file, feature_file); % Run calculations for network centralities using the distance-based interaction matrix
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function NetworkCentralities_DistBasedInteractionMatrix. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function NetworkCentralities_DistBasedInteractionMatrix: %s\n', ERROR.message);
            return
        end
    end
    
    if pred_type == 1
        feature_file = strcat(outputdir, '/', 'ntw_cen_bt_activesite.csv'); 
    else
        feature_file = strcat(outputdir, '/', 'ntw_cen_bt_allostericsite.csv'); 
    end
    
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for network centrality with BT potential matrix as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for network centrality with BT potential matrix as feature file already exists!\n');
    else    
        try
            NetworkCentralities_BT_Interactions(pdb_ca_file, feature_file, pred_type); % Run calculations for network centralities for the network created with BT potential
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function NetworkCentralities_BT_Interactions. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function NetworkCentralities_BT_Interactions: %s\n', ERROR.message);
            return
        end
    end    
    
    if pred_type == 1
        feature_file = strcat(outputdir, '/', 'ntw_cen_corr_activesite.csv'); 
    else
        feature_file = strcat(outputdir, '/', 'ntw_cen_corr_allostericsite.csv'); 
    end
    
    if exist(feature_file, 'file')
        fprintf('Skipping calculations for network centrality with dynamic correlation matrix as feature file already exists!\n');
        fprintf(LH, 'Skipping calculations for network centrality with dynamic correlation matrix as feature file already exists!\n');	
    else    
        try
            NetworkCentralities_DynamicCorr_DistTransformed(pdb_ca_file, feature_file, pred_type); % Run calculations for network centralities for the network of dynamic cross-correlations
        catch ERROR
            errorflag = 1;
            fprintf ('ERROR encountered in function NetworkCentralities_DynamicCorr_DistTransformed. Check log file for details\n');
            fprintf (LH, 'ERROR encountered in function NetworkCentralities_DynamicCorr_DistTransformed: %s\n', ERROR.message);
            return
        end    
    end
    
	fprintf('Done!\n');	
    fprintf(LH,'Done!\n');
    
    % For allosteric site prediction, we will run the calculations for
    % shortest path to active site residues.
    if pred_type ~= 1
        feature_file = strcat(outputdir, '/', 'shortest_path.csv');
        if exist(feature_file, 'file')
            fprintf('Skipping calculations for shortest dynamically correlated path between allosteric and active site residuesx as feature file already exists!\n');
            fprintf(LH, 'Skipping calculations for shortest dynamically correlated path between allosteric and active site residues as feature file already exists!\n');
        else    
            fprintf('Running calculations for shortest dynamically correlated path between allosteric and active site residues\n');
            fprintf(LH, 'Running calculations for shortest dynamically correlated path between allosteric and active site residues\n');
            try
                [active_site_ind] = ParseActiveSiteFile(active_site_file);
                CalculateShortestPathToActiveSiteRes(pdb_ca_file, feature_file, active_site_ind);
            catch ERROR
                errorflag = 1;
                fprintf('ERROR encountered in function CalculateShortestPathToActiveSiteRes. Check log file for details\n');
                fprintf(LH, 'ERROR encountered in function CalculateShortestPathToActiveSiteRes : %s\n', ERROR.message);
                return
            end    
            fprintf('Done!\n');	
            fprintf(LH,'Done!\n');
        end   
    end    
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


