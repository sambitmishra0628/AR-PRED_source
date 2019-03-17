function [] = NetworkCentralities_DynamicCorr_DistTransformed(pdb_ca_file, feature_file, pred_type)
    % Compute the following network properties by converting the 
    % matrix of dynamic cross-correlations derived from Elastic Network Models
    % into a distanc matrix. 
    % 1. Betweenness centrality
    % 2. Eigen centrality
    % 3. Closeness centrality
    % 4. Degree centrality
    % 5. Page rank centrality
	
    IND = readPDB(pdb_ca_file);
    enm_corr_mat = ENM_Dynamic_Correlation(IND);
    
	% Convert to a distance matrix
    enm_corr_mat = enm_corr_mat ./ max(abs(enm_corr_mat(:))); % Normalize to ensure the value lies between -1 and 1
    enm_corr_mat_dist = 1-enm_corr_mat;
    
    % Get the raw centrality values. We will use the inbuilt MATLAB
    % 2017 functions to calculate the betweenness, closeness, degree,
    % page rank and eigen centralities.

    G = graph(enm_corr_mat_dist);
	
	% We will consider the weights equivalent to the distance
	bc = G.centrality('betweenness','Cost',G.Edges.Weight);

	% Normalize the raw numbers returned in bc
	n_nodes = numnodes(G);
	norm_factor = ((n_nodes-2)*(n_nodes-1))/2; % normalize to make the bc as probabilities
	bc = bc ./ norm_factor;
    
	% Compute the closeness centrality
    cc = G.centrality('closeness','Cost',G.Edges.Weight);
    
    % Compute the eigen-vector centrality. For the eigen, degree and
    % page rank centralities however, we need to use the edge importance
    % which corresponds to the inverse of the distance
    ec = G.centrality('eigenvector','Importance',1./G.Edges.Weight);
    
    % Compute the degree centrality
    dc = G.centrality('degree','Importance',1./G.Edges.Weight);
    
    % Compute the page rank centrality
    pc = G.centrality('pagerank','Importance',1./G.Edges.Weight);
    

    % Write the calculated feature vectors into the respective feature
    % files based on the predictor type. 
    % pred_type = 1 : active site prediction
    % pred_type = 2 : allosteric site prediction
    if pred_type == 1
        feature_matrix = [bc';cc';ec';dc';pc'];
    else
        feature_matrix = [cc';ec';dc';pc'];
    end
    csvwrite(feature_file, feature_matrix);
end
   
    
function [CM] = ENM_Dynamic_Correlation(IND)
    pow = 2; % we will vary the spring strength between the residues by (1/dij)^pow
    num_modes = 26; % We will select the first 20 low frequency modes (7:26)
    [V,E] = pfANM(IND,pow);
    [~,CM]=calc_ANISO_BF(V,E,num_modes);
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


