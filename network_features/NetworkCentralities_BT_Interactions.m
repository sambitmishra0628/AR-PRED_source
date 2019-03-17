function [] = NetworkCentralities_BT_Interactions(pdb_ca_file, feature_file, pred_type)
    % Calculate the centralities for each node by simulating each protein
    % as a weighted network, with the edge weights corresponding to the
    % pairwise energies between amino acid pairs obtained from the
    % Betancourt and Thirumalai potential.
    % We calculate the following centralities.
    % 1. Betweenness
    % 2. Closeness
    % 3. Eigen vector
    % 4. Degree
    % 5. Page rank
	
	% The predictor type here corresponds to the predictor type.
	% 	1: Active site,  
	%	2: Allosteric site,
	
    pot_file = strcat('BT-full.csv');
    BT_POT = csvread(pot_file);
    RES_ORDER = {'CYS','PHE','LEU','TRP','VAL','ILE','MET','HIS','TYR','ALA','GLY','PRO','ASN','THR','SER','ARG','GLN','ASP','LYS','GLU'};

    P = readPDB(pdb_ca_file,1);
    resnames = P.resName;
    nres = numel(resnames);
    BT_pdb = zeros(nres,nres); % To store the network of contact energies for the current pdb
    for m=1:nres-1
        for n=m+1:nres
            resname_m = resnames{m};
            resname_n = resnames{n};
            ind_m = find(ismember(RES_ORDER,resname_m));
            ind_n = find(ismember(RES_ORDER,resname_n));
            BT_mn = BT_POT(ind_m,ind_n); % Find the contact energy between residue m and residue n

            % convert to positive by expressing as boltzmann factors.
            % Note that in the matrix, the favourable interactions are
            % more negative (less energy and hence, energetically
            % favourable).
            BT_mn = exp(-BT_mn);
            BT_pdb(m,n) = BT_mn;
            BT_pdb(n,m) = BT_mn;
        end
    end
   
    % Perform centrality calculatations on the weighted potential
    % network BT_pdb_i. For computing betweenness and closeness, we
    % will use edge weights as 1/BT(i,j) so that the
    % chosen shortest path (path with minimum weight) has the most favourable
    % interaction energy.
    G = graph(BT_pdb);

	% Calculate the betweenness centralities
	bc = G.centrality('betweenness','Cost',1./G.Edges.Weight); % treat as weighted graph in case corr cutoff is not specified

	% Normalize the betweenness values
	n_nodes = numnodes(G);
	norm_factor = ((n_nodes-2)*(n_nodes-1))/2; % normalize to make the bc as probabilities
	bc = bc ./ norm_factor;
 
    % Compute the closeness centrality
    cc = G.centrality('closeness','Cost',1./G.Edges.Weight);
    
    % Compute the eigen centrality. For eigen, degree and page rank centralities
    % we assign edge 'Importance' equal to the interaction energy.
    ec = G.centrality('eigenvector','Importance',G.Edges.Weight);
    
    dc = G.centrality('degree','Importance',G.Edges.Weight);
    
    pc = G.centrality('pagerank','Importance',G.Edges.Weight);
    
    % Write the feature vectors into the respective feature
    % files depending on the predictor type.
    if pred_type == 1
        feature_matrix = [bc';cc';ec';dc';pc'];
    else
        feature_matrix = [cc';ec';dc';pc'];
    end	
    csvwrite(feature_file, feature_matrix);
    
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
