function [] = CalculateShortestPathToActiveSiteRes(pdb_ca_file, feature_file, activesite_pdb_ind)
    % For a residue, calculate the shortest path to the active site 
    % residues. We will consider the shortest path to any one of the
    % residues in the active site from the graph of correlated
    % fluctuations expressed as a distance matrix.
    
    dist_cutoff = 13; % Distance cutoff for long range interactions
    PDB_struct = readPDB(pdb_ca_file, 1); 
    IND = PDB_struct.IND;
    
    nres = size(IND,1);
    
    % If no activesite residues are provided then there is nothing
    % to calculate. Assign 0s to all feature values.
    if ~exist('activesite_pdb_ind', 'var') || isempty(activesite_pdb_ind)
        shortest_path_vect = zeros(1,nres);
        median_all_shortest_paths = zeros(1,nres);
    else
        CM = ENM_Dynamic_Correlation(IND);
        CM = CM ./ max(abs(CM(:))); % normalize
        CM_dist = 1 - CM; % express as a distance matrix %%%% We need to also consider the distance
        CM_dist_adjusted = DistCutoffandCorrWeightedMat(IND,dist_cutoff,CM_dist);
        G = graph(CM_dist_adjusted);
        shortest_path_vect = zeros(nres,1); % Store the shortest distance for each residue to
        % on of the catalytic res.
        median_all_shortest_paths = zeros(nres,1); % Store the median of the shortest paths between a residue and all catalytic residues
        % Get the serial indices for the active site residues
        activesite_serial_ind = FindSerialIndexForActiveSite(PDB_struct, activesite_pdb_ind);
        
        % For each residue find the shortest path to one of the catalytic
        % residues
        for r = 1:nres
            shortest_path_r = [];
            all_shortest_paths = []; % store all the shortest paths between r and all catalytic residues
            for c = (activesite_serial_ind')
                path_nodes = G.shortestpath(r,c); % Get the nodes for the shortest path between r and c
                path_len = numel(path_nodes)-1; % Get the path length excluding the start node
                if isempty(shortest_path_r)
                    shortest_path_r = path_len;
                elseif path_len<shortest_path_r
                    shortest_path_r = path_len;
                end
                all_shortest_paths = [all_shortest_paths;path_len];
            end
            shortest_path_vect(r)=shortest_path_r; % include the shortest distance to one of the catalytic residues
            median_all_shortest_paths(r) = median(all_shortest_paths(:)); % include the median of all shortest paths from a given residue to a set of catalytic residues
        end
    end 

    % Write the shortest path vector to the feature file
    feature_matrix = [shortest_path_vect';median_all_shortest_paths'];
    csvwrite(feature_file, feature_matrix);
end

function [serial_index] = FindSerialIndexForActiveSite(PDB_struct, activesite_res_ind)
    % Convert the PDB index into serial index for the active site residues
    serial_index = [];
    PDB_index = PDB_struct.resSeq;
    for m = (activesite_res_ind)
        s = find(PDB_index == m);
        serial_index = [serial_index; s];
    end    
end

function [CM] = ENM_Dynamic_Correlation(IND)
    pow = 2; % we will vary the spring strength between the residues by (1/dij)^pow
    num_modes = 26; % We will select the first 20 low frequency modes (7:26)
    [V,E] = pfANM(IND,pow);
    [~,CM]=calc_ANISO_BF(V,E,num_modes);
end


function [CM_dist_adjusted] = DistCutoffandCorrWeightedMat(IND, dist_cutoff,CM_dist)
    % Use the distance cutoff to identify interacting residue pairs and
    % then only include the distance transformed correlations for these
    % pairs.
    nres = size(IND,1);
    CM_dist_adjusted = zeros(nres,nres);
    for i=1:nres-1
        for j=i+1:nres
            dist = norm(IND(i,:)-IND(j,:));
            if dist <= dist_cutoff
                CM_dist_adjusted(i,j) = CM_dist(i,j);
                CM_dist_adjusted(j,i) = CM_dist(j,i);
            else
                % If the distance is more than the distance cutoff of 13
                % Ang, we do not want to take this path. So assign a high
                % edge cost/weight
                CM_dist_adjusted(i,j) = 10000;
                CM_dist_adjusted(j,i) = 10000;
            end    
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

