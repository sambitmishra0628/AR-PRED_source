    function [] = NetworkCentralities_AdjMat(pdb_ca_file, feature_file)
  
    % Calculate network properties for the AA interaction
    % networks created based on the dist_cutoff from the CA PDB files.
    % We will create adjacency matrices using the distance cutoffs and
    % then calculate the following properties.
    % 1. Betweenness centrality
    % 2. Eigen centrality
    % 3. Closeness centrality
    % 4. Degree centrality
    % 5. Page rank centrality
    
	
	dist_cutoff = 13; % The optimal distance cutoff to model the protein as a network
    coord = readPDB(pdb_ca_file);
    adj_mat = MakeAdjacencyMatrix(coord, dist_cutoff);
    G = graph(adj_mat); % Create a graph object using the matlab graph module
            
    % Perform betweenness centrality calculations for the adjacency matrix
    bc = G.centrality('betweenness');
    % Convert the bc values from numbers to probability scores by
    % dividing by a factor of (n-2)(n-1)/2
    % (https://www.mathworks.com/help/matlab/ref/graph.centrality.html)
    n_nodes = numnodes(G);
    norm_factor = ((n_nodes-2)*(n_nodes-1))/2;
    bc = bc ./ norm_factor;
    
    
    % Perform closeness centrality calculations on the adjacency matrix
    cc = G.centrality('closeness');

    % Perform eigen centrality calculations on the adjacency matrix
    ec = G.centrality('eigenvector');
        
    % Perform degree centrality calculations on the adjacency
    % matrix
    dc = G.centrality('degree');
    
    
    % Perform page rank centrality calculations on the adjacency
    % matrix
    pc = G.centrality('pagerank');
       
    feature_matrix = [bc';cc';ec';dc';pc'];
    csvwrite(feature_file, feature_matrix);

end  
    
function [adj_mat] = MakeAdjacencyMatrix(coord, dist_cutoff)
    nres = size(coord,1);
    adj_mat = zeros(nres,nres);
    for i=1:nres-1
        for j=i+1:nres
            dist = norm(coord(i,:)-coord(j,:));
            if dist<=dist_cutoff
                adj_mat(i,j) = 1;
                adj_mat(j,i) = 1;
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
