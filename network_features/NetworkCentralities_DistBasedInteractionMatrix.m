function [] = NetworkCentralities_DistBasedInteractionMatrix(pdb_ca_file, feature_file)
    % Compute the following network properties by using the matrix of
    % distance between a residue pair which is equivalent the inverse of
    % the interaction strength. When using the distance matrix, we are
    % actually stating that the cost of the edge between two nodes is the
    % distance of separation. This means that when getting the shortest
    % path for node betweenness and closeness, the algorithm will consider
    % the path with least weight. However, while calculating degree, page
    % rank and eigen centralities we will use the inverse of the distance
    % of separation to weight the importance of edges and perform
    % calculations.
    % 1. Eigen centrality
    % 2. Closeness centrality
    % 3. Degree centrality
    % 4. Page rank centrality
    
    coord = readPDB(pdb_ca_file);
    dist_mat = MakeDistMatrix(coord); % Get the distance matrix from the CA coordinates

    % Get the raw centrality values. We will use the inbuilt MATLAB
    % 2017 graph functions to calculate the betweenness, closeness,
    % degree, page rank, eigen centralities and vertex eccentricity
    % for the weighted network.
    
    % Create a graph object from the matrix
    G = graph(dist_mat);
    
    % Compute the closeness centrality
    cc = G.centrality('closeness','Cost',G.Edges.Weight);
        
    % Compute the eigen centrality. For eigen centrality, we require
    % the importance of each edge. We will make the importance of each
    % edge proportionate to the strength of interaction which is 1/dij.
    ec = G.centrality('eigenvector','Importance',1./G.Edges.Weight);
    
    % Compute the degree centrality by using inverse distance of
    % separation as the edge weight.
    dc = G.centrality('degree','Importance',1./G.Edges.Weight);
        
    % Compute the page rank centrality by using inverse distance of
    % separation as the edge weight.
    pc = G.centrality('pagerank','Importance',1./G.Edges.Weight);
            
    % Write the calculated feature vectors into the respective feature
    % files.
    feature_matrix = [cc';ec';dc';pc'];
    csvwrite(feature_file, feature_matrix);
end

    
function [dist_mat] = MakeDistMatrix(coord)
    nres = size(coord,1);
    dist_mat = zeros(nres,nres);
    for i=1:nres-1
        for j=i+1:nres
            dist = norm(coord(i,:)-coord(j,:));
            dist_mat(i,j) = dist;
            dist_mat(j,i) = dist;
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

