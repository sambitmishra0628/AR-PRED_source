function [] = DFIandAlloResponse(pdb_ca_file,feature_file,activesite_res_ind)
    % For the set of c-alpha PDB files in the pdb_ca_dir, perform
    % perturbation response analysis and add the DFI and Allosteric
    % response scores.
    
    pow = 2; % For the distance dependent ENM, we will use an exponent of 2.
    iterations = 10; % Number of iterations for which we will repeat the
					 % calculations for DFI with different randomly selected points.
    
    pdb = readPDB(pdb_ca_file,1); % get all elements of the PDB file
    COORD = pdb.IND;
    activesite_res_ind = FindSerialIndexForActiveSite(pdb, activesite_res_ind); % Get the serial index of the active site residues   
    
    % get the residue names and total number for all residues in the structure
    %resnames = pdb.resName;
    total_res = size(COORD,1);
    %nres_sc = numel(find(I==1)); % number of residues in a single chain
    %num_perturb_res = numel(res_num); % The total number of residues to be perturbed
    
    % Compute the DFI profile of the protein upon perturbing all residues
    fluct_all = zeros(total_res,1);
    dR_all = zeros(total_res, total_res);
    dfi_all = zeros(total_res,1);
    
    % We will perform calculations for Hessian and Hessian inverse here
    % and not each time inside an iteration.
    [V,E] = pfANM(COORD,pow);
    rn = size(COORD,1); % number of residues
    
    % predefine the positional displacement matrix delta R. This matrix
    % will contain the displacements of a residue caused upon perturbing
    % all other residues, one at a time (the individual col vectors.
    
    n = size(V,1); % get the number of rows in the eigen vector. The inverse Hessian
    % matrix will have the same number of rows and columns
    num1 = 7;num2 = size(V,1);
    Hinv = sparse(n,n);
    % Understand that we are not computing the true inverse of the Hessian.
    % Rather, we are computing the pseudo inverse, because the first 6
    % eigen values of the Hessian matrix are zero and we have to discard
    % the first 6 eigen values and vectors
    for a=num1:num2
        Hinv = Hinv + (1/E(a))*(V(:,a)*V(:,a)');
    end
    % We want to randomize the perturbation process to eliminate bias from force magnitude or direction.
    % Hence, multiple iterations.
    for iter=1:iterations
        % Perturb all residues and get the response for individual
        % residues.
        [fluct,dR,dfi] = DynamicIndex(COORD,Hinv);
        fluct_all = fluct_all + fluct;
        dR_all = dR_all + dR;
        dfi_all = dfi_all + dfi;
    end
    % Get the average fluctuation and response for a single iteration
    % (Normalized by number of iterations)
    fluct_all = fluct_all ./ iterations;
    dR_all = dR_all ./ iterations;
    dfi_all = dfi_all ./ iterations;
    fluct_all_normalized = fluct_all ./ (total_res); % Normalize by (total res). Result is a column vector
    feature_matrix = [dfi_all'];
    
    % Now calculate the perturbation response on active site if residues of the active 
	% site are provided.
	if exist('activesite_res_ind', 'var') && ~isempty(activesite_res_ind)
    	num_active_site = numel(activesite_res_ind);
    	dR_activesite = dR_all(:,activesite_res_ind); % effect of perturbing activesite residues on remaining residues
	    ACTIVESITE_PERTURBATION_RESPONSE = (sum(dR_activesite,2)./num_active_site)./fluct_all_normalized; % Perturbation response of active site on remaining residues
        feature_matrix = [feature_matrix;ACTIVESITE_PERTURBATION_RESPONSE'];
 	end

    % Write the features (DFI and active site perturbations response) into the
    % specific feature file.
    csvwrite(feature_file, feature_matrix);
end


function [S, dR, dfi] = DynamicIndex(COORD,Hinv)
    % Given the coordinates of a PDB structure, compute the Dynamic
    % Flexibility Index matrix (refer
    % http://www.ncbi.nlm.nih.gov/pubmed/23745135)
    % COORD should be N by 3 matrix of residue coordinates where N is 
    % the number of residues
    rn = size(COORD,1); % number of residues
    % predefine the positional displacement matrix delta R. This matrix
    % will contain the displacements of a residue caused upon perturbing
    % all other residues, one at a time (the individual col vectors.
    
    dR = zeros(rn,rn); % to hold the fluctuation response of all the residues upon
                       % perturbing the residues in res_num
    n = size(Hinv,1); 
    for i=1:rn
        COORD_i = COORD(i,:); % Coordinates of the ith residue
        F_full = zeros(n,1); % The size of the force vector should be 3N by 1 
        
        s = 1;        
        % Generate 100 random points uniformly (sort of) around the ith
        % residue. We will then have 100 force vectors aligned along the vector
        % originating from each of these random points and ending at the
        % ith residue.
        randpoints = bsxfun(@plus,COORD_i,s.*randn(100,3)); 
        %plot3(randpoints(:,1),randpoints(:,2),randpoints(:,3),'.r')
        %plot3(COORD_i(1),COORD_i(2),COORD_i(3),'*g');
        %hold on;
        
        % Generate 100 forces with random magnitude ranging from 1 to 5
        % which, we will later use to perturb the given residue i. We restrict ourselves
        % only to a maximum magnitude of 5 to ensure that the response
        % upon perturbing each point is not highly biased by the magnitude
        % of force. 
        F_rand = randi(5,100,1);
        
        response = zeros(n,1); % The response vector for residue i
       % Generate response vector for each force vector. The response
       % vector will be a 3N by 1 vector having the positional
       % displacements for each residue in x,y and z directions upon
       % application of force on residue i
        for j=1:100
            v_pos = COORD_i-randpoints(j,:); % The position vector starting at a random point and ending at residue i
            v_pos = v_pos./norm(v_pos); % unit position vector
            
            % We will now get the angles for each component of the position
            % vector using the concept v_x = v cos theta, where v is the
            % magnitude of the vector. Note that the above formula comes
            % from the definition of dot product: a.b = |a||b|cos theta. If
            % b is the unit vector along x, then b = [1 0 0]. So the angle
            % between a and b i.e vector and x-axis would be arccos(a.b/|a||b|)
            theta_x = acosd(v_pos(1)/norm(v_pos)); % the angle with X-axis in degree
            theta_y = acosd(v_pos(2)/norm(v_pos)); % the angle with Y-axis in degree
            theta_z = acosd(v_pos(3)/norm(v_pos)); % the angle with Z-axis in degree
            
            Fj = F_rand(j); % The jth force 
            
            % Find the respective components of the force Fj. These
            % resultant force vector should be parallel to v_pos and act on
            % residue i.
            Fj_x = Fj*cosd(theta_x);
            Fj_y = Fj*cosd(theta_y);
            Fj_z = Fj*cosd(theta_z);
            
            % Update the full force vector with the individual force values
            % i.e. fx, fy and fz for residue i(x,y and z components). All
            % other residues will have force as 0.
            F_full((3*(i-1))+1) = Fj_x; F_full((3*(i-1))+2) = Fj_y; F_full((3*(i-1))+3) = Fj_z; 
            r = Hinv*F_full;
            response = response + r;
        end
        % Get the average response over the 100 force vectors
        response = response ./ 100;
        % reshape the response vector as a N by 3 vector
        response = reshape(response,3,rn)'; 
        % Get the magnitude of response for each residue upon perturbing
        % residue i as sqrt(x^2+y^2+z^2)
        response = sqrt(sum(response.^2,2));
        
        % add the computed response vector as a column vector in matrix dR
        % Here, dR is the perturbation response matrix in which each column
        % is the effect of perturbing a given residue on all other
        % residues.
        dR(:,i) = response; 
    end 
    % Calculate the average displacements (mean square fluctuations) for
    % all residues by only considering rows of perturbed residues
    % Sum the individual row vectors in the response matrix to get the
    % average fluctuation of each residue upon perturbing all the residues.
    S = sum(dR,2);
                                        
%   Compute dfi
    dfi = S./sum(S); 
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
