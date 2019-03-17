function [] = GetAAType(pdb_ca_file, featurefile)
    % Parse the given C-alpha PDB file and obtain the residue sequence.
	% Then, assign each residue to its respective category based on its
	% physico-chemical property.

	P= readPDB(pdb_ca_file,1);
    resnames = P.resName; % Get the residue names
    restypes = AAType(resnames); % Map the residue names to the specific residue type
                        
    % Write the residue type into the feature file 
	csvwrite(featurefile, restypes');
    
end

function [restypes] = AAType(resnames)
    Type1 = {'HIS','ARG','LYS','GLU','ASP'}; % charged residues
    Type2 = {'GLN','THR','SER','ASN','CYS','TYR','TRP'}; % polar residues
    Type3 = {'GLY','PHE','LEU','MET','ALA','ILE','PRO','VAL'}; % hydrophobic residues
    restypes = [];
    for i=1:numel(resnames)
        resi = resnames{i};
        if find(ismember(Type1,resi))
            restypes = [restypes;1];
        elseif find(ismember(Type2,resi))
            restypes = [restypes;2];
        elseif find(ismember(Type3,resi))    
            restypes = [restypes;3];
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

