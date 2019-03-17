function [] = GetSequenceHydropathy(pdb_ca_file, featurefile)
    % Based on the Kyte and Doolittle scale of amino acid Hydrophobicity,
    % parse through the sequence of the given PDB C-alpha file assigning the
    % hydrophobic values and writing into the corresponding feature file.
    KD_scale = struct('ILE', 4.5, 'VAL', 4.2, 'LEU', 3.8, 'PHE', 2.8, 'CYS', 2.5, 'MET', 1.9, 'ALA', 1.8, 'GLY', -0.4, 'THR', -0.7, 'TRP', -0.9, 'SER', -0.8, 'TYR', -1.3, 'PRO', -1.6, 'HIS', -3.2, 'GLU', -3.5, 'GLN', -3.5, 'ASP', -3.5, 'ASN', -3.5, 'LYS', -3.9, 'ARG', -4.5);
	
    P = readPDB(pdb_ca_file,1);
    resnames = P.resName;
    H_index = [];
    for j=1:length(resnames)
	    H_index = [H_index;KD_scale.(resnames{j})];
    end
        
    csvwrite(featurefile, H_index');
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

