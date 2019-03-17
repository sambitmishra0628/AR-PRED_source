function [] = AminoAcidIdentity(pdb_ca_file, feature_file)
    % Given the c-alpha file, get the corresponding amino acid sequence
	% and assign each residue to a numeric value between 1 through 20.
	
    P = readPDB(pdb_ca_file,1);
    resnames = P.resName; % Residue names for jth chain
    numeric_class_j = AAIdentity(resnames);
	csvwrite(feature_file, numeric_class_j');
	return;
end

function [numeric_class] = AAIdentity(resnames)
    numeric_class = [];
    aa_class = struct('ALA',1,'ARG',2,'ASN',3,'ASP',4,'CYS',5,'GLU',6,'GLN',7,'GLY',8,'HIS',9,'ILE',10,'LEU',11,'LYS',12,'MET',13,'PHE',14,'PRO',15,'SER',16,'THR',17,'TRP',18,'TYR',19,'VAL',20);
    for i=1:numel(resnames)
        resi = resnames{i};
        class_i = aa_class.(resi);
        numeric_class = [numeric_class;class_i];
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
