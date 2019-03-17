function [pdbindex] = GetResiduePDBIndex(pdb_ca_file)
    % Get the residue PDB indices for a given pdb c-alpha file
    currentdir = pwd();
    fprintf('Getting residue PDB index for file %s...',pdb_ca_file);
    P = readPDB(strcat(currentdir,'/',pdbfile),1);
    pdbindex = P.resSeq;
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

