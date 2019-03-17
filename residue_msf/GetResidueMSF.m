function [] = GetResidueMSF(pdb_ca_file, feature_file)
    % For the given C-alpha PDB file calculate the residue level
    % mean square fluctuations 

    pow = 2; % we will use the pfANM with exponent of 2
    
    IND = readPDB(pdb_ca_file);
    [V,E] = pfANM(IND, pow);
    [~,msf] = calc_ISO_BF(V,E);

    % normalize to unit vector
    msf = msf ./ max(msf(:));
    
    % write the msf into the feature file
    csvwrite(feature_file, msf');
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

