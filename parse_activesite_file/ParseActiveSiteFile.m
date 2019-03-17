function [active_site_ind] = ParseActiveSiteFile(activesitefile)
    % Given the PDB indices of residues which form the active site, parse
    % the active site residue indices into a vector. We expect the file to
    % be a .csv file i.e, having the residue indices separated by comma for the single
    % protein subunit/chain under consideration.
	if exist(activesitefile, 'file')
	    active_site_ind = csvread(activesitefile);
	else
		error ('ERROR! Could not find active site file ', activesitefile, '!');
		exit;
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

