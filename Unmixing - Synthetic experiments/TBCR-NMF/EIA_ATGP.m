% // ====================================================================
% // This file is part of the Endmember Induction Algorithms Toolbox for MATLAB 
% // Copyright (C) Grupo de Inteligencia Computacional, Universidad del 
% // País Vasco (UPV/EHU), Spain, released under the terms of the GNU 
% // General Public License.
% //
% // Endmember Induction Algorithms Toolbox is free software: you can redistribute 
% // it and/or modify it under the terms of the GNU General Public License 
% // as published by the Free Software Foundation, either version 3 of the 
% // License, or (at your option) any later version.
% //
% // Endmember Induction Algorithms Toolbox is distributed in the hope that it will
% // be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
% // of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% // General Public License for more details.
% //
% // You should have received a copy of the GNU General Public License
% // along with Endmember Induction Algorithms Toolbox. 
% // If not, see <http://www.gnu.org/licenses/>.
% // ====================================================================

function [E,C] = EIA_ATGP(data,p)
%% [E,C] = EIA_ATGP(data,p)
% 
% Manuel Grana <manuel.grana[AT]ehu.es>
% Miguel Angel Veganzones <miguelangel.veganzones[AT]ehu.es>
% Grupo de Inteligencia Computacional (GIC), Universidad del Pais Vasco /
% Euskal Herriko Unibertsitatea (UPV/EHU)
% http://www.ehu.es/computationalintelligence
% 
% Copyright (2011) Grupo de Inteligencia Computacional @ Universidad del Pais Vasco, Spain.
% Copyright (2007) GRNPS group @ University of Extremadura, Spain.
%
% ATGP endmembers induction algorithm.
% ------------------------------------------------------------------------------
% Input:   data      : column data matrix [nvariables x nsamples]
%          p         : number of endmembers to be induced. If not provided it is calculated by HFC method with tol=10^(-5).
%
% Output:  E         : set of induced endmembers [nvariables x p]
%          C         : induced endmembers indexes vector [nsamples] with {0,1} values, where '1' indicates that the corresponding sample has been identified as an endmember.
%
% Bibliographical references:
% [1] A. Plaza y C.-I. Chang, “Impact of Initialization on Design of Endmember Extraction Algorithms�?, Geoscience and Remote Sensing, IEEE Transactions on, vol. 44, nº. 11, págs. 3397-3407, 2006.

%% Parameters
if (nargin < 1)
    error('Insufficient parameters');
end
if (nargin < 2 || p <= 0)
    p = EIA_HFC(data,10^(-5));
end

%% data size
[nvariables,nsamples] = size(data);

%% Algorithm initialization
% the sample with max energy is selected as the initial endmember
max = -1;
idx = 1;
for i = 1:nsamples
    r = data(:,i);
    val = r'*r;
    if val > max
        max = val;
        idx = i;
    end
end
e_0 = data(:,idx);
% Initialization of the set of endmembers and the endmembers index vector
E = zeros(nvariables,p);
E(:,1) = e_0;
C = zeros(1,nsamples);
C(idx) = 1;
% Generate the identity matrix.
I = eye(nvariables);

%% Algorithm
for i = 1:p-1
    UC = E(:,1:i);
    % Calculate the orthogonal projection with respect to the pixels at present chosen. 
    % This part can be replaced with any other distance
    PU = I-UC*pinv(UC'*UC)*UC';
    max = -1;
    idx = 0;
    % Calculate the most different pixel from the already selected ones according to
    % the orthogonal projection (or any other distance selected)
    for j = 1:nsamples
        r = data(:,j);
        result = PU*r;
        val = result'*result;
        if (val > max) 
            max = val;
            idx = j;
        end
    end
    % The next chosen pixel is the most different from the already chosen ones
    e_i = data(:,idx);
    E(:,i+1) = e_i;
    C(idx) = 1;
end
