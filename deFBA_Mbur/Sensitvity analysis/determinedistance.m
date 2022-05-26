function d = determinedistance(r1,r2)
% This function determines the distance between two Morris trajectories
% INPUTS:
% - r1: trajectory 1
% - r2: trajectory 2
% OUTPUTS:
% - d: the distance between the two trajectories r1 and r2

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

k = size(r1,2);

d2 = 0;
d4 = 0;
d = 0;
for i = 1:k+1
    for j = 1:k+1
        for z = 1:k
            d1 = (r1(i,z)-r2(j,z))^2;
            d2 = d2 + d1;
        end
        d3 = sqrt(d2);
        d4 = d4 + d3;
    end
    d = d + d4;
end

end