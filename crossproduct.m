%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: crossproduct of 2 vectors

function crossprod = crossproduct(vect1, vect2)

crossprod = [vect1(2)*vect2(3)-vect1(3)*vect2(2), vect1(1)*vect2(3)-vect1(3)*vect2(1), vect1(1)*vect2(2)-vect1(2)*vect2(1)];