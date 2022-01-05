function [C] = stumpffC(z)
% A function to calculate C(z) coefficient for Kepler's Universal Variables
    if z > 0
        C = (1-cos(sqrt(z)))/z;
    elseif z == 0
        C = 1/2;
	elseif z < 0
        C = (cosh(sqrt(-z))-1)/-z;
    else
        disp('Error');
    end
end