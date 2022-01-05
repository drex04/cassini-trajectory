function [S] = stumpffS(z)
% A function to calculate S(z) coefficient for Kepler's Universal Variables
    if z > 0
        S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
    elseif z == 0
        S = 1/6;
	elseif z < 0
        S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z))^3;
    else
        disp('Error');
	end

end