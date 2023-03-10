% Author: Seyedali Mirjalili
% www.alimirjalili.com
% https://scholar.google.com/citations?user=TJHmrREAAAAJ&hl=en

function [ o ] = objectiveFunction(x, testfunctionNo)

switch testfunctionNo
    case 1
        o = 10 * -cos( 2* (x(1)^2+x(2)^2)^0.5 ) * exp( -0.5 * ((x(1)+1)^2+(x(2)-1)^2)^0.5 ) + 5.1 ;
    case 2
        o = -peaks(x(1),x(2)) + 8.5 ;
    case 3
        o =-1*( 0.2 + x(1).^2 + x(2).^2 - 0.1*cos(6*pi*x(1)) - 0.1*cos(6*pi*x(2)) );
    case 4
        o = -1 * sum(x.^2);
end

end