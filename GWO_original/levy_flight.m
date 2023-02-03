function length = levy_flight()
beta = 3/2;
sigma_u = power((gamma(1+beta)*sin(pi*beta/2))/(gamma((1+beta)/2)*beta*pow2((beta-1)/2)),1/beta);
v = normrnd(0,1);
u = normrnd(0,power(sigma_u,2));
length = u/abs(v).^(1/beta);