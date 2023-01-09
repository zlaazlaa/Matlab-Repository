function vi = boundConstraint (vi, lu)

% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang, jingqiao@gmail.com

[ps, D] = size(vi);  % the population size and the problem's dimension

vi=((vi>=lu(1, :))&(vi<=lu(2, :))).*vi...
        +(vi<lu(1, :)).*(lu(1, :)+0.25.*(lu(2, :)-lu(1, :)).*rand(ps,D))+(vi>lu(2, :)).*(lu(2, :)-0.25.*(lu(2, :)-lu(1, :)).*rand(ps,D));