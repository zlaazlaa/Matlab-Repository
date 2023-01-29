%___________________________________________________________________%
%  Grey Wolf Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% This function initialize the first population of search agents
function Positions=Logstic_initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    Z = zeros(SearchAgents_no, dim);
    
    % 随机生成一个d维向量
    Z(1, :) = rand(1, dim);
    
    % 利用logistic生成n_pop个向量
    for i=2:SearchAgents_no
        Z(i,:) = 4.0*Z(i-1,:).*(1-Z(i-1,:));
    end
    
    % 将z的各个分量载波到对应变量的取值区间
    Positions = zeros(SearchAgents_no, dim);
    for i=1:SearchAgents_no
        Positions(i,:) = lb + (ub - lb)*Z(i,:);
    end
end