%___________________________________________________________________________________%
%  MTDE: An effective multi-trial vector-based differential evolution algorithm     %
%  and its applications for engineering design problems                             %
%  source codes version 1.0                                                         %
%                                                                                   %
%  Developed in MATLAB R2018a                                                       %
%                                                                                   %
%  Author and programmer: M. H.Nadimi-Shahraki, S. Taghian, S. Mirjalili, H. Faris  %
%                                                                                   %
%  e-Mail:nadimi@ieee.org, shokooh.taghian94@gmail.com,                             %                                                                                          
%  ali.mirjalili@gmail.com, hossam.faris@ju.edu.jo                                  %
%                                                                                   %
%  Homepage:http://www.alimirjalili.com                                             %
%                                                                                   %
%  Main paper: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili, H. Faris            %
%  MTDE: An effective multi-trial vector-based differential evolution               %
%  algorithm and its applications for engineering design problems,                  %
%  Applied Soft Computing Journal 97 (2020) 106761                                  %
%                                                                                   %
%  DOI: 10.1016/j.asoc.2020.106761                                                  %
%___________________________________________________________________________________%

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_Functions_details(F)

switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=30;
end
end

% F1

function o = F1(x)
o=sum(x.^2);
end

% F2

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end

