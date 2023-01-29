%___________________________________________________________________%
%  An Improved Grey Wolf Optimizer for Solving Engineering          %
%  Problems (I-GWO) source codes version 1.0                        %
%                                                                   %
%  Developed in MATLAB R2018a                                       %
%                                                                   %
%  Author and programmer: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili                                                %
%                                                                   %
%         e-Mail:nadimi@ieee.org, shokooh.taghian94@gmail.com,                                         ali.mirjalili@gmail.com                                                   %
%                                                                   %
%                                                                   %
%       Homepage:http://www.alimirjalili.com                                                   %
%                                                                   %
%   Main paper: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili     %
%               An Improved Grey Wolf Optimizer for Solving         %
%               Engineering Problems , Expert Systems with          %
%               Applicationsins, in press,                          %
%               DOI: 10.1016/j.eswa.2020.113917                     %
%___________________________________________________________________%
%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
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

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% N = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run I-GWO: [Best_score,Best_pos,GWO_cg_curve]=IGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________

close all
clear
clc

Algorithm_Name = 'I-GWO';

N = 30; % Number of search agents
Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration = 500; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Fbest,Lbest,Convergence_curve]=IGWO(dim,N,Max_iteration,lb,ub,fobj);
display(['The best solution obtained by I-GWO is : ', num2str(Lbest)]);
display(['The best optimal value of the objective funciton found by I-GWO is : ', num2str(Fbest)]);

figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(Convergence_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('I-GWO')


