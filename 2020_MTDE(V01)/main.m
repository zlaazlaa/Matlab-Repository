%___________________________________________________________________________________%
%  MTDE: An effective multi-trial vector-based differential evolution algorithm     %
%  and its applications for engineering design problems                             %
%  source codes version 1.0                                                         %
%                                                                                   %
%  Developed in MATLAB R2018a                                                       %
%                                                                                   %
%  Author and programmer: M. H.Nadimi-Shahraki, S. Taghian, S. Mirjalili, H. Faris  %
%                                                                                   %
%  e-Mail:nadimi@ieee.org, shokooh.taghian94@gmail.com,                             %                                                                                          %
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

close all
clear
clc

Algorithm_Name = 'MTDE';
Function_name='F1';                          % Name of the test function that can be test by F1 or F2 

% Load details of the selected benchmark function
[lb, ub, D, fobj] = Get_Functions_details(Function_name);
N = 100;                                     % The population size
MaxFES = D * 10000;

disp(['MTDE is running for Func ', Function_name]);
[gbestval,gbestpos,Convergence_curve] = MTDE(D,N,MaxFES,lb,ub,fobj);
disp(['The best solution obtained by MTDE is : ', num2str(gbestpos)]);
disp(['The best optimal value of the objective funciton found by MTDE is : ', num2str(gbestval)]);
