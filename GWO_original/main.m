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

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run GWO: [Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________

clear all
for kw = 1:22
    if kw == 17
        continue
    end
    SearchAgents_no=100; % Number of search agents
    Function_name=strcat('F', num2str(kw)); % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

    Max_iteration=500; % Maximum numbef of iterations
    
    % Load details of the selected benchmark function
    [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
    
    [Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    [my_Best_score,my_Best_pos,my_GWO_cg_curve]=IGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    PSO_cg_curve=PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); % run PSO to compare to results
    [Best_score_SCA,Best_pos_SCA,cg_curve_SCA]=SCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    [my_Best_score2,my_Best_pos2,my_GWO_cg_curve2]=IGWO2(dim,SearchAgents_no,Max_iteration,lb,ub,fobj);
    figure('Position',[500 500 660 290])
    fig = figure(Visible="off");
    % Draw search space
    subplot(1,2,1);
    func_plot(Function_name);
    title('Parameter space')
    xlabel('x_1');
    ylabel('x_2');
    zlabel([Function_name,'( x_1 , x_2 )'])
    
    %Draw objective space
    subplot(1,2,2);
    semilogy(GWO_cg_curve,'Color','r')
    hold on
    semilogy(my_GWO_cg_curve, 'Color', 'g');
    semilogy(PSO_cg_curve, 'Color', 'b');
    semilogy(cg_curve_SCA, 'Color', 'magenta')
    semilogy(my_GWO_cg_curve2, 'Color', 'black');
    title('Objective space')
    xlabel('Iteration');
    ylabel('Best score obtained so far');
    
    axis tight
    grid on
    box on
    legend('GWO', 'myGWO', 'PSO', 'SCA', 'best-GWO')
    asdas = num2str(kw);
    disp(asdas)
    file_name = strcat('C:\Users\zlaa\Documents\MATLAB\GWO_original\all_fun\', num2str(kw), '.png');
    frame = getframe(fig);
    %imwrite(frame.cdata, file_name, 'png', 600);
    print(fig, '-dpng', '-r600', strcat('C:\Users\zlaa\Documents\MATLAB\GWO_original\all_fun\', num2str(kw), '.png'))  
    display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
    display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);
    disp('-------------------------------')
end