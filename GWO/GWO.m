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

% Grey Wolf Optimizer
function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj,tag)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter

% Main loop
while l<Max_iter
    %strArray(size(Positions,1))=struct('score', 'pos');
    for i=1:size(Positions,1)
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        
        % put data into struct
        strArray(i).score=fitness;
        strArray(i).pos=Positions(i,:);

        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score 
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    if tag == 0
        a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    end
    if tag == 1
        %a=2*(1-(l*l)/(Max_iter*Max_iter));
        if l < Max_iter / 2
             a = sin(pi * 2 * l / Max_iter + pi / 2) + 1;
            %a=2-2*l*((2)/Max_iter);
        else
             a = sin(pi * 2 * (l - Max_iter / 2) / Max_iter + pi / 2) + 1;
            %a=2-2*(l - Max_iter / 2)*((2)/Max_iter);
        end
    end

    % sort the struct
    [~, index] = sort([strArray.score]);
    XSum=zeros(1,dim);
    for i = 1:dim/2
        XSum = XSum+strArray(index(i)).pos;
    end
    XSum = XSum/dim-dim/2+1;

    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            if tag == 1 
                Positions(i,j)=(X1+X2+X3+XSum(j))/4;% Equation (3.7)
            else
                Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            end
        end
    end
    l=l+1;
    disp(Alpha_score);
    xx=zeros(1,size(Positions, 1));
    yy=zeros(1,size(Positions, 1));
    for i=1:size(Positions, 1)
        xx(i)=Positions(i,1);
        yy(i)=Positions(i,2);
    end
    %fig = figure; % 新建一个figure，并将图像句柄保存到fig
    %scatter(xx,yy);
    %axis([lb,ub,lb,ub]);
    %frame = getframe(fig); % 获取frame
    %img = frame2im(frame); % 将frame变换成imwrite函数可以识别的格式
    %aa = 'a';
    %pp = '.png';
    %nn = num2str(l);
    %imwrite(img, [aa,nn,pp]); % 保存到工作目录下，名字为"a.png"
    Convergence_curve(l)=Alpha_score;
end



