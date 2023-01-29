% I Grey Wolf Optimizer
function [Alpha_score,Alpha_pos,Convergence_curve]=IGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=Logstic_initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter

% Main loop
MIN_SCORE_ALL = inf;
while l<Max_iter
    f_average = 0.0;
    f_min = inf;
    f_max = -inf;
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        f_average = f_average + fitness;
        if f_min > fitness
            f_min = fitness;
        end
        if f_max < fitness
            f_max = fitness;
        end
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
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0


%     [Alpha_pos, Alpha_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Alpha_score, Alpha_pos);
%     [Beta_pos, Beta_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Beta_score, Beta_pos);
%     [Delta_pos, Delta_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Delta_score, Delta_pos);

    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1) % wolves
        for j=1:size(Positions,2) % dim
                       
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
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            %Positions(i,j) = (X1 * Alpha_score + X2 * Beta_score + X3 * Delta_score) / (Alpha_score + Beta_score + Delta_score);
        end
    end

 

    l=l+1;
    if (Alpha_score < MIN_SCORE_ALL)
        MIN_SCORE_ALL = Alpha_score;
    end
    Convergence_curve(l)=MIN_SCORE_ALL;
end