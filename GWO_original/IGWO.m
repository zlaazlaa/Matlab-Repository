% I Grey Wolf Optimizer
function [Alpha_score,Alpha_pos,Convergence_curve]=IGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

Worse_pos=zeros(1,dim);
Worse_score=-inf;

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

        if fitness>Worse_score
            Worse_score=fitness;
            Worse_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0

    range_x = (Max_iter - l)/Max_iter;
    range_y = (1-range_x)/3;

    normal_sum = round(range_x*SearchAgents_no);
    levy_sum = round(range_y*SearchAgents_no);
    brownian_sum = round(range_y*SearchAgents_no);
    de_sum = SearchAgents_no-normal_sum-levy_sum-brownian_sum;

%     [Alpha_pos, Alpha_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Alpha_score, Alpha_pos);
%     [Beta_pos, Beta_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Beta_score, Beta_pos);
%     [Delta_pos, Delta_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Delta_score, Delta_pos);


    % Update the Position of search agents
    for i=1:normal_sum % wolves
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

%     Alpha_pos_new = Alpha_pos + levy_flight()*(Worse_pos-Alpha_pos);
%     Beta_pos_new = Beta_pos + levy_flight()*(Worse_pos-Beta_pos);
%     Delta_pos_new = Delta_pos + levy_flight()*(Worse_pos-Delta_pos);
% 
%     if Alpha_score > fobj(Alpha_pos_new)
%         Alpha_pos = Alpha_pos_new;
%         Alpha_score = fobj(Alpha_pos_new);
%     end
%     if Beta_score > fobj(Beta_pos_new)
%         Beta_pos = Beta_pos_new;
%         Beta_score = fobj(Beta_pos_new);
%     end
%     if Delta_score > fobj(Delta_pos_new)
%         Delta_pos = Delta_pos_new;
%         Delta_score = fobj(Delta_pos_new);
%     end

    % use levy flight to create new wolves
    for i = normal_sum + 1:normal_sum+levy_sum
        rand_num = randi([1,SearchAgents_no]);
        Positions(i,:) = Positions(rand_num,:)+0.1*levy_flight()*(Worse_pos-Positions(rand_num,:));
        for j = 1:size(Positions,2)
            if Positions(i,j) > ub
                Positions(i,j) = ub;
            end
            if Positions(i,j) < lb
                Positions(i,j) = lb;
            end
        end
    end

    % use brownian to create new wolves
    for i = normal_sum+levy_sum + 1:normal_sum+levy_sum+brownian_sum
        Positions(i,:) = brownian(Positions,lb,ub,Alpha_pos);
    end

    % use DE+SA to create new wolves
    for i = normal_sum+levy_sum+brownian_sum+1 : SearchAgents_no
        rand_num = randi([1,SearchAgents_no]);
        Positions(i,:) = update_wolf(l,Max_iter,Positions,dim,fobj,lb,ub,fobj(Positions(rand_num,:)),Positions(rand_num,:));
    end

    if (Alpha_score > fobj(Alpha_pos))
        Alpha_score = fobj(Alpha_pos);
        %disp('yes')
    end

    
    if Beta_score < Alpha_score
        Alpha_score = Beta_score;
    end
    if Delta_score < Alpha_score
        Alpha_score = Delta_score;
    end
    l=l+1;
    if (Alpha_score < MIN_SCORE_ALL)
        MIN_SCORE_ALL = Alpha_score;
    end
    Convergence_curve(l)=MIN_SCORE_ALL;
end