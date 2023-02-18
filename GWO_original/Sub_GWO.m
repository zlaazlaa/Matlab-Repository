function [Alpha_score,Positions] = Sub_GWO(start_i,end_i,Positions,fobj,ub,lb,dim, Max_iter,l)
    % initialize alpha, beta, and delta_pos
    Alpha_pos=zeros(1,dim);
    Alpha_score=inf; %change this to -inf for maximization problems
    
    Beta_pos=zeros(1,dim);
    Beta_score=inf; %change this to -inf for maximization problems
    
    Delta_pos=zeros(1,dim);
    Delta_score=inf; %change this to -inf for maximization problems

    for i=start_i:end_i
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
%         if fitness < best_score_of_every
%             best_score_of_every = fitness;
%             best_of_every(i,:) = Positions(i,:);
%         end
%         if fitness < best_score_of_whole
%             best_score_of_whole = fitness;
%             best_of_whole = Positions(i,:);
%         end

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
    
    % Update the Position of search agents including omegas
    for i=start_i:end_i
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
            
        end
    end

    %% Myxomycetes information exchange
    old_position = Positions;
    for i = start_i:end_i % position to be changed
        %for j = 1:SearchAgents_no % neighbor selected
        j = randi([start_i,end_i]);
        diff = old_position(i) - old_position(j);
        trust_score = abs(fobj(old_position(i,:)) - fobj(old_position(j,:))) / (norm(diff,2)+0.01);
        if (fobj(old_position(i,:)) - fobj(old_position(j,:)) > 0)
            Positions(i,:) = Positions(i,:) - 0.1 * (1 - exp(-trust_score)) * diff;
        else
            Positions(i,:) = Positions(i,:) + 0.1 * (1 - exp(-trust_score)) * diff;
        end
        %end
        %Positions(i,:) = Positions(i,:) + 0.5 * rand() * (best_of_whole - Positions(i,:)) + 0.5 * rand() * (best_of_every(i,:) - Positions(i,:));
        for j = 1:dim % prevent error solution
            if Positions(i,j) > ub
                Positions(i,j) = ub;
            end
            if Positions(i,j) < lb
                Positions(i,j) = lb;
            end
        end
        if fobj(Positions(i,:)) > fobj(old_position(i,:))
            Positions(i,:) = old_position(i,:);
        end
    end
end