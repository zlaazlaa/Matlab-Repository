function [Alpha_score,Alpha_pos,Convergence_curve]=GWO_NianJun(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

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
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        
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
            
        end
    end

    %% Myxomycetes information exchange
    old_position = Positions;
    accept=0;
    reject=0;
    llb = 0.2*(Max_iter - l)/Max_iter;
    uub = 1-llb;
    for i = 1:SearchAgents_no % position to be changed
        for w = 1:round(SearchAgents_no/50) % neighbor selected
            j = randi([1,SearchAgents_no]);
            diff = old_position(i,:) - old_position(j,:);
            i_score = fobj(old_position(i,:));
            j_score = fobj(old_position(j,:));
            trust_score = abs(i_score - j_score) / (norm(diff,2)+0.01);
            %disp(trust_score);
            if (i_score - j_score > 0)
                Positions(i,:) = Positions(i,:) - (uub * (1 - exp(-trust_score/10)) + llb) * diff;
            else
                Positions(i,:) = Positions(i,:) + (uub * (1 - exp(-trust_score/10)) + llb) * diff;
            end
        end
        for j = 1:dim % prevent error solution
            if Positions(i,j) > ub
                Positions(i,j) = ub;
            end
            if Positions(i,j) < lb
                Positions(i,j) = lb;
            end
        end
        if fobj(Positions(i,:)) >= fobj(old_position(i,:))
            Positions(i,:) = old_position(i,:);
            reject=reject+1;
        else
            accept=accept+1;
        end
    end
    accept_curve(l+1)=accept/(accept+reject);
    if (accept/(accept+reject) == 0)
        aaa = 11;
    end

    %disp(Alpha_pos);

    %%
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
end
figure();
plot(accept_curve);
