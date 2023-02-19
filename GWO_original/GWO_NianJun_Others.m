function [Alpha_score,Alpha_pos,Convergence_curve]=GWO_NianJun_Others(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

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

best_of_every=Positions;
best_score_of_every=Inf(1,SearchAgents_no);
best_of_whole=Positions(1,:);
best_score_of_whole=inf;
figure();
plot(5,2);
% Main loop
while l<Max_iter
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        if (fitness < best_score_of_whole)
            best_score_of_whole = fitness;
            best_of_whole = Positions(i,:);
        end
        if (fitness < best_score_of_every(i))
            best_score_of_every(i) = fitness;
            best_of_every(i,:) = Positions(i,:);
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
            Delta_id = i;
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
            
            w1 = norm(X1,2)/(norm(X1,2) + norm(X2,2) + norm(X3,2));
            w2 = norm(X2,2)/(norm(X1,2) + norm(X2,2) + norm(X3,2));
            w3 = norm(X3,2)/(norm(X1,2) + norm(X2,2) + norm(X3,2));
            Positions(i,j)=(w1*X1+w2*X2+w3*X3)/3;% Equation (3.7)
            
        end
    end

    %% Myxomycetes information exchange
    old_position = Positions;
    accept=0;
    reject=0;
    llb = 0.4*(Max_iter - l)/Max_iter + 0.2;
    llb = llb * rand();
%     llb = 0.2*(Max_iter - l)/Max_iter;
    uub = 1-llb;
    for i = 1:SearchAgents_no % position to be changed
        for w = 1:round(SearchAgents_no/10) % neighbor selected
            j = randi([1,SearchAgents_no]);
            diff = old_position(i,:) - old_position(j,:);
            i_score = fobj(old_position(i,:));
            j_score = fobj(old_position(j,:));
            trust_score = abs(i_score - j_score) / (norm(diff,2)+0.01);
            %disp(trust_score);
            if (i_score - j_score > 0)
                Positions(i,:) = Positions(i,:) - (uub * (1 - exp(-trust_score/1)) + llb) * diff;
            else
                Positions(i,:) = Positions(i,:) + (uub * (1 - exp(-trust_score/1)) + llb) * diff;
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
    if (accept_curve(l+1) == 1)
        asdasd = 1231;
    end
    if (accept/(accept+reject) < 0.5)
%         disp('rebuild');
%         for w = 1:round(SearchAgents_no/10)
%             old_version = Positions(w,:);
%             rand_num = randi([1,SearchAgents_no]);
%             DIFF = Positions(w,:) - Positions(rand_num,:);
%             DIFF = DIFF / norm(DIFF,2);
% %             Positions(w,:) = Positions(rand_num,:)+0.1*levy_flight()*(Worse_pos-Positions(rand_num,:));
%             Positions(w,:) = Positions(w,:) + 0.5 * levy_flight() * DIFF;
%             for j = 1:size(Positions,2)
%                 if (Positions(w,j) > ub)
%                     Positions(w,j) = ub;
%                 end
%                 if (Positions(w,j) < lb)
%                     Positions(w,j) = lb;
%                 end
%             end
%             if(fobj(Positions(w,:)) > fobj(old_version))    
%                 Positions(w,:) = old_version;
%             end
%         end
% 
%         disp('rebuild2-PSO');
%         for w = 1:round(SearchAgents_no/5)
%             old_version = Positions(w,:);
%             Positions(w,:) = Positions(w,:) + 2 * rand() * (best_of_whole - Positions(w,:)) + 2 * rand() * (best_of_every(w,:) - Positions(w,:));
%             for j = 1:size(Positions,2)
%                 if (Positions(w,j) > ub)
%                     Positions(w,j) = ub;
%                 end
%                 if (Positions(w,j) < lb)
%                     Positions(w,j) = lb;
%                 end
%             end
%             if(fobj(Positions(w,:)) > fobj(old_version))    
%                 Positions(w,:) = old_version;
%             end
%         end

%         disp('rebuild3-levy-on-single-wolf')
% %         rand_wolf = randi([1,SearchAgents_no]);
%         rand_wolf = Delta_id;
%         rand_wolf_pos = Positions(rand_wolf,:);
%         Delta_score = fobj(rand_wolf_pos);
%         for t=1:Max_iter %%%%%%% 这里也许可以改
%             b=1;
%             a=0;
%             r = (b-a).*rand(dim,1) + a;
%             r = r/norm(r,2);
%             change_k = 0.5 * levy_flight();
%             if abs(change_k) > 0.6 * ub
%                 change_k = 0.5 * levy_flight();
%                 disp('up');
%             end
%             rand_wolf_pos = rand_wolf_pos + change_k * r.';
%             for j = 1:dim
%                 if rand_wolf_pos(j) > ub
%                     rand_wolf_pos(j) = ub;
%                 end
%                 if rand_wolf_pos(j) < lb
%                     rand_wolf_pos(j) = lb;
%                 end
%             end
%             new_score = fobj(rand_wolf_pos);
%             if new_score < Delta_score
%                 Delta_score = new_score;
%                 Positions(rand_wolf,:) = rand_wolf_pos;
%             end
%         end
%         plot(l+1,accept_curve(l+1),'.','Color','g','MarkerSize',30);
%         hold on;
% %         text(l+1,accept_curve(l+1),'q');
    end

    %disp(Alpha_pos);

    %%
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
end
plot(accept_curve);
% text(10,3,'qeeq');
