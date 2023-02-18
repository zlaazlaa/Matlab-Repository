function [Alpha_score,Alpha_pos,Convergence_curve]=my_new_GWO_master_slave(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

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


Alpha_pos1=zeros(1,dim);
Alpha_score1=inf; %change this to -inf for maximization problems

Beta_pos1=zeros(1,dim);
Beta_score1=inf; %change this to -inf for maximization problems

Delta_pos1=zeros(1,dim);
Delta_score1=inf; %change this to -inf for maximization problems


Alpha_pos2=zeros(1,dim);
Alpha_score2=inf; %change this to -inf for maximization problems

Beta_pos2=zeros(1,dim);
Beta_score2=inf; %change this to -inf for maximization problems

Delta_pos2=zeros(1,dim);
Delta_score2=inf; %change this to -inf for maximization problems
% Main loop
best_best=inf;
num_1=0;
num_2=0;
group_num=10;
Alpha_score_list = zeros(1,group_num);
while l<Max_iter
    group_size = round(SearchAgents_no/group_num);
    for i = 1:group_num
        [Alpha_score_list(i),Positions] = Sub_GWO(group_size*(i-1) + 1,group_size*(i-1)+group_size,Positions,fobj,ub,lb,dim,Max_iter,l);
    end
    for i=1:group_num
        if best_best > Alpha_score_list(i)
            best_best = Alpha_score_list(i);
        end
    end
    l=l+1;
    [~,index] = sort(Alpha_score_list);
    % update bad group
    for i = round(group_num/2):group_num
        orgin_i = index(i);
        for j = group_size*(orgin_i-1) + 1:group_size*orgin_i
            rand_num = randi([group_size*(orgin_i-1) + 1, group_size*orgin_i]);
            Positions(j,:) = update_wolf(l,Max_iter,Positions,dim,fobj,lb,ub,fobj(Positions(rand_num,:)),Positions(rand_num,:));
        end
    end
    
%     if (Alpha_score1 > Alpha_score2)
%         num_1 = num_1 + 1;
%         num_2 = 0;z
%         Alpha_score1 = Alpha_score2;
%     else
%         num_2 = num_2 + 1;
%         num_1 = 0;
%     end
%     if (best_best > Alpha_score1)
%         best_best = Alpha_score1;
%     end

%     disp(num_1);
    % update wolf if it sleep for a long time
%     if(num_1 > 2)
%         disp('sleep1');
%         for i = 1:round(SearchAgents_no/2)
%             rand_num = randi([1,round(SearchAgents_no/2)]);
%             Positions(i,:) = update_wolf(l,Max_iter,Positions,dim,fobj,lb,ub,fobj(Positions(rand_num,:)),Positions(rand_num,:));
%         end
%     end
%     if(num_2 > 2)
%         disp('sleep2');
%         for i = round(SearchAgents_no/2)+1:SearchAgents_no
%             rand_num = randi([round(SearchAgents_no/2)+1,SearchAgents_no]);
%             Positions(i,:) = update_wolf(l,Max_iter,Positions,dim,fobj,lb,ub,fobj(Positions(rand_num,:)),Positions(rand_num,:));
%         end
%     end


    
    Convergence_curve(l)=best_best;
end