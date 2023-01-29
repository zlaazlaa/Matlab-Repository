function [new_wolf, new_score] = update_wolf(l, Max_iter, Positions, dim, fobj, lb, ub, Alpha_score, Alpha_pos)
    % update alpha wolf according to Differential Evolution Algorithm(DE)
    mom = Positions(randi(size(Positions, 1)),:);
    dad = Positions(randi(size(Positions, 1)),:);
    mid_wolf = mom - dad;

    % calculate F
    F_l = 0.1;
    F_u = 0.9;
    F = F_l + (Max_iter - l) * (F_u - F_l) / Max_iter;
    mid_wolf = Alpha_pos + F * (mid_wolf); % 中间向量
    
    % calculate cr
    cr_l = 0.1;
    cr_u = 0.9;
    cr = cr_l + (Max_iter - l) * (cr_u - cr_l) / Max_iter;
    
    % update mid_wolf
    for i = 1: dim
        if rand() > cr
            mid_wolf(i) = Alpha_pos(i);
        end
    end

    % whether accept the new alpha wolf according to SA
    new_alpha_score = fobj(mid_wolf);
    Delta = new_alpha_score - Alpha_score;
%     if (Delta < 0) 
%         Alpha_score = new_alpha_score;
%         Alpha_pos = mid_wolf;
%     end
    if Delta < 0  % uphill move (good move)
        disp('better')
        Alpha_score = new_alpha_score;
        Alpha_pos = mid_wolf;
    else % downhill move (bad move)
        P=exp(-Delta/ (l + 1));
        disp('same')
        if rand()<=P
            disp('worse')
            Alpha_score = new_alpha_score;
            Alpha_pos = mid_wolf;
        end
    end

    for i=1:dim
        if Alpha_pos(i) < lb || Alpha_pos(i) > ub
            Alpha_pos(i) = lb + rand() * (ub - lb);
        end
    end
    % end the update of alpha wolf
    new_wolf = Alpha_pos;
    new_score = Alpha_score;
end