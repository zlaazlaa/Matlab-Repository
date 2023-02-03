function new_wolf_pos = brownian(Positions,lb, ub,alpha_pos)
y = randi([1,size(Positions,1)]);
z = randi([1,size(Positions,1)]);
sigma=normrnd(0,1);
new_wolf_pos = alpha_pos + sigma * (Positions(y,:)-Positions(z,:));
for j = 1:size(Positions,2)
    if new_wolf_pos(j) > ub
        new_wolf_pos(j) = ub;
    end
    if new_wolf_pos(j) < lb
        new_wolf_pos(j) = lb;
    end
end