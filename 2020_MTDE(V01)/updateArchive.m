function [LifeTime_archive, lifeTimeCounter] = updateArchive(LifeTime_archive, pop, fitVal, lifeTimeCounter)

if LifeTime_archive.N == 0, return; end

if size(pop, 1) ~= size(fitVal,1), error('check it'); end

%% Remove duplicate elements
popAll = [LifeTime_archive.X; pop];
funvalues = [LifeTime_archive.Fit; fitVal];
[~, IX]= unique(popAll, 'rows');
IX = sort(IX, 'ascend');
if length(IX) < size(popAll, 1) 
    popAll = popAll(IX, :);
    funvalues = funvalues(IX, :);
end

%% Add all new individuals or Remove extra individuals with higher lifeTime value
if size(popAll, 1) <= LifeTime_archive.N      
    LifeTime_archive.X = popAll;
    LifeTime_archive.Fit = funvalues;
    lifeTimeCounter(1:size(popAll, 1),:) = lifeTimeCounter(1:size(popAll, 1),:) + 1;
else                   
    size1 = size(popAll, 1);
    remove = size1 - LifeTime_archive.N;
    LifeTime_archive.X = popAll(remove+1:end , :);
    LifeTime_archive.Fit = funvalues(remove+1:end, :);
    lifeTimeCounter(1:LifeTime_archive.N-remove,:) = lifeTimeCounter(remove+1:end,:) + 1;
    lifeTimeCounter(LifeTime_archive.N-remove+1:end,:) = 1;
end
