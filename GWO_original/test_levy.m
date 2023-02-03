curvy = zeros(1000);
for i =1:1000
    curvy(i) = levy_flight();
end
plot(curvy)