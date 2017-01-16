A = load('work_odom');
B = load('work_loc');
T = sqrt((A(:,1) - B(:,1)).^2 + (A(:,2) - B(:,2)).^2);
sum(T)/length(T)