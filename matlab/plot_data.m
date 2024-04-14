data = load('../data/data.txt');
x = data(:,1);
y = data(:,2);

plot(x, y);
title('Plot of Y vs X');
xlabel('X-axis');
ylabel('Y-axis');
