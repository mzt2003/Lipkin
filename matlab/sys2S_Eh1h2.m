function plotEntropy()
% 打开文件
fid = fopen('../data/Energy_Entropy33.txt', 'rt');
if fid == -1
    error('File could not be opened.');
end

% 读取参数
params = str2num(fgetl(fid)); % 读取第一行并转换成数字
n1 = params(1);
n2 = params(2);
ep1 = params(3);
ep2 = params(4);
V1 = params(5);
V2 = params(6);
V12 = params(7);

% 读取数据
data = fscanf(fid, '%f %f %f %f', [4, inf])'; % 读取其余数据
fclose(fid); % 关闭文件
data
Energy = data(:, 1);
h1 = data(:, 3);
h2= data(:,4);
S= data(:, 2);


figure;
plot(S, Energy);
hold on
plot(S,h1);
plot(S, h2);
hold off
xlabel('Entropy')
ylabel('$E$, $\langle H_1 \rangle$ and $\langle H_2 \rangle$', 'Interpreter', 'latex');
legend({'Total Energy', '$\langle H_1 \rangle$', '$\langle H_2 \rangle$'}, 'Interpreter', 'latex');

saveas(gcf, strcat('../data/S_Eh1h2_',num2str(n1),num2str(n2),num2str(ep1),num2str(ep2),num2str(V1),num2str(V2),num2str(V12),'.png'));




end
