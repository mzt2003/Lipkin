% 打开文件
fid = fopen('../data/data.txt', 'rt');
if fid == -1
    error('File could not be opened.');
end

% 读取参数
params = str2num(fgetl(fid)); % 读取第一行并转换成数字
N = params(1);
chi = params(2);
n = params(3);

% 读取数据
data = fscanf(fid, '%f %f', [2 inf])'; % 读取其余数据
fclose(fid); % 关闭文件

Energy = data(:, 1);
Entropy = data(:, 2);

% 画图
figure;
plot(Energy, Entropy, '-o');
xlabel('Energy');
ylabel('Entropy');
title(sprintf('Energy vs Entropy - N=%d, chi=%.2f, n=%d', N, chi, n));

% 保存图形到文件
saveas(gcf, strcat('../data/energy_entropy_N‘',num2str(N),'chi',num2str(chi),'n',num2str(n),'.png'));

% 读取数据
data = load('../data/entropy_vs_chi.txt');

% 解析数据
chi = data(:, 1);
entropy = data(:, 2);

% 绘图
figure;
plot(chi, entropy, '-o');
xlabel('Chi');
ylabel('Ground State Entropy');
title('Ground State Entropy vs Chi for N=100, n=20');
grid on;

% 保存图像
saveas(gcf, '../data/entropy_vs_chi.png');
