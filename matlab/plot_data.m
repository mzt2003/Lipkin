% 打开文件
fid = fopen('../data/Energy_Entropy.txt', 'rt');
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
%hold on
%data2 = load("../data/Energy_Entropy22.txt");
%scatter(data2(:,1), data2(:,2));
%hold off
% 保存图形到文件
saveas(gcf, strcat('../data/energy_entropy_N',num2str(N),'chi',num2str(chi),'n',num2str(n),'.png'));

%%%%%%%%%%%%%%%%%%%%%%%%
% 打开文件
fid = fopen('../data/S_vs_chi.txt', 'rt');
if fid == -1
    error('File could not be opened.');
end

% 读取参数
params = str2num(fgetl(fid)); % 读取第一行并转换成数字
N = params(1);
n = params(2);

% 读取数据
data = fscanf(fid, '%f %f', [2 inf])'; % 读取其余数据
fclose(fid); % 关闭文件

% 解析数据
chi = data(:, 1);
entropy = data(:, 2);

% 绘图
figure;
plot(chi, entropy, '-o');
xlabel('Chi');
ylabel('Ground State Entropy');
title(sprintf('Ground State Entropy vs Chi for N=%d, n=%d',N,n));
grid on;

% 保存图像
saveas(gcf, '../data/entropy_vs_chi.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 打开文件
fid = fopen('../data/S_vs_n.txt', 'rt');
if fid == -1
    error('File could not be opened.');
end

% 读取参数
params = str2num(fgetl(fid)); % 读取第一行并转换成数字
N = params(1);
chi = params(2);

% 读取数据
data = fscanf(fid, '%f %f', [2 inf])'; % 读取其余数据
fclose(fid); % 关闭文件

% 解析数据
n = data(:, 1);
entropy = data(:, 2);

% 绘图
figure;
plot(n, entropy, '-o');
xlabel('n');
ylabel('Ground State Entropy');
title(sprintf('Ground State Entropy vs Chi for N=%d, chi=%.2f',N,chi));
grid on;

% 保存图像
saveas(gcf, '../data/entropy_vs_n.png');


%%%%%%%%%%%%%%%%%
% 读取文件
filename = "../data/S_vs_chi_n.txt";
fileID = fopen(filename, 'r');

N = fscanf(fileID, '%d', 1);
fgetl(fileID);  % 读取剩余部分并移动到下一行

% 读取第二行（chis数组）
chis_line = fgetl(fileID);
chis = str2double(strsplit(chis_line, ' ', 'CollapseDelimiters', true));  
  % 将字符串分割并转换为数字数组
chis(isnan(chis)) = [];  % 移除任何NaN元素

% 读取第三行（ns数组）
ns_line = fgetl(fileID);
ns = str2double(strsplit(ns_line, ' ', 'CollapseDelimiters', true));
ns(isnan(ns)) = [];  % 移除任何NaN元素

% 读取后面的行（gs_S矩阵）
gs_S = [];
while ~feof(fileID)
    line = fgetl(fileID);
    if ~ischar(line)  % 如果读取到文件末尾，则退出循环
        break;
    end
    currentRow = str2double(strsplit(line, ' ', 'CollapseDelimiters', true));
    currentRow(isnan(currentRow)) = [];
    gs_S = [gs_S; currentRow];  % 转换行为数字数组并追加到gs_S
end

fclose(fileID);  % 关闭文件
% 绘制三维图
[X, Y] = meshgrid(chis, ns);
surf(X, Y, gs_S);


xlabel('chi');
ylabel('n');
zlabel('S');
title(sprintf('ground state entanglement entropy for N = %d', N));

saveas(gcf, '../data/entropy_vs_chi_n.png');





%data2 = load("../data/Energy_Entropy22.txt");
%scatter(data2(:,1), data2(:,2));
%xlabel('Eigenvalues');
%ylabel('Entropy');
%title('Entropy vs Eigenvalues');
%saveas(gcf, '../data/Entropy_Eigenvalues22.png');

