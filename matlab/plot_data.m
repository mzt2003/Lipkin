%sys1_plotE_S('../data/sys1_Energy_Entropy.txt');
%sys1_plotS_chi("../data/sys1_S_vs_chi.txt");
%sys1_plotS_n('../data/sys1_S_vs_n.txt');
%sys1_S_chi_n("../data/sys1_S_vs_chi_n.txt");
%sys2_plotEnergy('../data/sys2_Energy_Entropy.txt');
%sys2_plotEntropy('../data/sys2_Energy_Entropy.txt');
%sys12_E_S('../data/sys1_Energy_Entropy.txt',"../data/sys2_Energy_Entropy_.txt");
%sys2_gsS_v_v12("../data/sys2_S_v_v12.txt");
%sys2_gsS_v12("../data/sys2_S_v12.txt")
%sys2_gsS_v12_multi("../data/sys2_S_v12_multi.txt")
sys2_h1h2_v12("../data/sys2_h1h2_v12.txt")




function sys1_plotE_S(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
N = params(1);
chi = params(2);
n = params(3);

data = fscanf(fid, '%f %f', [2 inf])'; 
fclose(fid); 

Energy = data(:, 1);
Entropy = data(:, 2);

figure;
plot(Energy, Entropy, '-o');
xlabel('Energy');
ylabel('Entropy');
title(sprintf('Entropy vs Energy - N=%d, $\\chi$=%.2f, n=%d', N, chi, n), 'Interpreter', 'latex');
saveas(gcf, '../plot/sys1_Entropy_Energy.png');
end





function sys1_plotS_chi(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
N = params(1);
n = params(2);

data = fscanf(fid, '%f %f', [2 inf])'; 
fclose(fid); 
chi = data(:, 1);
entropy = data(:, 2);

figure;
plot(chi, entropy, '-o');
xlabel('$\chi$', 'Interpreter', 'latex');
ylabel('Ground State Entropy');
title(sprintf('Ground State Entropy vs $\\chi$ for N=%d, n=%d', N, n), 'Interpreter', 'latex');
grid on;

saveas(gcf, '../plot/sys1_entropy_vs_chi.png');
end





function sys1_plotS_n(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end
params = str2num(fgetl(fid)); 
N = params(1);
chi = params(2);

data = fscanf(fid, '%f %f', [2 inf])'; 
fclose(fid);

n = data(:, 1);
entropy = data(:, 2);

figure;
plot(n, entropy, '-o');
xlabel('n');
ylabel('Ground State Entropy');
title(sprintf('Ground State Entropy vs n for N=%d, $\\chi$=%.2f',N,chi), 'Interpreter', 'latex');
grid on;

saveas(gcf, '../plot/sys1_entropy_vs_n.png');
end






function sys1_S_chi_n(filename)
fileID = fopen(filename, 'r');

N = fscanf(fileID, '%d', 1);
fgetl(fileID);  

chis_line = fgetl(fileID);
chis = str2double(strsplit(chis_line, ' ', 'CollapseDelimiters', true));  
chis(isnan(chis)) = []; 

ns_line = fgetl(fileID);
ns = str2double(strsplit(ns_line, ' ', 'CollapseDelimiters', true));
ns(isnan(ns)) = [];  

gs_S = [];
while ~feof(fileID)
    line = fgetl(fileID);
    if ~ischar(line)  
        break;
    end
    currentRow = str2double(strsplit(line, ' ', 'CollapseDelimiters', true));
    currentRow(isnan(currentRow)) = [];
    gs_S = [gs_S; currentRow];  
end

fclose(fileID);  
[X, Y] = meshgrid(chis, ns);
h= surf(X, Y, gs_S);
set(h, 'EdgeColor', 'none'); 

xlabel('$\chi$', 'Interpreter', 'latex');
ylabel('n', 'Interpreter', 'latex');
zlabel('S', 'Interpreter', 'latex');
title(sprintf('ground state entanglement entropy for N = %d', N), 'Interpreter', 'latex');
saveas(gcf, '../plot/sys1_entropy_vs_chi_n.png');
end






function sys2_plotEnergy(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
n1 = params(1);
n2 = params(2);
ep1 = params(3);
ep2 = params(4);
V1 = params(5);
V2 = params(6);
V12 = params(7);

data = fscanf(fid, '%f %f %f %f', [4, inf])'; 
fclose(fid); 
data
Energy = data(:, 1);
h1 = data(:, 3);
h2= data(:,4);
S= data(:, 2);

figure;
plot(Energy, h1, '-o');
hold on
plot(Energy, h2);
hold off
xlabel('Total Energy', 'Interpreter', 'latex');
ylabel('$\langle H_1 \rangle$ and $\langle H_2 \rangle$', 'Interpreter', 'latex');
title(' Total Energy vs $\langle H_1 \rangle$ and $\langle H_2 \rangle$', 'Interpreter', 'latex');
legend({'$\langle H_1 \rangle$', '$\langle H_2 \rangle$'}, 'Interpreter', 'latex');
str = sprintf('n_1=%d, n_2=%d\n\\epsilon_1=%.2f, \\epsilon_2=%.2f\nV_1=%.2f, V_2=%.2f, V_{12}=%.2f', ...
              n1, n2, ep1, ep2, V1, V2, V12);
text(max(Energy)*0.1+min(Energy)*0.9, max([h1; h2])*0.9, str, 'FontSize', 8); 
saveas(gcf, strcat('../plot/sys2_energy_h1h2',num2str(n1),num2str(n2),num2str(ep1),num2str(ep2),num2str(V1),num2str(V2),num2str(V12),'.png'));
end





function sys2_plotEntropy(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
n1 = params(1);
n2 = params(2);
ep1 = params(3);
ep2 = params(4);
V1 = params(5);
V2 = params(6);
V12 = params(7);

data = fscanf(fid, '%f %f %f %f', [4, inf])'; 
fclose(fid); 
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

saveas(gcf, strcat('../plot/sys2_S_Eh1h2_',num2str(n1),num2str(n2),num2str(ep1),num2str(ep2),num2str(V1),num2str(V2),num2str(V12),'.png'));
end




function sys12_E_S(filename,file2)

fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
N = params(1);
chi = params(2);
n = params(3);

data = fscanf(fid, '%f %f', [2 inf])'; 
fclose(fid); 

Energy = data(:, 1);
Entropy = data(:, 2);

figure;
plot(Energy, Entropy, '-o');
xlabel('Energy');
ylabel('Entropy');
title(sprintf('Entropy vs Energy - N=%d, $\\chi$=%.2f, n=%d', N, chi, n), 'Interpreter', 'latex');
hold on
fid = fopen(file2, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 

data2 = fscanf(fid, '%f %f', [4 inf])'; 
fclose(fid); 
data2
scatter(data2(:,1), data2(:,2));
hold off

saveas(gcf, strcat('../plot/sys12_energy_entropy_N_',num2str(N),'chi',num2str(chi),'n',num2str(n),'.png'));
end





function sys2_gsS_v_v12(filename)

fileID = fopen(filename, 'r');
N = fscanf(fileID, '%d', 1);
fgetl(fileID);  

chis_line = fgetl(fileID);
chis = str2double(strsplit(chis_line, ' ', 'CollapseDelimiters', true));  
chis(isnan(chis)) = [];  

ns_line = fgetl(fileID);
ns = str2double(strsplit(ns_line, ' ', 'CollapseDelimiters', true));
ns(isnan(ns)) = [];  

gs_S = [];
while ~feof(fileID)
    line = fgetl(fileID);
    if ~ischar(line)  
        break;
    end
    currentRow = str2double(strsplit(line, ' ', 'CollapseDelimiters', true));
    currentRow(isnan(currentRow)) = [];
    gs_S = [gs_S; currentRow];  
end

fclose(fileID);  
[X, Y] = meshgrid(chis, ns);
h= surf(X, Y, gs_S);
set(h, 'EdgeColor', 'none'); 

xlabel('$V\times N$', 'Interpreter', 'latex');
ylabel('$V_{12} \times N$', 'Interpreter', 'latex');
zlabel('$S$', 'Interpreter', 'latex');
title(sprintf('ground state entanglement entropy for N = %d', N));
saveas(gcf, strcat('../plot/sys2_gsS_v_v12_N',num2str(N),'.png'));

end




function sys2_gsS_v12(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
n1 = params(1);
n2 = params(2);
ep1 = params(3);
ep2 = params(4);
V1 = params(5);
V2 = params(6);

data = fscanf(fid, '%f %f', [2 inf])'; 
fclose(fid); 
chi = data(:, 1);
entropy = data(:, 2);

figure;
plot(chi, entropy, '-o');
xlabel('$V_{12}$', 'Interpreter', 'latex');
ylabel('Ground State Entropy');
grid on;
str = sprintf('n_1=%d, n_2=%d\n\\epsilon_1=%.2f, \\epsilon_2=%.2f\nV_1=%.4f, V_2=%.4f', ...
              n1, n2, ep1, ep2, V1, V2);
text(max(chi)*0.7+min(chi)*0.3, max(entropy)*0.5, str, 'FontSize', 8); 
saveas(gcf, strcat('../plot/sys2_entropy_vs_v12',num2str(n1),num2str(n2),num2str(ep1),num2str(ep2),num2str(V1),num2str(V2),'.png'));
end





function sys2_gsS_v12_multi(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end


params = str2num(fgetl(fid));
n1 = params(1);
n2 = params(2);
ep1 = params(3);
ep2 = params(4);
V_params = params(5:end);
num_v = length(V_params) / 2;
V1s = V_params(1:num_v);
V2s = V_params(num_v+1:end);

data = fscanf(fid, '%f %f', [1+num_v inf])';
fclose(fid); 
chi = data(:, 1);
gs_entropies = data(:, 2:end);

figure;
hold on;
colors = lines(num_v); 
legends = cell(1, num_v);

for i = 1:num_v
    plot(chi*(n1+n2) , gs_entropies(:, i), '-o', 'Color', colors(i, :));
    legends{i} = sprintf('V_1\\times N = %.4f, V_2\\times N = %.4f', V1s(i)*(n1+n2), V2s(i)*(n1+n2));
end

xlabel('$V_{12}\times N $', 'Interpreter', 'latex');
ylabel('Ground State Entropy', 'Interpreter', 'latex');
title(sprintf('$n_1=%d$, $n_2=%d$, $\\epsilon_1=%.2f$, $\\epsilon_2=%.2f$', n1, n2, ep1, ep2), 'Interpreter', 'latex');
legend(legends, 'Location', 'best');
grid on;
hold off
%str = sprintf('n_1=%d, n_2=%d\n\\epsilon_1=%.2f, \\epsilon_2=%.2f', ...
%              n1, n2, ep1, ep2);
%text(max(chi*(n1+n2))*0.1+min(chi*(n1+n2))*0.9, max(max(gs_entropies))*0.9, str, 'FontSize', 8); 
saveas(gcf, strcat('../plot/sys2_entropy_vs_v12_multi', num2str(n1), num2str(n2), num2str(ep1), num2str(ep2), '.png'));
end










function sys2_h1h2_v12(filename)
fid = fopen(filename, 'rt');
if fid == -1
    error('File could not be opened.');
end

params = str2num(fgetl(fid)); 
n1 = params(1);
n2 = params(2);
ep1 = params(3);
ep2 = params(4);
V1 = params(5);
V2 = params(6);

data = fscanf(fid, '%f %f', [3 inf])'; 
fclose(fid); 
chi = data(:, 1);
h1 = data(:, 2);
h2 = data(:, 3);


figure;
plot(chi.*(n1+n2), h1, '-o');
hold on
plot(chi.*(n1+n2), h2)
legend({'$\langle H_1 \rangle$','$\langle H_2 \rangle$'},'Interpreter', 'latex')
xlabel('$V_{12}\times N $', 'Interpreter', 'latex');
ylabel('Energy')
title('\langle H_1 \rangle and \langle H_2 \rangle of Ground State');
str=sprintf('n_1=%d, n_2=%d,\n \\epsilon_1=%.2f, \\epsilon_2=%.2f,\n V_1\\times N =%.2f, V_2\\times N=%.2f', n1, n2, ep1, ep2, V1*(n1+n2), V2*(n1+n2));
text(max(chi*(n1+n2))*0.1+min(chi*(n1+n2))*0.9, max([max(h1),max(h2)])*0.9+min([min(h1),min(h2)])*0.1, str, 'FontSize', 8); 
grid on;
saveas(gcf, strcat('../plot/sys2_h1h2_vs_v12',num2str(n1),num2str(n2),num2str(ep1),num2str(ep2),num2str(V1),num2str(V2),'.png'));
end
