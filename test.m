%% 1. 假设已有数据
dof    = [10, 20, 40, 80, 160, 320, 640, 1280];
h1err  = [0.5, 0.28, 0.15, 0.08, 0.042, 0.021, 0.01, 0.005];
h1rerr = [0.4, 0.15, 0.07, 0.035, 0.018, 0.009, 0.004, 0.002];

%% 2. 对数坐标下绘制两条曲线

figure('Color','white');
loglog(dof, h1err,  'ro-','LineWidth',1.5,'MarkerSize',6);  
hold on;
loglog(dof, h1rerr, 'b*-','LineWidth',1.5,'MarkerSize',6);  
grid on;
xlabel('Number of Dof','FontSize',12);
ylabel('Error','FontSize',12);
legend('De','De^{r_2}','Location','SouthWest');
title('Error Convergence','FontSize',13);

% 适度限定坐标范围（看自己实际数据）
xlim([1e1,1e4]);
ylim([1e-3,1]);

%% 3. 在图中添加一个“slope = 1”的三角形

% (1) 选定一个参考点(x1,y1)（在数据坐标下）。
%     例如令x1=1000, y1=1e-2（可根据自己曲线大致中部区域选取）。
x1 = 1e3;
y1 = 1e-2;

% (2) 决定水平延伸倍数 alpha（通常可以取2、3等）。此处示例取 alpha=2。
alpha = 2;

% (3) 若想在loglog图上指示“slope = 1”，
%     则从(x1,y1)到(x2,y2)应满足：log10(y2)-log10(y1) = -1*( log10(x2)-log10(x1) ).
%     由于 x2 = alpha*x1, 所以 y2 = y1 / alpha^(slope).
slopeVal = 1;                       % 想要标识的斜率
x2 = x1 * alpha;
y2 = y1 / (alpha^slopeVal);

% (4) 画一个“小三角形”。做法之一是依次画三条线：
%     - 水平边： (x1,y1) 到 (x2,y1)
%     - 竖直边： (x2,y1) 到 (x2,y2)
%     - 斜边  ： (x1,y1) 到 (x2,y2)  (可选)
plot([x1, x2],[y1, y1],'k-','LineWidth',1.2,'HandleVisibility','off');    % 底边
plot([x2, x2],[y1, y2],'k-','LineWidth',1.2,'HandleVisibility','off');    % 竖边

% 可以只保留这两条边，让读者直观看到长宽之比；也可加斜边：
plot([x1, x2],[y1, y2],'k-','LineWidth',1.2,'HandleVisibility','off');    % 斜边(可选)

% 在三角形旁边写文字，例如 slope = 1
text(1.2*x1, 0.9*y1, '1', 'FontSize',10, 'Color','k');  

%% 4. 如果还想标 slope = 0.5 的三角形，再添加一组类似操作
% 只要更改 slopeVal = 0.5 即可：
slopeVal2 = 0.5;
x3 = 2000;              % 另一个参考点
y3 = 2e-2;              % 自己再选一下
alpha2 = 2;
x4 = x3 * alpha2;
y4 = y3 / (alpha2^slopeVal2);

plot([x3, x4],[y3, y3],'k-','LineWidth',1.2,'HandleVisibility','off'); 
plot([x4, x4],[y3, y4],'k-','LineWidth',1.2,'HandleVisibility','off'); 
plot([x3, x4],[y3, y4],'k-','LineWidth',1.2,'HandleVisibility','off'); 
text(1.2*x3,0.9*y3,'0.5','FontSize',10,'Color','k');
