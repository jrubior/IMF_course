%% store results and plot IRFs


Ltildeq5vararrwIDEN150=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 50th
Ltildeq5vararrwIDEN116=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 16th
Ltildeq5vararrwIDEN184=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 84th


for ii=1:size(Ltilde,1)
    for jj=1:size(Ltilde,2)
        for kk=1:size(Ltilde,3)
            Ltildeq5vararrwIDEN150(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.5);
            Ltildeq5vararrwIDEN116(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.16);
            Ltildeq5vararrwIDEN184(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.84);
        end
    end
end

scale=1;
fontsizetitle  = 8;%9;
fontsizeaxes   = 8;%9;
fontsizelegend = 8;%9;
axiswidth      = 1;


hFig = figure('name','Panel (b) Figure 2','NumberTitle','off');
set(hFig, 'Position', [0 20 500 250])
subplot(2,3,1)
plot(0:1:horizon,squeeze(Ltildeq5vararrwIDEN150(1:horizon+1,1,1))*scale,'LineWidth',2,'Color',[0    1.0000    0.4961])
hline(0,'-r')
xlabel('Quarters')
ylabel('Percent')
hold on
a=(squeeze(Ltildeq5vararrwIDEN116(1:horizon+1,1,1))*scale)';
b=(squeeze(Ltildeq5vararrwIDEN184(1:horizon+1,1,1))*scale)';
x = 0:1:horizon;
[~,~]=jbfill(x,a,b,[0.6602    0.6602    0.6602],[0.6602    0.6602    0.6602],0,0.5);
set(gca,'XTick',[0;10;20;30;40])
set(gca,'XTickLabel',['0 ';'10';'20';'30';'40'])
set(gca,'YTick',[-0.35 0 0.35 0.7 1.05])
axis([0 40 -0.35 1.05])
set(gca,'LineWidth',axiswidth)
set(gca,'FontSize',fontsizeaxes)
grid on
box off
H=gca;
set(H,'Gridalpha',0.05);
title('Adjusted TFP','Interpreter','tex','FontSize',fontsizetitle)

subplot(2,3,2)
plot(0:1:horizon,squeeze(Ltildeq5vararrwIDEN150(1:horizon+1,2,1))*scale,'LineWidth',2,'Color',[0    1.0000    0.4961])
hline(0,'-r')
xlabel('Quarters')
ylabel('Percent')
hold on
a=(squeeze(Ltildeq5vararrwIDEN116(1:horizon+1,2,1))*scale)';
b=(squeeze(Ltildeq5vararrwIDEN184(1:horizon+1,2,1))*scale)';
x = 0:1:horizon;
[~,~]=jbfill(x,a,b,[0.6602    0.6602    0.6602],[0.6602    0.6602    0.6602],0,0.5);
set(gca,'XTick',[0;10;20;30;40])
set(gca,'XTickLabel',['0 ';'10';'20';'30';'40'])
set(gca,'YTick',[-3.5 0 3.5 7 10.5])
axis([0 40 -3.5 10.5])
set(gca,'LineWidth',axiswidth)
set(gca,'FontSize',fontsizeaxes)
grid on
box off
H=gca;
set(H,'Gridalpha',0.05);
title('Stock Prices','Interpreter','tex','FontSize',fontsizetitle)

subplot(2,3,3)
plot(0:1:horizon,squeeze(Ltildeq5vararrwIDEN150(1:horizon+1,3,1))*scale,'LineWidth',2,'Color',[0    1.0000    0.4961])
hline(0,'-r')
xlabel('Quarters')
ylabel('Percent')
hold on
a=(squeeze(Ltildeq5vararrwIDEN116(1:horizon+1,3,1))*scale)';
b=(squeeze(Ltildeq5vararrwIDEN184(1:horizon+1,3,1))*scale)';
x = 0:1:horizon;
[~,~]=jbfill(x,a,b,[0.6602    0.6602    0.6602],[0.6602    0.6602    0.6602],0,0.5);
set(gca,'XTick',[0;10;20;30;40])
set(gca,'XTickLabel',['0 ';'10';'20';'30';'40'])
set(gca,'YTick',[-0.4 0 0.40 0.8 1.2])
axis([0 40 -0.4 1.2])
set(gca,'LineWidth',axiswidth)
set(gca,'FontSize',fontsizeaxes)
grid on
box off
H=gca;
set(H,'Gridalpha',0.05);
title('Consumption','Interpreter','tex','FontSize',fontsizetitle)

subplot(2,3,4)
plot(0:1:horizon,squeeze(Ltildeq5vararrwIDEN150(1:horizon+1,4,1))*scale,'LineWidth',2,'Color',[0    1.0000    0.4961])
hline(0,'-r')
xlabel('Quarters')
ylabel('Percentage points')
hold on
a=(squeeze(Ltildeq5vararrwIDEN116(1:horizon+1,4,1))*scale)';
b=(squeeze(Ltildeq5vararrwIDEN184(1:horizon+1,4,1))*scale)';
x = 0:1:horizon;
[~,~]=jbfill(x,a,b,[0.6602    0.6602    0.6602],[0.6602    0.6602    0.6602],0,0.5);
set(gca,'XTick',[0;10;20;30;40])
set(gca,'XTickLabel',['0 ';'10';'20';'30';'40'])
set(gca,'YTick',[-1.1 -0.55 0 0.55 1.1])
axis([0 40 -1.1 1.1])
set(gca,'LineWidth',axiswidth)
set(gca,'FontSize',fontsizeaxes)
grid on
box off
H=gca;
set(H,'Gridalpha',0.05);
title('Real Interest Rate','Interpreter','tex','FontSize',fontsizetitle)

subplot(2,3,5)
plot(0:1:horizon,squeeze(Ltildeq5vararrwIDEN150(1:horizon+1,5,1))*scale,'LineWidth',2,'Color',[0    1.0000    0.4961])
hline(0,'-r')
xlabel('Quarters')
ylabel('Percent')
hold on
a=(squeeze(Ltildeq5vararrwIDEN116(1:horizon+1,5,1))*scale)';
b=(squeeze(Ltildeq5vararrwIDEN184(1:horizon+1,5,1))*scale)';
x = 0:1:horizon;
[~,~]=jbfill(x,a,b,[0.6602    0.6602    0.6602],[0.6602    0.6602    0.6602],1,0.5);
set(gca,'XTick',[0;10;20;30;40])
set(gca,'XTickLabel',['0 ';'10';'20';'30';'40'])
set(gca,'YTick',[-0.5 0 0.5 1 1.5])
axis([0 40 -0.5 1.5])
set(gca,'LineWidth',axiswidth)
set(gca,'FontSize',fontsizeaxes)
grid on
box off
H=gca;
set(H,'Gridalpha',0.05);
title('Hours Worked','Interpreter','tex','FontSize',fontsizetitle)

set(gcf, 'PaperPositionMode', 'auto');
print('results/pngfiles/figure_1_panel_b.png','-dpng');
print('results/epsfiles/figure_1_panel_b.eps','-depsc');

savefile='results/matfiles/results.mat';
save(savefile,'Ltilde','FEVD','horizon','imp_w','uw','ne');





