close all;
clear all;
load('snap25.mat');
e_pod25 = sqrt(sum((w_test-wr_test_interpPROM).^2,1))./sqrt(sum((w_test).^2,1));
e_dmd25 = sqrt(sum((w_test-wr_test_interpPROM_dmd).^2,1))./sqrt(sum((w_test).^2,1));
load('snap50.mat');
e_pod50 = sqrt(sum((w_test-wr_test_interpPROM).^2,1))./sqrt(sum((w_test).^2,1));
e_dmd50 = sqrt(sum((w_test-wr_test_interpPROM_dmd).^2,1))./sqrt(sum((w_test).^2,1));
load('snap100.mat');
e_pod100 = sqrt(sum((w_test-wr_test_interpPROM).^2,1))./sqrt(sum((w_test).^2,1));
e_dmd100 = sqrt(sum((w_test-wr_test_interpPROM_dmd).^2,1))./sqrt(sum((w_test).^2,1));

figure
width = 5;     % Width in inches
height = 4;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
plot(t,e_pod25,'k-.','LineWidth',1);
plot(t,e_dmd25,'k','LineWidth',1);
plot(t,e_pod50,'b-.','LineWidth',1);
plot(t,e_dmd50,'b','LineWidth',1);
plot(t,e_pod100,'r-.','LineWidth',1);
plot(t,e_dmd100,'r','LineWidth',1);
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$L_2$ error','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'POD-PROM,$N_{snap} = 25$','DMD-PROM,$N_{snap} = 25$',...
    'POD-PROM,$N_{snap} = 50$','DMD-PROM,$N_{snap} = 50$',...
    'POD-PROM,$N_{snap} = 100$','DMD-PROM,$N_{snap} = 100$'},...
    'FontUnits','points','interpreter','latex','FontSize',10,'Location','northwest');
legend boxoff;

clear all;
load('rank2.mat');
e_pod2 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd2 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank4.mat');
e_pod4 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd4 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank6.mat');
e_pod6 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd6 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank8.mat');
e_pod8 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd8 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank10.mat');
e_pod10 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd10 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank12.mat');
e_pod12 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd12 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank14.mat');
e_pod14 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd14 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
load('rank16.mat');
e_pod16 = norm(w_test-wr_test_interpPROM)./norm(w_test);
e_dmd16 = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);

saveas(gcf,'compare1','epsc'); 

figure
width = 5;     % Width in inches
height = 4;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
plot([2,4,6,8,10,12,14,16],[e_pod2,e_pod4,e_pod6,e_pod8,e_pod10,e_pod12,e_pod14,e_pod16],'ko-.','LineWidth',1);
plot([2,4,6,8,10,12,14,16],[e_dmd2,e_dmd4,e_dmd6,e_dmd8,e_dmd10,e_dmd12,e_dmd14,e_dmd16],'rd-','LineWidth',1);

xlabel('rank $r$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$L_2$ error','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'POD-PROM,$N_{snap} = 100$','DMD-PROM,$N_{snap} = 100$'},...
    'FontUnits','points','interpreter','latex','FontSize',10,'Location','northeast');
legend boxoff;

saveas(gcf,'compare2','epsc'); 











