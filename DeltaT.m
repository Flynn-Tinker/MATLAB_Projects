clear all
close all

load Ta_2009.mat

% What is in the Ta data?

Datenum = Ta(:,1);
Delta_T = Ta(:,2);

%dn = datestr(Ta(:,1));
%dn = dn(:,4:6)

figure(1)
plot(Datenum,Delta_T)

%Set Ticks
% labels = datestr(Datenum, 20);
datetick('x','m')
% set(gca, 'XTick', Datenum);
% set(gca, 'XTickLabel', labels);

%Label Axes and Set Title
xlabel('Date','Fontsize',20)
ylabel('Delta T','Fontsize',20)
title('SST 2009','Fontsize',20)
ylim([-1 6]);
set(gca,'Fontsize',20);

% orientation landscape

