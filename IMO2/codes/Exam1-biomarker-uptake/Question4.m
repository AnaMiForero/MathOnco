
warning off
path = '/Users/ana/Desktop/Exam1_IMO2/';

dfa = readtable(strcat(path,'uptake_condition_0.csv'), 'Delimiter',','); % min(L*dt, biomark(Ny+iy+1,Nx+ix+1))
df1 = table2array(dfa);
total_cells_0 = size(df1,1);


dfb = readtable(strcat(path,'uptake_condition_1.csv'), 'Delimiter',','); % L*dt/(1+exp(-kappa*(biomark(x,y) - 50)))
df2 = table2array(dfb);
total_cells_1 = size(df2,1);

percentage = @(con, total_cells) sum(con >= 200.0)*100/total_cells; %con is the concentration of biomarker inside cell

percentage_0 = [];
percentage_1 = [];
for i = 1:size(df1,2)
    percentage_0 = [percentage_0, percentage(df1(:,i), total_cells_0)];
    percentage_1 = [percentage_1, percentage(df2(:,i), total_cells_1)];
end

y = [percentage_0;percentage_1];

%Plots
figure(1)
histogram(df2(:,4),10, 'FaceColor','#EDB120');
title('\textbf{Initial concentration: $50$}','Interpreter','latex')
xlabel('\textbf{Concentration inside cell [$\mu g/ \mu m^2$]}', 'Interpreter','latex'); 
ylabel('\textbf{Number of cells}', 'Interpreter','latex');
set(gca,'fontsize',17)

figure(2)
histogram(df2(:,5),10, 'FaceColor','#D95319');
title('\textbf{Initial concentration: $75$}','Interpreter','latex')
xlabel('\textbf{Concentration inside cell [$\mu g/ \mu m^2$]}', 'Interpreter','latex'); 
ylabel('\textbf{Number of cells}', 'Interpreter','latex');
set(gca,'fontsize',17)

figure(3)
h = bar(1:size(df1,2), y, 1);
set(h(1),'FaceColor','#7E2F8E')
set(h(2),'FaceColor','#77AC30')
hold on 
plot(0:size(df1,2)+1,75*ones(size(df1,2)+2,1), 'k--', 'LineWidth',2)
axis([0 size(df1,2)+1 0 110])
ylabel('\textbf{\% of cells with concentration inside $\ge 200$}', 'Interpreter','latex'); 
xlabel('\textbf{Initial concentration}', 'Interpreter','latex'); 
set(gca,'xticklabel',{'1', '10', '25', '50', '75', '100', '200', '300', '400', '500'}, 'fontsize',13)
leg = legend({'\textbf{constant: $\mathcal{L}$}','$\mathcal{L}/[1+e^{-\kappa(\gamma(\textbf{x},t) - 50)}]$'}, ...
                    'Interpreter','latex','Location','northeastoutside','EdgeColor','none');
t = title(leg,'\textbf{\underline{Uptake rates}}', 'Interpreter','latex');
