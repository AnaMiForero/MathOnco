%--------------------------------------------------------------
%---------------------------------------------------------------
%                         Data
%---------------------------------------------------------------
%---------------------------------------------------------------
warning off
path = '/Users/ana/Desktop/Exam1_IMO2/';
df = readtable(strcat(path,'Data.csv'), 'Delimiter',',');
df1 = table2array(df);
con = df1(:,1);
uptake_data = df1(:,2);


%--------------------------------------------------------------
%---------------------------------------------------------------
%                         Main  
%---------------------------------------------------------------
%---------------------------------------------------------------

params_guess = [0.5,0.5];
lb = [0,0];
ub = [1,1];
[solution,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(params,x) Model(params,x),params_guess,con,uptake_data, lb, ub)

% %Plots

%Plot for fit
f1 = figure(1);
plot(con,uptake_data,'bo', con, Model(solution,con),'g--', 'MarkerSize',10, 'LineWidth',3);
set(gca,'fontsize',18)
xlabel('\textbf{Concentration [$\mu g/ \mu m^3$]}', 'Interpreter','latex'); 
ylabel('\textbf{Uptake rate per cell [$\mu g/ s$]}', 'Interpreter','latex');
legend({'Data','Fit'}, 'Location','best', 'Position',[0.65 0.25 0.2 0.15],'Interpreter','latex','FontSize',14); %[left bottom width height]
grid on
exportgraphics(f1,'fit.pdf','BackgroundColor','none');

% Plot for L parameter behavior
L_1 = [1.3, 1, 0.8];
k_1 = 0.5;
con1 = 0:200;

f2 = figure(2);
plot(con1,Model([L_1(1),k_1], con1),'b-', con1, Model([L_1(2),k_1], con1),'g-',con1, Model([L_1(3),k_1],con1),'m-', 'MarkerSize',10, 'LineWidth',3);
title('$\kappa = 0.5$','Interpreter','latex')
xlabel('\textbf{Concentration [$\mu g/ \mu m^3$]}', 'Interpreter','latex'); 
ylabel('\textbf{Uptake rate per cell [$\mu g/ s$]}', 'Interpreter','latex');
legend({'$\mathcal{L} = 1.3$','$\mathcal{L} =1.0$','$\mathcal{L} = 0.8$'}, 'Location','best', 'Position',[0.65 0.25 0.2 0.15],'Interpreter','latex','FontSize',14);
set(gca,'fontsize',18)
grid on
exportgraphics(f2,'kconst.pdf','BackgroundColor','none');


%Plot for k parameter behavior
L_2 = 0.5;
k_2 = [1, 0.3, 0.1]; 

f3 = figure(3);
plot(con1,Model([L_2,k_2(1)], con1),'r-', con1, Model([L_2,k_2(2)], con1),'k-',con1, Model([L_2,k_2(3)],con1),'c-', 'MarkerSize',10, 'LineWidth',3);
title('$\mathcal{L} = 0.5$','Interpreter','latex')
legend({'$\kappa = 1$','$\kappa = 0.3$','$\kappa =0.1$'}, 'Location','best', 'Position',[0.65 0.25 0.2 0.15],'Interpreter','latex','FontSize',14);
set(gca,'fontsize',18)
xlabel('\textbf{Concentration [$\mu g/ \mu m^3$]}', 'Interpreter','latex'); 
ylabel('\textbf{Uptake rate per cell [$\mu g/ s$]}', 'Interpreter','latex');
grid on
exportgraphics(f3,'Lconst.pdf','BackgroundColor','none');

 
% %--------------------------------------------------------------
% %---------------------------------------------------------------
% %                         Functions
% %---------------------------------------------------------------
% %---------------------------------------------------------------


function y = Model(params,x)
    L = params(1);
    k = params(2);
    y = L./(1+exp(-k*(x - 50)));
end

