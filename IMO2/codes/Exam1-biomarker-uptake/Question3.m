
concentration = [1, 10, 25, 50, 75, 100, 200, 300, 400, 500];      % [micrograms/micrometer^2]
export = [];

for g = 1:size(concentration,2)
rng(7)
%----------- COMPUTATIONAL DOMAIN ------------------%

xmin = -500; xmax = 500; ymin = xmin; ymax = xmax; 
hg = 6;                                     % width of the grid 
[xx,yy] = meshgrid(xmin:hg:xmax,ymin:hg:ymax);
[Ngy, Ngx] = size(xx);
area = (Ngx-1)*(Ngy-1);                     % area of the lattice 

%----------- PARAMETERS------------------%

time = 2;                                       % [h]
dt = 0.6;                                       % time step [s]
Niter = round(time*3600/dt);                    % number of iterations 
diff  =10;                                      % diffusion coefficient 
stability = diff*dt/(hg^2);                     % stability condition  
Rcell= hg;                                      % cell radius [micrometer]
areacell = pi*Rcell^2;                          % cell area [micrometer]
Ncells = round(0.35*area/areacell);             % number of cells
cellMaturation = 20*3600;                       % maturation age (20 h) [s]
ageMat = ones(Ncells,1)*cellMaturation;         % maturation age for each cell
maxNeigh = 6;                                   % number of neighbours needed before division can happen
nu = 120;                                       % medium viscosity 
L = 0.0201;                                     % [micrograms/s]   
kappa = 0.0937;                                 % [micrometer^2/microgram]   
stif = 100;                                     % tumor cell springs stiffness
step = 1+floor(Rcell/hg);
uptake_cond = 1;                                % if 0, min(L*dt, biomark(Ny+iy+1,Nx+ix+1)). If 1, L*dt/(1+exp(-kappa*(biomark(x,y) - 50))) 
to_save = 0;                                    % if 1, save csv. if 0 not to save.

%----------- INITIAL CONDITIONS ---------%

cells = 2*(xmax-2*Rcell)*(rand(Ncells,2) - 0.5);    % cells' coordinates
age = rand(Ncells,1)*cellMaturation;                % initialize ages
cell_biomark = zeros(Ncells,1);                     % biomarker concentration inside cell
gamma0 = concentration(g);                         % Initial value of biomarker [micrograms]
biomark = gamma0*ones(Ngy,Ngx); 

% set zero value of biomarker inside cells
for ii = 1:Ncells
        angle = linspace(0,2*pi,50);  
        xc = cells(ii,1)+(Rcell*cos(angle));
        yc = cells(ii,2)+(Rcell*sin(angle));
        in = inpolygon(xx,yy,xc,yc);
        biomark(in) = 0; %This works because "in" is a logical array
end

%-------------IF VISUALIZATION-----------------%

plot = 1; %1 to plot, 0 to not plot. 
concentration_cond = 75; %choose the concentration you want to plot 
if (plot == 1 && gamma0 == concentration_cond)
    figure
end

%-------------SIMULATION------------------------%


for iter = 0:Niter

    age = age + dt; % tumor cell age 
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % tumor cell neighbours
    %%%%%%%%%%%%%%%%%%%%%%%
    
    cellNeigh = zeros(Ncells,1);  
    
    for ii = 1:Ncells
       for jj = 1:Ncells
         dx = cells(ii,1)-cells(jj,1); 
         dy = cells(ii,2)-cells(jj,2);
         dxy = sqrt(dx^2+dy^2);     
         if (dxy <= 3.5*Rcell) && (dxy > 0)
           cellNeigh(ii) = cellNeigh(ii) + 1;
         end
       end 
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % tumor cell division
    %%%%%%%%%%%%%%%%%%%%%%%
    
    for cc = 1:Ncells
          if (age(cc) > ageMat(cc)) && (cellNeigh(cc) < maxNeigh)
            biom_mother = cell_biomark(cc);
            Ncells = Ncells + 1;  % cell division 
            theta = rand*2*pi; 
            cells(Ncells,1:2) = cells(cc,1:2)+0.5*Rcell*[cos(theta),sin(theta)];
            age(Ncells) = 0;
            age(cc) = 0;
            ageMat(Ncells) = ageMat(cc)+2*(rand-0.5)*0.2*ageMat(cc); 
            %Division of biomarker:
            cell_biomark(Ncells) = 0.5*biom_mother;
            cell_biomark(cc) = 0.5*biom_mother;
         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % Forces
    %%%%%%%%%%%%%%%%%%%%%%%
    % tumor cell-tumor cell repulsive and adhesive forces
    RepForce = zeros(Ncells,2); %Respulsive force
    for ii = 1:Ncells-1
        for jj = ii+1:Ncells
            dx = cells(ii,1)-cells(jj,1);
            dy = cells(ii,2)-cells(jj,2);
            dxy = sqrt(dx^2+dy^2);
            if (dxy<2*Rcell)&&(dxy>0)
                RepForce(ii,1:2)=RepForce(ii,1:2)+stif*(2*Rcell-dxy)*[dx,dy]/dxy;
                RepForce(jj,1:2)=RepForce(jj,1:2)-stif*(2*Rcell-dxy)*[dx,dy]/dxy;
            end

        end
    end
    cells = cells+dt*RepForce/nu;

    % clean up 
    % check outside tumor cells
    ind = find((cells(:,1)>xmin)&(cells(:,1)<xmax)&...
             (cells(:,2)>ymin)&(cells(:,2)<ymax));
    cells = cells(ind,:);
    Ncells = size(cells,1);
    RepForce = RepForce(ind,:);
    age = age(ind,:);
    ageMat = ageMat(ind,:);
    cell_biomark = cell_biomark(ind,:); 

    %%%%%%%%%%%%%%%%%%%%%%%
    % Biomarker diffussion
    %%%%%%%%%%%%%%%%%%%%%%%
    % Neumann boundary condition in a central difference approximation.

    biomark(1,:) = biomark(2,:);
    biomark(Ngy,:) = biomark(Ngy-1,:);
    biomark(:,1) = biomark(:,2);
    biomark(:,Ngx) = biomark(:,Ngx-1); 
    
   
    % biomarker diffusion
    biomarkL = biomark(1:Ngy-2,2:Ngx-1);  % diffusion from left
    biomarkR = biomark(3:Ngy,2:Ngx-1);  % diffusion from right
    biomarkT = biomark(2:Ngy-1,3:Ngx);  % diffusion from top
    biomarkB = biomark(2:Ngy-1,1:Ngx-2);  % diffusion from bottom

    biomark(2:Ngy-1,2:Ngx-1) = biomark(2:Ngy-1,2:Ngx-1) +(diff*dt/(hg*hg))*...
                          (biomarkL+biomarkR+biomarkT+biomarkB-4*biomark(2:Ngy-1,2:Ngx-1));

    %%%%%%%%%%%%%%%%%%%%%%%
    % Biomarker uptake
    %%%%%%%%%%%%%%%%%%%%%%%
     for ii = 1:Ncells %Iterate over cells
       Nx = 1+floor((cells(ii,1)-xmin)/hg); % closest grid point to the cell
       Ny = 1+floor((cells(ii,2)-ymin)/hg);
        
       for ix = -step:step  %Iterate over grid points near the cell
        for iy = -step:step  
            if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)
                ixy = sqrt((cells(ii,1)-(xmin+(Nx+ix)*hg))^2+(cells(ii,2)-(ymin+(Ny+iy)*hg))^2);
                if (ixy<Rcell) %If grid point is inside cell
                    if uptake_cond == 0
                        uptake = min(L*dt, biomark(Ny+iy+1,Nx+ix+1));
                        biomark(Ny+iy+1,Nx+ix+1) = max(0, biomark(Ny+iy+1,Nx+ix+1) - uptake);
                    elseif uptake_cond == 1
                        uptake = L*dt/(1+exp(-kappa*(biomark(Ny+iy+1,Nx+ix+1) - 50))) ;
                        biomark(Ny+iy+1,Nx+ix+1) = max(0, biomark(Ny+iy+1,Nx+ix+1) - uptake);
                    end

                 %Store cell biomark uptake per each cell
                 cell_biomark(ii) = cell_biomark(ii) + uptake;
                end 
                
            end 
        end 
       end
     end 


    %%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%
    if (plot == 1 && gamma0 == concentration_cond)
        if (mod(iter,100) == 0)
            clf
            ax1 = axes;
            contourf(ax1,xx,yy,biomark,'edgecolor','none');
            axis equal
            ax2 = axes;
            %plot cells 
            for ii = 1:Ncells
                angle = linspace(0,2*pi,50);  
                xc = cells(ii,1)+(Rcell*cos(angle));
                yc = cells(ii,2)+(Rcell*sin(angle));
                patch(ax2,xc,yc,cell_biomark(ii))
            end
            axis equal
            axis([xmin-0.5,xmax+0.5,ymin-0.5,ymax+0.5])
    
            % Hide the top axes
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];
            colormap(ax1,"cool")
            colormap(ax2,flipud(hot))
    
    
            % colorbars
            cb1 = colorbar(ax1,'Position',[0.1 0.1 0.05 0.815]); % Position [left bottom width height]
            cb2 = colorbar(ax2,'Position',[0.85 0.1 0.05 0.815]);
            cb1.Label.String = 'Biomarker concentration';
            cb2.Label.String = 'Uptake concentration';
            cb1.Label.FontSize = 14;
            cb2.Label.FontSize = 14;
            cb1.Label.FontWeight = 'bold';
            cb2.Label.FontWeight = 'bold';
            title(ax1, ['Iteration=',num2str(iter),'   Time=',num2str(round(iter*dt/3600,2)),' [h]' ,'   Ncells=',num2str(Ncells) ], "FontSize",14)
    
        pause(0.1)
        end
    end
end
disp(strcat('Concentration: ',num2str(gamma0), ' done'))
export = [export, cell_biomark];
end

if to_save == 1
    to_csv = array2table(export);
    to_csv.Properties.VariableNames(1:size(concentration,2)) = {'con_1', 'con_10', 'con_25', 'con_50', 'con_75', 'con_100', 'con_200', 'con_300', 'con_400', 'con_500'};
    writetable(to_csv, strcat('uptake_condition_',num2str(uptake_cond),'.csv'))
end

