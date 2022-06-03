%%% NOTE: il modello funziona e conserva la massa (fa i loop spaziali).
%%% calcola la sezione del flusso come proiezione nella direzione che
%%% collega il punto medio tra due prismi bagnati del tirante della cella
%%% con il potenziale maggiore (tenendo conto della quota del fondo). Il
%%% flusso è definito tramite una funzione alla fine



%%%%%%
%NOTE: 
%%%%%%

% lunghezze: metri  
% tempo:     minuti 

% Assunzione: terreno inizialmente all field capacity --> tutta la pioggia diventa runoff. Porosità efficace : porosità (i.e. contenuto d'acqua a saturazione) - contenuto d'acqua alla field capacity 

% Propagazione dei flussi in 8 direzioni

% DTM: Deve essere presente in formato .asc come da esempio nella stessa cartella dell'eseguibile

% Timestep fissato in minuti (eventualmente implementare il timestep variabile) 

% Nel caso in cui il modello non convergesse o fosse instabile bisogna
% passare ai secondi (facendo attenzione a cambiare le unità di misura coerentemente). Vedere condizione di Courant

% Passare alle ore nel caso in cui i tempi di calcolo fossero troppo lunghi (attenzione a fare le dovute modifice, vedi sopra)

%> Accetta conducibilità idrauliche e spessore di suolo spazialmente variabili

%%%%%%%%%%
%%%% START
%%%%%%%%%%


L=5                           %   <<< INPUT   dimensione lato cella DTM (metri) - può essere diversa da quella nell'header dels ascii del DTM   
timesteps=3000                %   <<< INPUT   numero di timesteps: lunghezza massima della simulazione in MINUTI (o secondi, vedi sopra)
intervallo_salvataggi=30      %   <<< INPUT   ongi quanti timesteps (i.e. minuti) salva gli output (divisore di timesteps)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rain characteristics (volume, duration, time start)           <<< INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rain_depth= 0.1;        % in metri 
rain_duration=120;      % minuti (verificare coerenza con il K e il tempo della simulazione) 
rain_start=5;           % minuti (dopo quanto tempo dall'inizio della simulazione inizia a piovere)

rain_intensity=rain_depth/rain_duration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETRI MODELLO%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

soil_depth=1                                      %   <<< INPUT   profondità dello strato roccioso impermeabile (in metri)
k=0.01                                            %   <<< INPUT   conducibilità idraulica in metri/minuto 
porosity=0.3                                      %   <<< INPUT   porosità (-)
water_content_at_FC=0.1                           %   <<< INPUT   contenuto d'acqua alla capacità di campo (-)
effective_porosity=porosity-water_content_at_FC   %   <<< INPUT   POROSITA' EFFICACE  (equivale all'acqua che può essere drenata) : porosità (i.e. contenuto d'acqua a saturazione) - contenuto d'acqua alla field capacity 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TOPO TOOLBOX %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identifica la root directory e assegna la posizione di topotoolbox (deve essere nella stessa cartella)
root_dir=fileparts(mfilename('fullpath'));
addpath(genpath(strcat(root_dir,'\TopoToolbox-2')));  % assicurarsi che il path sia giusto

DEM = GRIDobj('DTM.asc');    %%% file DTM.asc deve essere in formato ASCII con l'header  
DTM=DEM.Z;


% %%%%%%%%%%%%%%%
% %%% inporta DTM 
% %%%%%%%%%%%%%%%
% 
% fileID = fopen('DTM.asc'); 
% DTM_t = textscan(fileID, '%s', 'Delimiter', '\n','headerLines', 6);
% fclose(fileID);
% DTM_t = DTM_t{1,1};
% DTM=[]
% 
% for i=1:length(DTM_t)
%     DTM(i,:)= str2double( strsplit(DTM_t{i,1},' '));
% end
% %DTM(:,end)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% allocazione della memoria  - inizializzazione delle variabili
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_row=length(DTM(:,1));
num_col=length(DTM(1,:));
 
num_row=num_row+2;
num_col=num_col+2;

temp=nan(num_row,num_col);
temp(2:end-1,2:end-1)=DTM;
DTM_bordo=temp;
DTM_bordo(DTM_bordo==-9999)=NaN;

%conducibilità idraulica (eventualmente spazialmente variabile)
if isfile('K.asc')
    K = GRIDobj('K.asc');               % i valori di K nel file K.asc (se presente) vengono usati per la K spazialmente esplicita 
else
    K = k*ones(num_row,num_col);        % se K non è assegnata come sopra, viene usato il valore di k omogeneo dato in input all'inizio 
end 


%profondità suolo  (eventualmente spazialmente variabile) - come per K, vedi sopra
if isfile('SOIL_DEPTH.asc')
    SOIL_DEPTH = GRIDobj('SOIL_DEPTH.asc');               
else
    SOIL_DEPTH = soil_depth*ones(num_row,num_col);        
end 


flusso_su=zeros(num_row,num_col);
flusso_giu=zeros(num_row,num_col);
flusso_dx=zeros(num_row,num_col);
flusso_sx=zeros(num_row,num_col);
flusso_su_sx=zeros(num_row,num_col);
flusso_giu_sx=zeros(num_row,num_col);
flusso_su_dx=zeros(num_row,num_col);
flusso_giu_dx=zeros(num_row,num_col);

Array_portate=zeros(1,2);
distance_to_outlet=zeros(num_row,num_col);
dist_to_outlet_matrix=[];

u_max=timesteps/intervallo_salvataggi; % deve essere intero
u=1;                                   % contatore

rain_dept_dt=zeros(timesteps,1);
Q=zeros(timesteps,1);
volume_invasato=zeros(timesteps,1);
SAT_SOIL_THICK=nan(num_row,num_col,u_max);
OUTFLOW=nan(num_row,num_col,u_max);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   forse vuole il DTM quadrato
% %%%flow accumulation  --- definisce i canali di drenaggio   %%%  !!! IN QUESTA VERSIONE DEL MODELLO IL CANALE NON VIENE ESTRATTO IN QUESTO MODO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,Y] = meshgrid([L/2:L:max_size*L-L/2],[L/2:L:max_size*L-L/2]);
% [flow_accumulation,flowdir,slope,runs] = wflowacc(X,Y,DTM_square,'type' ,'single');       %%%%%   <----------- DA CAMBIARE possiblie mettere 'single' 'multi'  ( algoritmo usato per calcolare le direzioni di drenaggio)
% celle_canale(flow_accumulation>2500)=1;                                                    %%%%% (2500)  <----------- DA CAMBIARE - numero di celle necessario avere a monte per definire un pixel canale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TOPO TOOLBOX %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()

% visualizza il dem
subplot(2,4,1)
imagesc(DEM)  

% calcola gradiente
G = gradient8(DEM);
subplot(2,4,2)
imageschs(DEM,G,'ticklabel','nice','colorbarylabel','Slope [-]','caxis',[0 1])

% fill sinks 
DEMf = fillsinks(DEM);             

% flow accumulation 
FD = FLOWobj(DEMf);
A  = flowacc(FD);
subplot(2,4,3)
imageschs(DEM,dilate(sqrt(A),ones(5)),'colormap',flowcolor,...
    'colorbarylabel','Flow accumulation [sqrt(# of pixels)]',...
    'ticklabel','nice');

% drainage basin
DB = drainagebasins(FD);
DB = shufflelabel(DB);


% area and plot 
nrDB = numel(unique(DB.Z(:)))-1; % nr of drainage basins
STATS = regionprops(DB.Z,'PixelIdxList','Area','Centroid');
imageschs(DEM,DB,'colorbar',false,'ticklabel','nice');
hold on
for run = 1:nrDB
    if STATS(run).Area*DB.cellsize^2 > 10e6
        [x,y] = ind2coord(DB,...
            sub2ind(DB.size,...
            round(STATS(run).Centroid(2)),...
            round(STATS(run).Centroid(1))));
        text(x,y,...
            num2str(round(STATS(run).Area * DB.cellsize^2/1e6)),...
            'BackgroundColor',[1 1 1]);
    end
end
hold off


% flow distance 
D = flowdistance(FD);
D = D/1000;
subplot(2,4,4)
imageschs(DEM,D,'ticklabel','nice','colorbarylabel','Flow distance [km]')
  

% frquency of distance from outlet
[~,IX] = max([STATS.Area]);
subplot(2,4,5)
histogram(D.Z(DB.Z == IX),'Normalization','pdf');
xlabel('Distance to outlet [km]');
ylabel('fraction of area');


% calculate flow accumulation
A = flowacc(FD);
% Note that flowacc returns the number of cells draining in a cell
W = A>500;   % <<< INPUT  area cumulata necessaria per inizializzare un canale

% create an instance of STREAMobj
S = STREAMobj(FD,W);
subplot(2,4,6)
plot(S);
axis image


% extract the largest subnetwork
S = klargestconncomps(S,1);
%figure()
%plot(S); axis image

% plot distance vs elevation
subplot(2,4,7)
plotdz(S,DEM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% calcola caratteristiche morfologiche del bacino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catchment_area=sum(sum(isnan(DTM_bordo)==0))*L^2  / 1000^2;    % area bacino in km^2
slope_average=mean(mean(G.Z(isnan(G.Z)==0)));
min_catch_elevation=min(DEM);
max_catch_elevation=max(DEM);
distance_to_outlet(2:end-1 , 2:end-1)=D.Z;  % in km

%caratteristiche flusso in stato canale
manning= 0.05;          %% <<< INPUT   manning =  0.04 - 0.07   per torrenti di montagna molto irregolari    m/s                                                             %<----- scabrezza da utilizzare nel trasporto in stato canale
average_cahannel_slope = (max_catch_elevation-min_catch_elevation)/(max(max(D.Z))*1000);
flow_velocity_channel  = 0.317*(catchment_area)^0.125*average_cahannel_slope^0.375/manning^0.75; % quantificazione empirica velocità media della  nei canali


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INIZIO COMPUTAZIONI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z_impervious = DTM_bordo - SOIL_DEPTH;   % profondità dello strato impermeabile
potenziale   = Z_impervious;  

rain_dept_dt(rain_start:rain_start+rain_duration-1)=rain_intensity;

tic

for t=1:timesteps

     potenziale=potenziale+rain_dept_dt(t)/effective_porosity; % aggiunge la pioggia 
    
    
      % calcolo dei flussi (vedere funzioni alla fine)

      for i=2:num_row-1                                 
          
        for j=2:num_col-1
                       
            flusso_su(i,j)     =  flusso_fun(i,j,-1,0,K,potenziale,Z_impervious,L);   
            flusso_giu(i,j)    =  flusso_fun(i,j,1,0,K,potenziale,Z_impervious,L);   
            flusso_dx(i,j)     =  flusso_fun(i,j,0,1,K,potenziale,Z_impervious,L);
            flusso_sx(i,j)     =  flusso_fun(i,j,0,-1,K,potenziale,Z_impervious,L);
            
            flusso_su_dx(i,j)  =  flusso_fun(i,j,-1,1,K,potenziale,Z_impervious,L);
            flusso_su_sx(i,j)  =  flusso_fun(i,j,-1,-1,K,potenziale,Z_impervious,L);
            flusso_giu_dx(i,j) =  flusso_fun(i,j,1,1,K,potenziale,Z_impervious,L);
            flusso_giu_sx(i,j) =  flusso_fun(i,j,1,-1,K,potenziale,Z_impervious,L);

        end
        
      end
             
            
      potenziale = potenziale + ( flusso_su + flusso_giu +   flusso_dx + flusso_sx +  flusso_giu_sx + flusso_giu_dx + flusso_su_sx + flusso_su_dx) / effective_porosity / (L^2) ;
      
      potenziale(potenziale<=Z_impervious)=Z_impervious(potenziale<=Z_impervious); 
    
      pot=potenziale-Z_impervious;  % potenziale relativo al fondo 
        
      volume_invasato(t)=sum(sum(pot(isnan(pot)==0)))*effective_porosity*L^2;

      delta_pot = pot - SOIL_DEPTH;
      delta_pot(delta_pot<0)=0;
%     delta_pot(isnan(delta_pot))=0;
        
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% convoluzione per calcolare la portata nella sezione di chiusura
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      Q_pixel=delta_pot(delta_pot>0)*L^2*effective_porosity;
      dist_to_outlet_matrix=distance_to_outlet(delta_pot>0);
      delta_t_pixel =   dist_to_outlet_matrix*1000/flow_velocity_channel / 60;    %% tempo necessario per raggiungere l'outlet in minuti
        
      delta_t_pixel(delta_t_pixel<=0.5)=1;
      delta_t_pixel=round(delta_t_pixel);
             
      Array_portate(:,2)=Array_portate(:,2)-1;
      Q(t)= sum(Array_portate( Array_portate(:,2)==0,1)); 
      Array_portate(Array_portate(:,2)==0,:)=[];                           
      Array_portate=[Array_portate;[Q_pixel,delta_t_pixel]];
        
      
      potenziale = potenziale - delta_pot; % rimuove l'acqua emersa
        
        
      if ismember(t,[1:intervallo_salvataggi:timesteps])
       SAT_SOIL_THICK(:,:,u)=pot;    
       OUTFLOW(:,:,u)=delta_pot*effective_porosity;   % metri / minuto
       percent_completed = round(u/u_max * 100)
       u=u+1;
      end
        
        
        
end
 
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% calcola un pò di cose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q=Q/60*1000;                        % converte la portata da metri cubi al minuto in litri al secondo
[Q_max,peak_time] = max(Q)          % portata in litri al secondo e tempo in minuti


% salva i risultati nella stessa directory dell'eseguibile
save('risultati.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(Q)
title('Idrogramma')
xlabel('Tempo (minuti)') 
ylabel('Portata (litri/secondo)')  


%%% Numero di celle con outflow > 1 e loro distanza dall'outlet
outflow_number_cells=[];
outflow_number_cells_time=[];

for i=1:u-1
    outflow_number_cells(i)=sum(sum((OUTFLOW(:,:,i)>0)));
    outflow_number_cells_time(i)= i * intervallo_salvataggi;
end

outflow_distance_to_outlet = NaN(max(outflow_number_cells),u-1);
for i=1:u-1
    temp = distance_to_outlet(OUTFLOW(:,:,i)>0);
    outflow_distance_to_outlet(1:length(temp),i)=temp(:);
end

figure()
subplot(1,2,1)
plot([1:1:u-1]*intervallo_salvataggi,nanmean(outflow_distance_to_outlet))
xlabel('time (minutes)')
ylabel('average distance to outlet of outflow cells')

subplot(1,2,2)
plot([1:1:u-1]*intervallo_salvataggi,outflow_number_cells)
xlabel('time (minutes)')
ylabel('number of outflow cells')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% crea la GIF dello spessore di suolo saturo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Saturated_soil_thickness.gif';
max_colorbar=max(max(max(SAT_SOIL_THICK(:,:,:))));

for w=1:u-1

    k=pcolor(SAT_SOIL_THICK(:,:,w));
    set(k,'EdgeColor', 'none');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    axis image;
    caxis([0 max_colorbar])
    c=colorbar;
    c.Label.String = 'Thickness of saturated soil (m)'
    legend(['Time (hours) = ' num2str(round(w*intervallo_salvataggi/ 60))]);
    
    drawnow 
    
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if w == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% crea la GIF dell OUTFLOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Outflow.gif';
max_colorbar=max(max(max(OUTFLOW(:,:,:).^0.5)));  %%%% <--- eventualmente cambiare ,il valore è trasformato per metterlo meglio in evidenza

for w=1:u-1

    k=pcolor(OUTFLOW(:,:,w).^0.5);     %%%% <---- eventualmente cambiare , il valore è trasformato per metterlo meglio in evidenza 
    set(k,'EdgeColor', 'none');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    axis image;
    caxis([0 max_colorbar])
    c=colorbar;
    c.Label.String = 'Outflow (m/min)'
    legend(['Time (hours) = ' num2str(round(w*intervallo_salvataggi/ 60))]);
    
    drawnow 
    
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if w == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FUNCTIONS%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%calcola flussi 

function flusso = flusso_fun(i,j,a,b,K,potenziale,Z_impervious,L)

  if a~=0 && b~=0
    dist= sqrt(2)*L;
  else
    dist=L;
  end


if isnan(potenziale (i+a,j+b)) || isnan(potenziale (i,j))
    flusso=0;
elseif potenziale(i+a, j+b)<=Z_impervious(i+a, j+b) && potenziale (i,j)< Z_impervious(i+a,j+b)
    flusso=0;
    
elseif potenziale(i, j)<=Z_impervious(i, j) && potenziale (i+a,j+b)< Z_impervious(i,j)
    flusso=0;
    
elseif potenziale(i, j)<=Z_impervious(i, j) && potenziale(i+a, j+b)<=Z_impervious(i+a, j+b) % forse inutile ???
    flusso=0;
    
elseif potenziale(i,j) >potenziale(i+a,j+b) 
    flusso = (2* K(i+a,j+b)* K(i,j))/(K(i+a,j+b)+ K(i,j))    *   L/2  *  ( potenziale(i+a, j+b) - potenziale(i,j) )  /   sqrt(dist^2   + ( (potenziale(i, j)+ Z_impervious(i, j))/2   -  (potenziale(i+a, j+b)  + Z_impervious(i+a, j+b))/2)^2)    *  (  potenziale(i, j)     - max(Z_impervious(i, j), Z_impervious(i+a, j+b))   )   *  cos(  atan(  abs( ((potenziale(i,j)+Z_impervious(i,j))/2  -  (potenziale(i+a,j+b)+Z_impervious(i+a,j+b))/2)/dist)  )  )     ;
else
    flusso = (2* K(i+a,j+b)* K(i,j))/(K(i+a,j+b)+ K(i,j))    *   L/2  *  ( potenziale(i+a, j+b) - potenziale(i,j) )  /   sqrt(dist^2   + ( (potenziale(i, j)+ Z_impervious(i, j))/2   -  (potenziale(i+a, j+b)  + Z_impervious(i+a, j+b))/2)^2)    *  (  potenziale(i+a, j+b) - max(Z_impervious(i, j), Z_impervious(i+a, j+b))   )   *  cos(  atan(  abs( ((potenziale(i,j)+Z_impervious(i,j))/2  -  (potenziale(i+a,j+b)+Z_impervious(i+a,j+b))/2)/dist)  )  )     ;   
    

    
 
%%%%%% versione originale

%     if potenziale(i,j) >potenziale(i+a,j+b)
%         flusso = K(i+a,j+b)    *   L/2  *  ( potenziale(i+a, j+b) - potenziale(i,j) )  /   sqrt(dist^2+ ( (potenziale(i, j)+ Z_impervious(i, j))/2   -(potenziale(i+a, j+b)  + Z_impervious(i+a, j+b))/2)^2)    *  (  potenziale(i, j)     - max(Z_impervious(i, j), Z_impervious(i+a, j+b))      ) ;
%     else
%         flusso = K(i,j)        *   L/2  *  ( potenziale(i+a, j+b) - potenziale(i,j) )  /   sqrt(dist^2+ ( (potenziale(i, j)+ Z_impervious(i, j))/2   -(potenziale(i+a, j+b)  + Z_impervious(i+a, j+b))/2)^2)    *  (  potenziale(i+a, j+b) - max(Z_impervious(i, j), Z_impervious(i+a, j+b))      ) ;    
%     end


%%%%% versione modificata che tiene conto anche del tirante di valle (non
%%%%% conserva la massa

%     flusso = (2* K(i+a,j+b)* K(i,j))/(K(i+a,j+b)+ K(i,j))    *   L/2  *  ( potenziale(i+a, j+b) - potenziale(i,j) )  /   sqrt(dist^2   + ( (potenziale(i, j)+ Z_impervious(i, j))/2   -  (potenziale(i+a, j+b)  + Z_impervious(i+a, j+b))/2)^2)    *  (   (potenziale(i,j) + potenziale(i+a,j+b))/2  -  (Z_impervious(i,j) + Z_impervious(i+a,j+b))/2     )  ;% *  cos(  atan(  abs( ((potenziale(i,j)+Z_impervious(i,j))/2  -  (potenziale(i+a,j+b)+Z_impervious(i+a,j+b))/2)/dist)  )  )     ;

end


end



