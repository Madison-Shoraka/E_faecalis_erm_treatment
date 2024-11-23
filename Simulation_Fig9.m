%This code takes in a 2D tiff biofilm image (y,z), analyzes initial
%structure, runs antibiotic treatment simulation, then analyze post
%antibiotic treatment simulation results using same parameters as
%pre-treatment analysis


% structure finding algorithm, IC condtions

%TWO modifiable parameters pre simulation. Modify them to properly read in
%the tiff image. Due to biological variablility, staining strength can change from
%day to day. Further rigid base layer must be properly removed for
%clustering algorithm to work
threshold=230; %above what grey scale do we want to use, change this depending of strength of expiremental staining
density_threshold=.25; %modify between .25-.35 to such that it is approx inflection point of density_vector plot

connectivity=26; %connectivity parameter: can be 6,18,26
threshol_vol=10; %looking for strucutres above a certain number of voxels

%read in tiff files (each slice needs to be a seperate file)
files = dir('*.tif');
%creates binary matrix taking out grey values below threshold
k=3; %which of your tiff files do you want to read  
t = Tiff(files(k).name,'r'); %reads in tiff file
A = read(t);
u=double(A);

%condenses image
myMeanFunction = @(block_struct) mean(block_struct.data);
blockMeans = blockproc(u, [1 4], myMeanFunction);
blockMeans=blockMeans';
blockMeans = blockproc(blockMeans, [1 4], myMeanFunction);
blockMeans=blockMeans''';
u=blockMeans;

sizes=size(u);

%binarizes image
u=(u>threshold);
connected = bwareaopen(u,threshol_vol,4); %gets rid of strucutres below volume threshold
%cuts out rigid baselayer 
break_param=0;
density_vector=[];
for i=1:sizes(1)
    density_vector=[density_vector,sum(sum(u(i,:)))/(sizes(2))];
    if break_param==0
        if sum(sum(u(i,:)))/(sizes(2))<= density_threshold
            break_param=i; 
        end 
    end 
end
bacteria_cut=connected(break_param:end,:);
break_param
pts=[];
marker=1;

%density based clustering
T_d=15; %note, this parameter is equal to the mu_d in model
minpts=50; 
sizes=size(bacteria_cut);    
 
for i=1:sizes(1)
    for j=1:sizes(2)
        if bacteria_cut(i,j)>0
            pts(marker,2)=sizes(1)-i;
            pts(marker,1)=j;
            marker=marker+1;
        end
    end
end
        
idx = dbscan(pts,T_d,minpts);
meanx=[];
meany=[];
uni=unique(idx);
        
sums_specific=[];
for num=1:length(uni)
    if uni(num)~=-1
        M=(idx==uni(num));
        B=double(M);
        new_pts=B.*pts;
        new_pts( ~any(new_pts,2), : ) = [];
        sums_specific=[sums_specific, sum(sum(B))];
       
            meanx(uni(num)) = mean(new_pts(:,2));
            meany(uni(num)) = mean(new_pts(:,1));
  
    end 

end 
%distances between structures
dist=pdist([meanx;meany]');
sums_specific   %volumes of the structures

figure(3)
hold on
gscatter(pts(:,2),pts(:,1),idx);
plot(meanx,meany,'kx',...
'MarkerSize',15,'LineWidth',3) 
hold off

figure(4)
hold on
plot(density_vector)
hold off

%create IC from tiff file
sizes=size(u); 
bacteria=zeros(sizes(1)+40,sizes(2)+40);
bacteria(20:sizes(1)+19,20:sizes(2)+19)=rot90(u,2);
sizes=size(bacteria);

antibiotics=zeros(size(bacteria));
stress_mat=zeros(size(bacteria));
L1=5; %distance at which less dense cells grow

%populate less dense cells that are thresholded from above
for i=20:sizes(1)-20
    for j=20:sizes(2)-20
        
        if mod(i,L1)==0 && mod(j,L1)==0
            if bacteria(i,j)==0
                bacteria(i,j)=1;
            end 
        end
        
    end 
end 

%connection algorithm to define where initial stress is
BW2 = bwareaopen(bacteria,threshol_vol,6); %gets rid of strucutres below volume threshold
L=bwlabeln(BW2,connectivity); %labels/finds connected strucutures
z=bwconncomp(BW2,connectivity); %stores infor about connected strucutures
num_objects=z.NumObjects;
sticky=L>0;
figure(2)
imagesc(bacteria)
bacteria=double(bacteria);

%MAIN SCRIPT-> Run model
domain_x=[0,1];
domain_y=[0,1];

x_pts=sizes(2);
y_pts=sizes(1);
antibiotics(:,x_pts-20)=ones(size(antibiotics(:,1)))*.7392; 
x = linspace(domain_x(1),domain_x(2),x_pts); % x-values
y = linspace(domain_y(1),domain_y(2),y_pts); % y-values
[X,Y] = meshgrid(x,y); % 2d grid


T_g=1200; %number of antibiotics causing cell to stop growing
T_d=10000;
T_a=66; %antibiotics necessary for bacteria to become sticky
g_max=.028; %maximum chance of bacteria doubling
agg_prob=0.02; %probability a nearby cell will become sticky if next to a sticky cell
death_prob=g_max/3; %probability of bacteria dying a natural death
signal_stress_prob_max=1; %maximize for bacterial survival
stress_distance=15; %mu_d
T=1; %end time
dt=.000007; %everything particle has a 1 chance of moving, dt_a
nt=T/dt; %time steps
t=0;
neighborhood_L1=[];
neighborhood_L1_nobounds=[];
counter=1;
L1=5; %distance at which less dense cells grow

for i=-L1:L1
    for j=-L1:L1
        if sqrt(i^2+j^2)<L1 && sqrt(i^2+j^2)>0
             neighborhood_L1_nobounds(counter,1)=i;
            neighborhood_L1_nobounds(counter,2)=j;
           
            counter=counter+1;
        
        end
        
    end
end
neighborhood_1=[0,0;1,0;-1,0;0,1;0,-1];

%_______________________________________________________________
%where it can do work, antibiotics 
k=1;
where_a_moves=[0,0;
    k,-k;
    k,0;
    k,k;
    0,k;
    -k,k;
    -k,0;
    -k,-k;
    0,-k];

can_move=nan(y_pts,x_pts,9);
cell_counter=zeros(29,1);
cell_counter(1)=sum(sum(bacteria));
marker1=2;

 %where it can do work, antibiotics 
 for ii=20:y_pts-20
     for jj=20:x_pts-20
        %can_move(ii,jj,1)=1;
        marker=1;
    
        
        
        
        if ii+k<=y_pts-20 
            
                can_move(ii,jj,marker)=3;
                marker=marker+1;
             
            
            
        end 
        
        
        if jj+k<=x_pts-20
            
                can_move(ii,jj,marker)=5;
                marker=marker+1;
            
            
        end 
        
        if ii-k>=20
            
                can_move(ii,jj,marker)=7;
                marker=marker+1;
            
            
        end 
        
        
        
        if jj-k>=20
            
                can_move(ii,jj,marker)=9;
                marker=marker+1;
            
           
        end 
     end 
 end 



%____________________________________________________________



for k=1:nt
    change_bacteria=zeros(size(bacteria));
    change_antibiotics=zeros(size(antibiotics));
    change_sticky=zeros(size(sticky));
    for i=20:y_pts-20
        for j=20:x_pts-20

            if mod(k,1000)==0
    %____________________________________________________________
    %bacteria growth
    
                if bacteria(i,j)==0
                     space=1;
                    if sticky(i,j)==0
                        
                        for ii=1:length(neighborhood_L1_nobounds)
                            neighborhood_y=i+neighborhood_L1_nobounds(ii,1);
                            neighborhood_x=j+neighborhood_L1_nobounds(ii,2);
                            
                            if bacteria(neighborhood_y,neighborhood_x)==1
                                space=0;
                            end 
                            
                        end 
                  
                    end
                      %up and down is first coord
                    prob_0to0=1-sticky(i+1,j)*growth_prob(antibiotics(i+1,j))*(bacteria(i+1,j)/(2-bacteria(i+1,j)));
                    prob_0to0=prob_0to0*(1-sticky(i-1,j)*growth_prob(antibiotics(i-1,j))*(bacteria(i-1,j)/(2-bacteria(i-2,j))));
                    prob_0to0=prob_0to0*(1-space*(1-sticky(i+L1,j))*growth_prob(antibiotics(i+L1,j))*(bacteria(i+L1,j)/(2-bacteria(i+2*L1,j))));
                    prob_0to0=prob_0to0*(1-space*(1-sticky(i-L1,j))*growth_prob(antibiotics(i-L1,j))*(bacteria(i-L1,j)/(2-bacteria(i-2*L1,j))));
        
                    %left and right is second coord
                    prob_0to0=prob_0to0*(1-sticky(i,j+1)*growth_prob(antibiotics(i,j+1))*(bacteria(i,j+1)/(2-bacteria(i,j+2)))*bacteria(i-1,j+1)*bacteria(i+1,j+1));
                    prob_0to0=prob_0to0*(1-sticky(i,j-1)*growth_prob(antibiotics(i,j-1))*(bacteria(i,j-1)/(2-bacteria(i,j-2)))*bacteria(i-1,j-1)*bacteria(i+1,j-1));
                    prob_0to0=prob_0to0*(1-space*(1-sticky(i,j+L1))*growth_prob(antibiotics(i,j+L1))*(bacteria(i,j+L1)/(2-bacteria(i,j+2*L1)))*bacteria(i-L1,j+L1)*bacteria(i+L1,j+L1));
                    prob_0to0=prob_0to0*(1-space*(1-sticky(i,j-L1))*growth_prob(antibiotics(i,j-L1))*(bacteria(i,j-L1)/(2-bacteria(i,j-2*L1)))*bacteria(i-L1,j-L1)*bacteria(i+L1,j-L1));
                    
                    prob_0to1=1-prob_0to0;
                    if rand(1)<= prob_0to1
                        change_bacteria(i,j)=1; %goes to 1
                    end 
                else 
                     prob_1to0=death_prob;
           
                     prob_1to1=1-prob_1to0;
                     if rand(1)<= prob_1to0 
                       
                        change_bacteria(i,j)=-1; %goes to 0
                        if sticky(i,j)==1
                            change_sticky(i,j)=-1;
                        end 
                        
                        change_antibiotics(i,j)=-antibiotics(i,j);
               
                    end
    
                
    
                end 
                 
                
                
               
               
      %_________________________________________________________________
     
                %sticky switch
    
                if sticky(i,j)==0 && bacteria(i,j)==1 && change_bacteria(i,j)==0
                    
                   
    
                     %becoming sticky by antibiotics
                    if antibiotics(i,j)>=T_a 
                        change_sticky(i,j)=1;    
                    end 
    
                     %probabilities of becoming sticky via aggregation
                    sticky_prob_0to1=bacteria(i,j)*(1-(1-agg_prob*sticky(i+1,j))*(1-agg_prob*sticky(i-1,j))*(1-agg_prob*sticky(i,j+1))*(1-agg_prob*sticky(i,j-1)));
                    
    
    
                    %if become sticky by signaling OR aggregation
                    % else become sticky by aggregation alone
                    if sum(sum((antibiotics(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))>=T_a).*(bacteria(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))==1)))>=1

                        [y_distance_points,x_distance_points]=find(sticky(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))==1 & bacteria(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))==1);
                            
                         if ~isempty(y_distance_points)
                             
                             %to compute distance between our cell and sticky cells
                            points=cat(2,y_distance_points,x_distance_points);
                            current_y=round((min(i+stress_distance,y_pts)-max(i-stress_distance,1)+1)/2);
                            current_x=round((min(j+stress_distance,x_pts)-max(j-stress_distance,1)+1)/2);
                            
                            current_point=cat(2,current_x*ones(size(points,1),1),current_y*ones(size(points,1),1));
                            distances=sqrt((current_point(:,1)-points(:,1)).^2 +(current_point(:,2)-points(:,2)).^2);
                            min_distance=min(distances);
                          
                   
                        
                            %alpha_1 measures how close the cells are that are sticky
                             alpha_1=(-1/16)*min_distance+1; 
                            %alpha_2 is the density of sticky cells in the area
                            alpha_2=sum((sticky(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))),'all')/numel(ones(length(max(i-stress_distance,1):min(i+stress_distance,y_pts)),length(max(j-stress_distance,1):min(j+stress_distance,x_pts)))); 
                            prob_signal_0to1=.5*signal_stress_prob_max*alpha_1+.5*signal_stress_prob_max* alpha_2;
                            
    
                            if rand(1)<=prob_signal_0to1 || rand(1)<= sticky_prob_0to1
                                change_sticky(i,j)=1;
                            end 
                         end 
                       
                     
                        
                        
                    else 
                        if rand(1)<=sticky_prob_0to1
                            change_sticky(i,j)=1;
                        end
    
                    end 
                    
                end
            end
 %____________________________________________________________________
        %antibiotics move
   
            if (antibiotics(i,j)>.1 && bacteria(i,j)==0 && change_bacteria(i,j)~=-1) || (bacteria(i,j)==1 && antibiotics(i,j)>= T_g+T_d && change_bacteria(i,j)~=-1)
                nans=~isnan(can_move(i,j,:));
                nans=sum(nans);
        
               spot_a=can_move(i,j,1:nans);

               if bacteria(i,j)==1
                   number_of_jumps=(antibiotics(i,j)-(T_g+T_d))/nans; %sink that holds as much as V+epsilon
               else
                   number_of_jumps=antibiotics(i,j)/nans; %not sink
               end
               
               

                %choose where to do work
               
                if number_of_jumps>0
                   for p=1:length(spot_a)
                        
                        
                        move_y=where_a_moves(spot_a(p),1);
                        move_x=where_a_moves(spot_a(p),2);
                   
                        change_antibiotics(i+move_y,j+move_x)=number_of_jumps+change_antibiotics(i+move_y,j+move_x);
                        change_antibiotics(i,j)=-number_of_jumps+change_antibiotics(i,j);
                     
                    end 
                end
          
               
                
            end
            

            
            

        end 
    end 
    bacteria=bacteria+change_bacteria;
    antibiotics=antibiotics+change_antibiotics;
    sticky=sticky+change_sticky;
   
       antibiotics(:,20)=antibiotics(:,20)+ones(size(antibiotics(:,20)))*1.3; %set boundary conditions 1/5 100XMIC, 1 hour flux
       antibiotics(:,x_pts-20)=antibiotics(:,x_pts-20)+ones(size(antibiotics(:,20)))*1.3;
    
  
    t=t+dt;
    
 

end 

%Analysis of output

sizes=size(bacteria);
break_param=sizes(1)-20-break_param; %use the same break parameter as above, cutting off the rigid base

connected = bwareaopen(bacteria,threshol_vol,4); %gets rid of strucutres below volume threshold

break_param 

bacteria_cut=connected(1:break_param,:);
pts=[];
marker=1;
cloud_covered_pts=[];
T_d=15;
minpts=50;
sizes=size(bacteria_cut);    
 
for i=1:sizes(1)
    for j=1:sizes(2)
        if bacteria_cut(i,j)>0
            pts(marker,2)=sizes(1)-i;
            pts(marker,1)=j;
            marker=marker+1;
        end
    end
end
        
idx = dbscan(pts,T_d,minpts);
meanx=[];
meany=[];
uni=unique(idx);
        
sums_specific=[];
for num=1:length(uni)
    if uni(num)~=-1
        M=(idx==uni(num));
        B=double(M);
        new_pts=B.*pts;
        new_pts( ~any(new_pts,2), : ) = [];
        sums_specific=[sums_specific, sum(sum(B))];
       
            meanx(uni(num)) = mean(new_pts(:,2));
            meany(uni(num)) = mean(new_pts(:,1));
  
    end 

end 
%distances between structures
dist=pdist([meanx;meany]');
sums_specific   


figure(3)
hold on
gscatter(pts(:,2),pts(:,1),idx);
plot(meanx,meany,'kx',...
'MarkerSize',15,'LineWidth',3) 
hold off

%probability of a bacterial cell doubling
function g=growth_prob(a)
    T_g=1200;
    g_max=4/3*.028; %maximum chance of bacteria doubling
    g= g_max*(1-a/T_g).*(a<=T_g); %chance of bacteria doubling
end 
   






                
                  
           
                    
         
    
    

