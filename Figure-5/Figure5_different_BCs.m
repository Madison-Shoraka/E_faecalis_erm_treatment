%Run this to simulate antibiotic treatment, Note change lines 270-278
%depending on what type of BC you want

T=1;
dt=.000007; %everything particle has a 1 chance of moving
nt=T/dt;
t=0;
neighborhood_L1=[];
neighborhood_L1_nobounds=[];
counter=1;

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
vid = VideoWriter("pretty_picture.mp4",'MPEG-4');
open(vid)
writeVideo(vid,bacteria)
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
                    if antibiotics(i,j)>=a_s 
                        change_sticky(i,j)=1;    
                    end 
    
                     %probabilities of becoming sticky via aggregation
                    sticky_prob_0to1=bacteria(i,j)*(1-(1-agg_prob*sticky(i+1,j))*(1-agg_prob*sticky(i-1,j))*(1-agg_prob*sticky(i,j+1))*(1-agg_prob*sticky(i,j-1)));
                    
    
    
                    %if become sticky by signaling OR aggregation
                    % else become sticky by aggregation alone
                    if sum(sum((antibiotics(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))>=a_s).*(bacteria(max(i-stress_distance,1):min(i+stress_distance,y_pts),max(j-stress_distance,1):min(j+stress_distance,x_pts))==1)))>=1

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
                            prob_signal_0to1=signal_stress_prob_max*alpha_1+signal_stress_prob_max*alpha_2;
                            
    
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
   
            if (antibiotics(i,j)>.1 && bacteria(i,j)==0 && change_bacteria(i,j)~=-1) || (bacteria(i,j)==1 && antibiotics(i,j)>= V+epsilon && change_bacteria(i,j)~=-1)
                nans=~isnan(can_move(i,j,:));
                nans=sum(nans);
        
               spot_a=can_move(i,j,1:nans);

               if bacteria(i,j)==1
                   number_of_jumps=(antibiotics(i,j)-(V+epsilon))/nans; %sink that holds as much as V+epsilon
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
    %(A) bulk BC's
    %antibiotics(20,:)=ones(size(antibiotics(20,:)))*6600; %set boundary conditions bulk
        %(B): Flux no reduction
       %antibiotics(:,20)=antibiotics(:,20)+ones(size(antibiotics(:,20)))*2.9568; %set boundary conditions 100XMIC, 2 hour flux
       %antibiotics(:,x_pts-20)=antibiotics(:,x_pts-20)+ones(size(antibiotics(:,20)))*2.9568;
      %(C): Flux with reduction
       antibiotics(:,20)=antibiotics(:,20)+ones(size(antibiotics(:,20)))*.989; %set boundary conditions 1/5 100XMIC, 1 hour flux
       antibiotics(:,x_pts-20)=antibiotics(:,x_pts-20)+ones(size(antibiotics(:,20)))*.989;
       
  
    t=t+dt;
    if mod(k,200)==0
        writeVideo(vid,bacteria)
        cell_counter(marker1)=sum(sum(bacteria));
        marker1=marker1+1;
    end 
 

end 
close(vid)

%computes growth probability from cell given intracellular anitbiotic
%concentration, a
function g=growth_prob(a)
    V=1200;
    g_max=4/3*.028; %maximum chance of bacteria doubling
    g= g_max*(1-a/V).*(a<=V); %chance of bacteria doubling
end 
   

