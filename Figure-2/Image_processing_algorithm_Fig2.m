%Image processing algorithm 
% slice finding algorithm for sufficient initial conditions 
%reslice tiff files starting from 'left' x-slices

dimensions=256;
% Modify the following grey-scale threshold so that the base of the image
% is approximately the observed percentage of density as computed in the 
% density plot, normally treated biofilms will need a lower threshold
% because they are taller,larger structures and have problems absorbing
% signal
threshold=200; %above what grey scale do we want to use
connectivity=26; %connectivity parameter: can be 6,18,26
threshol_vol=100; %looking for strucutres above a certain number of voxels

%read in tiff files (each slice needs to be a seperate file)
files = dir('*.tif');
pts=[];
%creates binary matrix taking out grey values below threshold
A = tiffreadVolume(files(1).name);
u=(A>=threshold);

%condenses to 1um^3 squares
averaging = @(block_struct) mean(block_struct.data);
blockMeans=permute(u,[1 3 2]);
blockMeans = blockproc(blockMeans, [1 4], averaging);
compressed=permute(blockMeans,[1,3,2]); % als run below code if using tiff of 1024x1024
%{
blockMeans = blockproc(u, [1 4], averaging);
blockMeans = blockproc(blockMeans, [4 1], averaging);
%}

%removes some noise and re-binarizes
compressed=(compressed>.25); 
compressed=double(compressed);
density_threshold=.25; %cellular density that defines no longer rigid base
sizes=size(compressed);
break_param=0;
density=[];



%remove rigid layers with density
for k=1:sizes(3)
    density(k)=sum(sum(compressed(:,:,k)))/(sizes(1)*sizes(2));
    if break_param==0
        if sum(sum(compressed(:,:,k)))/(sizes(1)*sizes(2))<= density_threshold
            break_param=k; 
        end 
    end 
end 
figure(1)
hold on
plot(density,'*')
title('density of z-slices T0.0')
hold off
no_rigid_base=compressed(:,:,break_param:end);
BW2 = bwareaopen(no_rigid_base,100,connectivity); %gets rid of strucutres below volume threshold
L=bwlabeln(BW2,connectivity); %labels/finds connected strucutures
V=categorical(L);
no_rigid_base=(L>0);


%searches for y-slices where distance and number thresholds are met
sizes=size(no_rigid_base);
good_x=zeros(sizes(2),1);
vols_x=zeros(sizes(2),1);
dist_threshold=30; %sufficient protected region
number_threshold=50; %sufficient number of bacteria cells per structure
b=permute(no_rigid_base,[1 3 2]);
epsilon=15; %radius of search for points to include within structure
minpts=number_threshold;
dist_edge_threshold=20; %sufficient distance from boundary
sums_avg=[]; %average structure size of good slices by dbscan

%going through slices to conduct dbscan
for j=1:sizes(2)
    pts=[];
    marker=1;
   

    
    vols_x(j)=sum(sum(b(:,:,j)));

   
    if sum(sum(b(:,:,j)))>=number_threshold
        for i=1:sizes(1)
            for k=1:sizes(3)
                if b(i,k,j)>0
                    pts(marker,2)=sizes(1)-i;
                    pts(marker,1)=k;
                    marker=marker+1;
                end
            end
        end
        
        idx = dbscan(pts,epsilon,minpts);
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
                if mean(new_pts(:,2))>= dist_edge_threshold && (sizes(1)-mean(new_pts(:,2)))>=dist_edge_threshold
                    meanx(uni(num)) = mean(new_pts(:,2));
                    meany(uni(num)) = mean(new_pts(:,1));
                end 
            end 

        end 
        %distances between structures
        dist=pdist([meanx;meany]');
       

        if isempty(dist)==0
            %determine if slice meets criteria 
            if max(dist)>= dist_threshold 
                
                    good_x(j)=1;
                    sums_avg=[sums_avg,mean(sums_specific)];
                    
                
            end 
        end
    end 
    
    
end 

valid_slices_x=find(good_x==1);
%j values where good==1 (multiplied by 4) gives x value of which slice is a sufficient initial condition
percentagegood_x=sum(good_x)/256
mean(sums_avg)

figure(2)
hold on
title('good x-slices T0.0')
imagesc(good_x)
hold off

figure(3)
hold on
title('volume per slice')
colorbar()
imagesc(vols_x)
hold off

%{
figure(1)
hold on
gscatter(pts(:,2),pts(:,1),idx);
plot(meanx,meany,'kx',...
     'MarkerSize',15,'LineWidth',3) 
hold off
%}

