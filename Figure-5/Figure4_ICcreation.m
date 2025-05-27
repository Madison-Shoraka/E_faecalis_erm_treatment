%Run this file first to read in tiff file and create IC, only 2D tiff (y,z)
%microscopy tiff file allowed

%read in tiff files (each slice needs to be a seperate file)
files = dir('*.tif');
pts=[];
%creates binary matrix taking out grey values below threshold
k=1; %choose which tiff file to run
t = Tiff(files(k).name,'r');
A = read(t);

sizes=size(A);
% structure finding algorithm 
len=1024;
domain_length=len/4;
height=sizes(1); %set height of tiff file here!
L1=5;
T=1;
dt=.000007; %everything particle has a 1 chance of moving
nt=T/dt;

threshold=240; %above what grey scale do we want to use
num_tiffs=1023; %number of tiff images in the stack
u=zeros(len,height);
size(u)
connectivity=26; %connectivity parameter: can be 6,18,26
threshol_vol=10; %looking for strucutres above a certain number of voxels


    

for i=1:height
    for j=1:len
        if A(i,j)<threshold
            u(i,j)=0;
        end
        if A(i,j)>=threshold
            u(i,j)=1;
        end
    end 
end
u=u(1:height,:);
  
myMeanFunction = @(block_struct) mean(block_struct.data);
blockMeans = blockproc(u, [1 4], myMeanFunction);
blockMeans = blockproc(blockMeans, [4 1], myMeanFunction);
u=blockMeans;
sizes=size(u);



bacteria=zeros(sizes(1)+60,sizes(2)+40);
bacteria(40:sizes(1)+60-21,20:sizes(2)+19)=rot90(u,2);

sizes=size(bacteria);
antibiotics=zeros(size(bacteria));
stress_mat=zeros(size(bacteria));


imagesc(bacteria)


domain_x=[0,1];
domain_y=[0,1];

x_pts=sizes(2);
y_pts=sizes(1);
antibiotics(:,x_pts-20)=ones(size(antibiotics(:,1)))*.7392; 
x = linspace(domain_x(1),domain_x(2),x_pts); % x-values
y = linspace(domain_y(1),domain_y(2),y_pts); % y-values
[X,Y] = meshgrid(x,y); % 2d grid


V=1200; %number of antibiotics causing cell to stop growing
epsilon=10000;
a_s=66; %antibiotics necessary for bacteria to become sticky
g_max=.028; %maximum chance of bacteria doubling
%growth_prob=@(a) g_max*(1-a/V).*(a<=V); %chance of bacteria doubling
agg_prob=0.02; %probability a nearby cell will become sticky if next to a sticky cell
death_prob=g_max/3; %probability of bacteria dying a natural death
signal_stress_prob_max=.25;
stress_distance=15;

bacteria=bacteria>.5;
bw2=bwareaopen(bacteria,threshol_vol);
L=bwlabeln(bw2,connectivity); %labels/finds connected strucutures

for i=30:sizes(1)-20
    for j=20:sizes(2)-20
        if mod(i,L1)==0 && mod(j,L1)==0
            if bacteria(i,j)==0
                bacteria(i,j)=1;
            end 
        end
        if sum(sum(bacteria(i-L1+1:i+L1-1,j-L1+1:j+L1-1)))==0
             bacteria(i,j)=1;
        end
    end 
end 

sticky=L>0;
imagesc(bacteria)
bacteria=double(bacteria);

