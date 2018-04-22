close all
clc, clear
%% Cell Counting Script

% TeamIntes : Evan, Marien, & Denzel

%% Reading Intensity Measurements
% Extract intensity from OCT file
% intensity = load('D0S1_intensity.mat'); % Output is in form [z,x,y]
% intensity = intensity.intensity;

%% Extract data from .oct to .dat
% for y=1:size(intensity,3)
%     im=intensity(:,:,y);
%     fname=sprintf('%s_%d','day0',y);
%     fname=sprintf('%s%s',fname,'.dat');
%     csvwrite(fname,im);
% end

% %% Reading file and preprocessing
file=dir('C:\Users\evan_\Desktop\CellCounts\d1_test\*.dat');
data=[];
for i=1:length(file)
    fname=strcat('C:\Users\evan_\Desktop\CellCounts\d1_test\',file(i).name);
    data(:,:,i)=csvread(fname);
end

% Find target region based on middle image
im=data(:,:,round(size(data,3)/2));
figure, imshow(uint8(im))
[c,r]=ginput(2);
close all 

% Resize all images to target region and apply equalization, sharpening
rsz=[];
h=ones(5)/25;
k=1.3;
for i=1:size(data,3)
    im=data(:,:,i);
    im=im(r(1):r(2),:);
    imb=filter2(h,im);
    mask=double(im)-imb;
    im=double(im)+k*mask;
    im=im+abs(min(im(:)));
    im=floor(255.*(im./max(im(:))));
    rsz(:,:,i)=imgaussfilt(im,1);
%     figure, imshow(uint8(rsz(:,:,i)))
end
rsz=uint8(rsz);

%% Display preprocessed images and select training data windows
% pts=[];
% for i=1:size(rsz,3)
%     figure,imshow(rsz(:,:,i))
%     [c,r]=ginput(1);
%     pts=[pts; r c];
%     close all
% end
% 
% for i=1:size(pts,1)
%     imn=rsz(pts(i,1):pts(i,1)+7,pts(i,2):pts(i,2)+7);
% %     figure,imshow(imn)
%     fname=sprintf('%s_%d','dish',i);
%     fname=sprintf('%s%s',fname,'.dat');
%     csvwrite(fname,imn);
% end


%% Segment Dish
m=mean(rsz(:));
stdev=std(double(rsz(:)));
dishfilt=zeros(size(rsz));
sto=strel('disk',3);
stc=strel('disk',4);
% Set all dish values to be binary 1 
for i=1:size(rsz,3)
%     figure, histogram(rsz(:,:,i))
%     figure, imshow(rsz(:,:,i))
%     pause(0.5)
    dishfilt(:,:,i)=rsz(:,:,i)>uint8(m+2.5*stdev);
    dishfilt(:,:,i)=imdilate(dishfilt(:,:,i),sto);
    dishfilt(:,:,i)=imerode(dishfilt(:,:,i),stc);
%     figure,imshow(dishfilt(:,:,i))
end


%% K-means
% 
% for i=1:size(rsz,3)
%     % Generate seed points
%     seedbg=[20 20]; % Top left corner
%     seedcl=[round(size(rsz,1)/2) round(size(rsz,2)/2)];
%     
%     % Check that cell seed is not a dish element, changes if needed
%     notdish=0;
%     if dishfilt(seedcl(1),seedcl(2))==0
%         notdish=1; 
%     end
%     while notdish==0
%         seedcl=[seedcl(1)-10,seedcl(2)]
%         if dishfilt(seedcl(1),seedcl(2),i)==0
%             notdish=1;
%         end
%     end
%     
%     % Check that cell seed is in a cell region, changes if needed
%     cpos=0;
%     counter=0;
%     cont=[];
%     while cpos==0
%         % x = [ contrast correlation energy homogeneity ]
%         % contrast > 0.1 --> cell
%         % else --> background
%         x=texture_test(rsz(:,:,i),seedcl(1),seedcl(2));
%         if x(1)>0.1
%             cpos=1;
%         elseif counter<=5
%             counter=counter+1;
%             cont=[cont; seedcl(1) seedcl(2) x(1)];
%             seedcl=[seedcl(1)+8 seedcl(1)+8];
%         elseif counter>5
%             bc=[0 0];
%             for i=1:size(cont,1)
%                 if cont(i,3)>bc(1,2)
%                     bc=[i cont(i,3)];
%                 end
%             end
%             seedcl=[cont(bc(1,1),1) cont(bc(1,2),2)]
%             cpos=1;
%         end
%     end
% %     figure, imshow(rsz(:,:,i))
% %     hold on
% %     plot(seedcl(2),seedcl(1),'r*','MarkerSize',10)
%     cellcoords=[];
%     bgcoords=[];
%     mcell=texture_test(rsz(:,:,i),seedcl(1),seedcl(2));
%     mcell=mcell(1);
%     mbg=texture_test(rsz(:,:,i),seedbg(1),seedbg(2));
%     mbg=mbg(1);
%     for r=1:5:size(rsz,1)
%         for c=1:5:size(rsz,2)
%             if dishfilt(r,c)==0
%                 if r+7>size(rsz,1) | r-7<=0 | c+7>size(rsz,2) | c-7<=0
%                     bgcoords=[bgcoords; r c];
%                 else
%                     ctrt=texture_test(rsz(:,:,i),r,c);
%                     ctrt=ctrt(1);
%                     if abs(abs(ctrt(1))-abs(mcell))<abs(abs(ctrt(1))-abs(mbg))
%                         cellcoords=[cellcoords; r c];
% %                     elseif abs(abs(ctrt(1))-abs(mcell))>abs(abs(ctrt(1))-abs(mbg))
% %                         bgcoords=[bgcoords; r c];
%                     end
%                 end
%             end
%         end
%     end
%     figure,imshow(rsz(:,:,i))
%     hold on
%     for z=1:size(cellcoords,1)
%         plot(cellcoords(z,2),cellcoords(z,1),'r*','MarkerSize',1)
%         hold on
%     end
% %     for z=1:size(bgcoords,1)
% %         plot(bgcoords(z,1),bgcoords(z,2),'b*','MarkerSize',10)
% %         hold on
% %     end
% end



% 
rsz = double(rsz); % Converting Resized Image Datatype to double
[thickness,x,y] = size(rsz);
% Initialize Matrix
cluster = zeros(thickness,x,y);
    
% Engaging the K-means function on each slice
for i = 1:size(rsz,3)
    aa = squeeze(rsz(:,:,i)); % Reduce dimensions to 2
    aa = aa(:); % K-means fucntion takes vector inputs
    [idx,C] = kmeans(aa,3); % Run K-means with 3 clusters: cells, well, bckgrnd
    idx = reshape(idx,[thickness,512]); % Reshaping data into original matrix format
    % This nested for-loop enables the grouping of multiple clusters
    %         for j = 1:thickness
    %             for k = 1:512
    %                 if idx(j,k) == 2
    %                     idx(j,k) = 1; % Turning the well into the background
    %                 end
    %             end
    %         end
    cluster(:,:,i) = idx; % Storing the clusters results into a 3D output
end

figure('Name','K-Means Clustering Results');
count = 1;
increment = 1;
for i = increment:increment:size(rsz,3)
    figure
    imagesc(cluster(:,:,i)); axis image
    title(sprintf('Image of Slice #%d',i));
    count = count+1;
end
    
















