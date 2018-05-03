%% Cell Counting Script
% BME 6440: Final Project
%
% TeamIntes : Evan, Marien, & Denzel
%
close all;clc
%% Modulating Code
% All values of 1 run the corresponding sections
Modulate.Day = 0;
Modulate.ReadMat = 0;
Modulate.Conversion = 0;
Modulate.ReadDat = 1;
Modulate.Preprocessing = 1;
Modulate.K_Means = 1;
Modulate.Save = 0;
%% Reading Intensity Measurements from .mat files
% Extract intensity from OCT file
if Modulate.ReadMat == 1
    
    if Modulate.Day == 1
        intensity = load('Hi_res_S1D1.mat'); % Output is in form [z,x,y]
%         intensity = load('Hi_res_S1D1.mat'); % Output is in form [z,x,y]
    elseif Modulate.Day == 2
        intensity = load('Hi_res_S1D2.mat'); % Output is in form [z,x,y]
    elseif Modulate.Day == 3
        intensity = load('Hi_res_S1D3.mat'); % Output is in form [z,x,y]
    elseif Modulate.Day == 4
        intensity = load('Hi_res_S1D4.mat'); % Output is in form [z,x,y]
    end
    
    intensity = intensity.intensity;
    
end
%
%% Converting Data from .mat to .dat
if Modulate.Conversion == 1
    
    FolderToSave = cd;
    % Generating the folder to save
    cd(FolderToSave);
    datetime=datestr(now,30);
    mkdir(num2str(datetime))
    filename3 = num2str(datetime);
    cd(filename3)
    %
    % Saving as .dat
    for y=1:size(intensity,3)
        im=intensity(:,:,y);
        fname=sprintf('%s_%d','day1',y);
        fname=sprintf('%s%s',fname,'.dat');
        csvwrite(fname,im);
    end
    
end
%% Reading .dat File and Storing into a 3D Matrix
if Modulate.ReadDat == 1
    data=[];
    file=dir('C:\Users\evan_\Desktop\HighRes\lowres\*.dat');
    for i=1:length(file)
        name=sprintf('%s_%d','day1',i);
        name=sprintf('%s%s',name,'.dat');
        fname=strcat('C:\Users\evan_\Desktop\HighRes\lowres\',name);
        data(:,:,i)=csvread(fname);
    end
end

%% Normalize Raw Data & Apply Unsharp Masking (Preprocessing)
%
if Modulate.Preprocessing == 1
    % Find target region based on middle image
    im=data(:,:,round(size(data,3)/2));
    figure; imshow(uint8(im))
    % User input to determine ROI
    [c,r] = ginput(2);  % First click above sample, second click below
    close all

    % Resize all images to ROI and apply normalization and sharpening
    rsz=[]; % Initialzing 3D Data Cube of Resized Image
    h=ones(5)/25; 
    k=1;  
    [z,x,y] = size(data);
    depth = abs(round(r(2))-round(r(1))) + 1;
    rsz = zeros(depth,x,y);
    hboost = zeros(depth,x,y);
    for i=1:size(data,3)
        im=data(:,:,i);
        % Window selection
        im=im(round(r(1)):round(r(2)),:);
        % H-boost/Unsharp masking 
        imb=filter2(h,im);
        mask=double(im)-imb;
        im=double(im)+k*mask;
        im=im+abs(min(im(:)));
        % Contrast equalization
        im=floor(255.*(im./max(im(:))));
        hboost(:,:,i) = im;
        % Blurred version copy
        rsz(:,:,i)=imgaussfilt(im,3); 
    end

    % Segment Dish
    % Generate Volume Statistics
    dishfilt=zeros(size(rsz));
    sto=strel('disk',4); % Opening/dilation
    stc=strel('square',3); % Erode/erosion
    % Set all dish values to be binary 1
    for i=1:size(hboost,3)
        a=rsz(:,:,i);
        m = mean(a(:));
        stdev=std(double(a(:)));
        dishfilt(:,:,i)=rsz(:,:,i)>uint8(m+1.8*stdev); 
        % Expand initial values to be more representative of dish
        dishfilt(:,:,i)=imdilate(dishfilt(:,:,i),sto);
        dishfilt(:,:,i)=imerode(dishfilt(:,:,i),stc);
        % Display dish segmentation
%         figure
%         imshow(dishfilt(:,:,i))
    end
    
    % Set intensity for background and dish values
    [z,x,y] = size(dishfilt);
    for k = 1:y
        a=double(hboost(:,:,y));
%         figure, histogram(rsz(:,:,k))
        m=rsz(round(depth/2),round(size(hboost,2)/2),k);
        stdev=std(a(:));
        for i = 1:z
            for j = 1:x
                if dishfilt(i,j,k) == 1
                    hboost(i,j,k) = m-2.5*stdev;
                elseif rsz(i,j,k)<= m-0.3*stdev
                    hboost(i,j,k) = m-2.5*stdev; 
                end
            end
        end
    end
    
    % Display initial segmentation
    for i=1:size(hboost,3)
        figure
        imshow(uint8(hboost(:,:,i)))
        pause(0.1)
    end
        
%     Enable for high res to get rid of speckles
%     stc=strel('square',2); % Erode/erosion
%     sto=strel('square',2); % Opening/dilation
%     for i=1:size(hboost,3)
%         hboost(:,:,i)=imerode(hboost(:,:,i),stc);
%         hboost(:,:,i)=imdilate(hboost(:,:,i),sto);
%     end
%     figure
%     for i=1:size(hboost,3)
%         pair=[uint8(data(r(1):r(2),:,i));hboost(:,:,i)];
%         imshow(pair,[])
%         pause(0.5)
%     end
    
end
%% K-means Clustering
%
if Modulate.K_Means == 1
    
    hboost = double(hboost); % Converting Resized Image Datatype to double
    [thickness,x,y] = size(hboost);
    % Initialize Matrix
    cluster = zeros(thickness,x,y);
    
    % Engaging the K-means function on each slice
    for i = 1:size(hboost,3)
        aa = squeeze(hboost(:,:,i)); % Reduce dimensions to 2
        aa = aa(:); % K-means fucntion takes vector inputs
        m=rsz(round(depth/2),round(size(hboost,2)/2),k);
        cond=0; rb=30; cb=30;
        while cond==0
            if hboost(rb,cb)<=m
                cond=1;
            elseif hboost(rb,cb)>m & cb==size(hboost,2)
                rb=rb+1;
                cb=1;
            else
                cb=cb+1;
            end
        end
        [idx,C] = kmeans(aa,2,'Start',[(rb*z + cb);round((z*x/2)-1 + y/2)]); % Run K-means with 3 clusters: cells, well, bckgrnd
        idx = reshape(idx,[thickness,size(hboost,2)]); % Reshaping data into original matrix format
        
        cluster(:,:,i) = idx; % Storing the clusters results into a 3D output
    end
    

    
    % Semi-Binary Final Image Output
    cluster = floor(cluster./2);
    win = zeros(z,x,y);
    for i = 1:z
        for j = 1:x
            for k = 1:y
                if cluster(i,j,k) == 0  % Day 1 == 1
                    cluster(i,j,k) = 0;
                elseif cluster(i,j,k) == 1 % Day 1 == 0
                    cluster(i,j,k) = 1;
                end
                win(i,j,k) = hboost(i,j,k).*cluster(i,j,k);  
            end
        end
    end
    win=uint8(win);
%     stc=strel('square',2); % Erode/erosion
%     sto=strel('square',2); % Opening/dilation
%     for i=1:size(win,3)
% %         win(:,:,i)=imerode(win(:,:,i),stc);
%         win(:,:,i)=imdilate(win(:,:,i),sto);
%     end
    
%     Display K-means results
    for i=1:size(win,3)
        figure
        imshow(uint8(win(:,:,i)))
    end
    
    
end

if Modulate.Save == 1
    
    if Modulate.Day == 1
%          cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D1_Data');
    elseif Modulate.Day == 2
%          cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D2_Data');
    elseif Modulate.Day == 3
%          cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D3_Data');
    elseif Modulate.Day == 4
%          cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D4_Data');
    end
    
    
    win = uint8(win);
    
    for i = 1:size(cluster,3)
        bb = squeeze(win(:,:,i));
        name = sprintf('%s_%d%s','Day',i,'.tif');
        imwrite(bb,name);
    end
   
end