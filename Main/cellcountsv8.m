%% Cell Counting Script
% BME 6440: Final Project
%
% TeamIntes : Evan, Marien, & Denzel
%
close all;clc
%% Modulating Code
% All values of 1 run the corresponding sections
Modulate.Day = 2;
Modulate.ReadMat = 1;
Modulate.Conversion = 1;
Modulate.ReadDat = 1;
Modulate.Preprocessing = 1;
Modulate.K_Means = 1;
Modulate.Save = 1;
%% Reading Intensity Measurements from .mat files
% Extract intensity from OCT file
if Modulate.ReadMat == 1
    
    if Modulate.Day == 1
        intensity = load('Hi_res_S1D1.mat'); % Output is in form [z,x,y]
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
    
    if Modulate.Day == 1
        file=dir('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D1_Data\*.dat');
        for i=1:length(file)
            fname=strcat('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D1_Data\',file(i).name);
            data(:,:,i)=csvread(fname);
        end
    elseif Modulate.Day == 2
        file=dir('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D2_Data\*.dat');
        for i=1:length(file)
            fname=strcat('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D2_Data\',file(i).name);
            data(:,:,i)=csvread(fname);
        end
    elseif Modulate.Day == 3
        file=dir('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D3_Data\*.dat');
        for i=1:length(file)
            fname=strcat('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D3_Data\',file(i).name);
            data(:,:,i)=csvread(fname);
        end
    elseif Modulate.Day == 4
        file=dir('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D4_Data\*.dat');
        for i=1:length(file)
            fname=strcat('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D4_Data\',file(i).name);
            data(:,:,i)=csvread(fname);
        end
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
    h=ones(5)/25; % Sharpening Kernel
    k=0.8;
    [z,x,y] = size(data);
    depth = abs(round(r(2))-round(r(1))) + 1;
    rsz = zeros(depth,x,y);
    hboost = zeros(depth,x,y);
    for i=1:size(data,3)
        im=data(:,:,i);
        im=im(round(r(1)):round(r(2)),:);
        imb=filter2(h,im);
        mask=double(im)-imb;
        im=double(im)+k*mask;
        im=im+abs(min(im(:)));
        im=floor(255.*(im./max(im(:))));
        hboost(:,:,i) = im;
        rsz(:,:,i)=imgaussfilt(im,3);
    end
    
    % Middle Slice of Sharpened Image
    [z,x,y] = size(hboost);
    figure('Name','Image Check of Middle Slice of Sharpened Image')
    imagesc(hboost(:,:,round(y/2))); axis image
    
    % Segment Dish
    % Generate Volume Statistics
    m = mean(rsz(:));
    stdev=std(double(rsz(:)));
    dishfilt=zeros(size(rsz));
    sto=strel('disk',4); % Opening/dilation
    stc=strel('disk',3); % Erode/erosion
    % Set all dish values to be binary 1
    %     figure
    for i=1:size(hboost,3)
        dishfilt(:,:,i)=rsz(:,:,i)>uint8(m+2*stdev);
        dishfilt(:,:,i)=imdilate(dishfilt(:,:,i),sto);
        dishfilt(:,:,i)=imerode(dishfilt(:,:,i),stc);
        %         imshow(dishfilt(:,:,i))
        %         i
    end
    
    %Subtract Dish from Images
    m = mean(hboost(:));
    stdev=std(double(hboost(:)));
    [z,x,y] = size(dishfilt);
    for i = 1:z
        for j = 1:x
            for k = 1:y
                if dishfilt(i,j,k) == 1
                    hboost(i,j,k) = m - 0.5*stdev;
                elseif hboost(i,j,k) <=  m + 0.75*stdev
                    hboost(i,j,k) = m - 0.5*stdev;
                end
            end
        end
    end
    
    
    stc=strel('square',2); % Erode/erosion
    sto=strel('square',1); % Opening/dilation
    for i=1:size(hboost,3)
        hboost(:,:,i)=imerode(hboost(:,:,i),stc);
        hboost(:,:,i)=imdilate(hboost(:,:,i),sto);
    end
    
    % Display Preprocessed images
    hboost=uint8(hboost);
    for i=1:round(size(hboost,3)/8):size(hboost,3)
        %         figure, imshow(hboost(:,:,i))
        %         figure, histogram(hboost(:,:,i))
        figure('Name','Image and Histogram of Slices');
        suptitle(sprintf('Slice #%d : Image & Histogram ',i));
        subplot(2,1,1)
        imagesc(squeeze(hboost(:,:,i))); colormap jet; axis image;
        title(sprintf('Image of Slice #%d',i));
        subplot(2,1,2)
        histogram(hboost(:,:,i));  title('Histogram of Pixel Intensities');
    end
    %
end
%% K-means Clustering
%
if Modulate.K_Means == 1
    
    hboost = double(hboost); % Converting Resized Image Datatype to double
    [thickness,x,y] = size(hboost);
    % Initialize Matrix
    cluster = zeros(thickness,x,y);
    
    % Engaging the K-means function on each slice
    for i = 1:size(hboost,2)
        aa = squeeze(hboost(:,:,i)); % Reduce dimensions to 2
        aa = aa(:); % K-means fucntion takes vector inputs
        [idx,C] = kmeans(aa,2,'Start',[(29*z + 30);round((z*x/2)-1 + y/2)]); % Run K-means with 3 clusters: cells, well, bckgrnd
        idx = reshape(idx,[thickness,size(hboost,2)]); % Reshaping data into original matrix format
        
        cluster(:,:,i) = idx; % Storing the clusters results into a 3D output
    end
    
    figure('Name','K-Means Clustering Results');
    count = 1;
    increment = y/8;
    for i = increment:increment:size(hboost,2)
        subplot(3,3,count)
        imagesc(cluster(:,:,i)); axis image;
        title(sprintf('Image of Slice #%d',i));
        count = count+1;
    end
    
    % Semi-Binary Final Image Output
    cluster = floor(cluster./2);
    win = zeros(z,x,y);
    for i = 1:z
        for j = 1:x
            for k = 1:y
                if cluster(i,j,k) == 1
                    cluster(i,j,k) = 0;
                elseif cluster(i,j,k) == 0
                    cluster(i,j,k) = 1;
                end
                win(i,j,k) = hboost(i,j,k).*cluster(i,j,k);
            end
        end
    end

end

%% Save Image
if Modulate.Save == 1
    
    if Modulate.Day == 1
         cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D1_Data');
    elseif Modulate.Day == 2
         cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D2_Data');
    elseif Modulate.Day == 3
         cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D3_Data');
    elseif Modulate.Day == 4
         cd('C:\Users\OCT\Desktop\BMED6640\Sample Matrices\D4_Data');
    end
    
    
    win = uint8(win);
    
    for i = 1:size(cluster,2)
        bb = squeeze(win(:,:,i));
        name = sprintf('%s_%d%s','Day',i,'.png');
        imwrite(bb,name);
    end
   
end