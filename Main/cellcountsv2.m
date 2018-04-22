%% Cell Counting Script
% BME 6440: Final Project
%
% TeamIntes : Evan, Marien, & Denzel
%
close all;clear;clc
%% Modulating Code
% All values of 1 run the corresponding sections
Modulate.ReadMat = 0;
Modulate.Conversion = 0;
Modulate.ReadDat = 1;
Modulate.Preprocessing = 1;
Modulate.K_Means = 1;
%% Reading Intensity Measurements from .mat files
% Extract intensity from OCT file
if Modulate.ReadMat == 1
    
    intensity = load('D1S1_intensity.mat'); % Output is in form [z,x,y]
    intensity = intensity.intensity;
    
end
%
%% Converting Data from .mat to .dat
if Modulate.Conversion == 1
    
    FolderToSave = 'G:\Denzel\RPI\Year 2\Spring 2018\Classes - Spring 2018\BMED 6460 - Biological Image Analysis\Project';
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
    
    file=dir('G:\Denzel\RPI\Year 2\Spring 2018\Classes - Spring 2018\BMED 6460 - Biological Image Analysis\Project\20180417T175923\*.dat');
    data=zeros(1024,512,512);
    for i=1:length(file)
        fname=strcat('G:\Denzel\RPI\Year 2\Spring 2018\Classes - Spring 2018\BMED 6460 - Biological Image Analysis\Project\20180417T175923\',file(i).name);
        data(:,:,i)=csvread(fname);
    end
    
end
%% Normalize Raw Data & Apply Unsharp Masking (Preprocessing)
%
if Modulate.Preprocessing == 1
    % Find target region based on middle image
    im=data(:,:,round(size(data,3)/2));
    figure, imshow(uint8(im))
    % User input to determine ROI
    [c,r] = ginput(2)  % First click above sample, second click below
    close all
    
    % Resize all images to ROI and apply normalization and sharpening
    rsz=[]; % Initialzing 3D Data Cube of Resized Image
    h=ones(5)/25; % Sharpening Kernel
    k=1.3;
    for i=1:size(data,3)
        im=data(:,:,i);
        im=im(round(r(1)):round(r(2)),:);
        imb=filter2(h,im);
        mask=double(im)-imb;
        im=double(im)+k*mask;
        im=im+abs(min(im(:)));
        im=floor(255.*(im./max(im(:))));
        rsz(:,:,i)=im;
    end
    
    % Display Preprocessed images
    rsz=uint8(rsz);
    for i=1:64:size(rsz,3)
        %         figure, imshow(rsz(:,:,i))
        %         figure, histogram(rsz(:,:,i))
        figure('Name','Image and Histogram of Slices');
        suptitle(sprintf('Slice #%d : Image & Histogram ',i));
        subplot(2,1,1)
        imagesc(squeeze(rsz(:,:,i))); colormap jet; axis image;
        title(sprintf('Image of Slice #%d',i));
        subplot(2,1,2)
        histogram(rsz(:,:,i));  title('Histogram of Pixel Intensities');
    end
    %
    
end
%% K-means Clustering
%
if Modulate.K_Means == 1
    
    rsz = double(rsz); % Converting Resized Image Datatype to double
    [thickness,x,y] = size(rsz);
    % Initialize Matrix
    cluster = zeros(thickness,x,y);
    
    % Engaging the K-means function on each slice
    for i = 1:512
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
    increment = 32;
    for i = increment:increment:512
        subplot(8,2,count)
        imagesc(cluster(:,:,i)); axis image
        title(sprintf('Image of Slice #%d',i));
        count = count+1;
    end
    
end
