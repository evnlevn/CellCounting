close all
clc, clear
%% Cell Counting Script

% TeamIntes : Evan, Marien, & Denzel

%% Reading Intensity Measurements
% Extract intensity from OCT file
% intensity = load('D1S1_intensity.mat'); % Output is in form [z,x,y]
% intensity = intensity.intensity;

%% Extract data from .oct to .dat
% for y=1:size(intensity,3)
%     im=intensity(:,:,y);
%     fname=sprintf('%s_%d','day1',y);
%     fname=sprintf('%s%s',fname,'.dat');
%     csvwrite(fname,im);
% end

%% Reading file and preprocessing
file=dir('C:\Users\evan_\Desktop\CellCounts\d1_test\*.dat');
data=[];
for i=1:length(file)
    fname=strcat('C:\Users\evan_\Desktop\CellCounts\d1_test\',file(i).name);
    data(:,:,i)=csvread(fname);
end

% Find target region based on middle image
im=data(:,:,round(size(data,3)/2));
figure, imshow(uint8(im))
[c,r]=ginput(2)
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
    rsz(:,:,i)=im;
end

% Display preprocessed images
rsz=uint8(rsz);
for i=1:size(rsz,3)
    figure, imshow(rsz(:,:,i))
%     figure, histogram(rsz(:,:,i))
end

for i=1:length(rsz)
    fname=sprintf('%s_%d','day1',y);
    fname=sprintf('%s%s',fname,'.png');
    imwrite(rsz(:,:,i),fname)
end



% 
% imn=imread('all.png');
% for r=1:size(imn,1)
%     for c=1:size(imn,2)
%         if imn(r,c)==0
%             imn(r,c)=im(r,c);
%         end
%     end
% end

% figure, imshow(imn)
% figure; histogram(imn)
% 
% mx=max(imn(:));
% mn=min(imn(:));
% imn = 255./(mx-mn).*(imn-mn);
% figure; imshow(imn)
% figure; histogram(imn)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
