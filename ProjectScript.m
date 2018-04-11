%% Cell Counting Script
% BME 6440: Final Project
%
% TeamIntes : Evan, Marien, & Denzel
% 
close all 
clear all
clc
%% Reading Intensity Measurements
% Add Files & Functions to path
path1 = 'G:\Denzel\RPI\Year 2\Spring 2018\Classes - Spring 2018\BMED 6460 - Biological Image Analysis\Project';
cd(path1);
addpath(genpath(path1));
%
% Load intensity matrix and convert to Magnitude
intensity = load('D1S1_intensity.mat'); % Output is in form [z,x,y]
intensity = intensity.intensity;

figure('Name','Default Image: B-Scan 256')
imagesc(intensity(:,:,256)); axis image

%% Unsharp Masking
lo_pass = ones(5)./25;
inten_filt = imfilter(intensity,lo_pass);
mask = intensity-inten_filt;
k = 1.2; % Sharpening parameter
intensity = intensity - k*mask;
figure('Name','Sharpened Image: B-Scan 256')
intensity = 255*intensity./max(intensity(:));
imagesc(intensity(:,:,256)); axis image

% Crop Image
intensity = intensity(550:650,:,:);
figure;histogram(intensity)
% Thresholding
C = 4; D0 = intensity;
D0 = D0 > (mean(D0) - C) & D0 < (mean(D0) + C);
figure('Name','Thresholded Cropped Image'); 
imagesc(D0(:,:,256)); axis image
