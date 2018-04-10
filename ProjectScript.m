%% Cell Counting Script
% BME 6440: Final Project
%
% TeamIntes : Evan, Marien, & Denzel
% 
%% Reading Intensity Measurements
% Add Files & Functions to path
path1 = 'G:\Denzel\RPI\Year 2\Spring 2018\Classes - Spring 2018\BMED 6460 - Biological Image Analysis\Project';
cd(path1);
addpath(genpath(path1));
%
%
tic
% Extracting OCT Information
% 
% Open OCT file with 7zip. 
% Note: Name will need to be changed every time analysis is conducted
%
addpath(genpath('C:\Users\OCT\Desktop\OCT Scripts\'));

octFile = OCTFileOpen('D0-Sample1.oct');
%
% Extract intensity from OCT file
intensity = OCTFileGetIntensity(octFile); % Output is in form [z,x,y]