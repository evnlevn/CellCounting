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
% Extract intensity from OCT file
intensity = load('D1S1_intensity.mat'); % Output is in form [z,x,y]
intensity = permute(intensity,[3,1,2])
