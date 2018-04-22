clc
clear all
close all

offsets=[1 1; 1 2; 2 1; 2 2; 1 3; 3 1; 2 3; 3 2; 3 3];

%% Cells
file=dir('C:\Users\evan_\Desktop\Test\cells\*.dat');
contrast=[];correlation=[];energy=[];homogeneity=[];

for x=1:size(offsets,1)
    for i=1:length(file)
        if i==1
            ccon=[]; ccor=[]; cen=[]; cho=[];
        end
        fname=strcat('C:\Users\evan_\Desktop\Test\cells\',file(i).name);
        im=csvread(fname);
        glcm=graycomatrix(uint8(im),'Offset',[offsets(x,1) offsets(x,2)]);
        stats=graycoprops(glcm,'all');
        ccon=[ccon;stats.Contrast()];
        ccor=[ccor;stats.Correlation()];
        cen=[cen;stats.Energy()];    
        cho=[cho;stats.Homogeneity()];
        if i==length(file)
            contrast(:,:,x)=ccon; correlation(:,:,x)=ccor;
            energy(:,:,x)=cen; homogeneity(:,:,x)=cho;
        end
    end
end

cell_vals=[];
for i=1:9
%     avg=[contrast(:,:,i) correlation(:,:,i) energy(:,:,i) homogeneity(:,:,i)];
    avg=[mean(contrast(:,:,i)) mean(correlation(:,:,i)) mean(energy(:,:,i)) mean(homogeneity(:,:,i))];
    cell_vals=[cell_vals;avg];
end
figure, plot(cell_vals(:,1),cell_vals(:,3),'r*','MarkerSize',10)
hold on
cell_vals=[mean(cell_vals(:,1)) mean(cell_vals(:,2)) mean(cell_vals(:,3)) mean(cell_vals(:,4))]


%% Dish
file=dir('C:\Users\evan_\Desktop\Test\dish\*.dat');
contrast=[];correlation=[];energy=[];homogeneity=[];

for x=1:size(offsets,1)
    for i=1:length(file)
        if i==1
            ccon=[]; ccor=[]; cen=[]; cho=[];
        end
        fname=strcat('C:\Users\evan_\Desktop\Test\dish\',file(i).name);
        im=csvread(fname);
        glcm=graycomatrix(uint8(im),'Offset',[offsets(x,1) offsets(x,2)]);
        stats=graycoprops(glcm,'all');
        ccon=[ccon;stats.Contrast()];
        ccor=[ccor;stats.Correlation()];
        cen=[cen;stats.Energy()];    
        cho=[cho;stats.Homogeneity()];
        if i==length(file)
            contrast(:,:,x)=ccon; correlation(:,:,x)=ccor;
            energy(:,:,x)=cen; homogeneity(:,:,x)=cho;
        end
    end
end

dish_vals=[];
for i=1:9
%     avg=[contrast(:,:,i) correlation(:,:,i) energy(:,:,i) homogeneity(:,:,i)];
    avg=[mean(contrast(:,:,i)) mean(correlation(:,:,i)) mean(energy(:,:,i)) mean(homogeneity(:,:,i))];
    dish_vals=[dish_vals;avg];
end
plot(dish_vals(:,1),dish_vals(:,3),'g*','MarkerSize',10)
hold on
dish_vals=[mean(dish_vals(:,1)) mean(dish_vals(:,2)) mean(dish_vals(:,3)) mean(dish_vals(:,4))]

%% Background
file=dir('C:\Users\evan_\Desktop\Test\bg\*.dat');
contrast=[];correlation=[];energy=[];homogeneity=[];

for x=1:size(offsets,1)
    for i=1:length(file)
        if i==1
            ccon=[]; ccor=[]; cen=[]; cho=[];
        end
        fname=strcat('C:\Users\evan_\Desktop\Test\bg\',file(i).name);
        im=csvread(fname);
        glcm=graycomatrix(uint8(im),'Offset',[offsets(x,1) offsets(x,2)]);
        stats=graycoprops(glcm,'all');
        ccon=[ccon;stats.Contrast()];
        ccor=[ccor;stats.Correlation()];
        cen=[cen;stats.Energy()];    
        cho=[cho;stats.Homogeneity()];
        if i==length(file)
            contrast(:,:,x)=ccon; correlation(:,:,x)=ccor;
            energy(:,:,x)=cen; homogeneity(:,:,x)=cho;
        end
    end
end

bg_vals=[];
for i=1:9
%     avg=[contrast(:,:,i) correlation(:,:,i) energy(:,:,i) homogeneity(:,:,i)];
    avg=[mean(contrast(:,:,i)) mean(correlation(:,:,i)) mean(energy(:,:,i)) mean(homogeneity(:,:,i))];
    bg_vals=[bg_vals;avg];
end
plot(bg_vals(:,1),bg_vals(:,3),'b*','MarkerSize',10)
hold on
bg_vals=[mean(bg_vals(:,1)) mean(bg_vals(:,2)) mean(bg_vals(:,3)) mean(bg_vals(:,4))]