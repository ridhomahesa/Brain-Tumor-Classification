clc;clear;close all;warning off all;
 
image_folder = 'BPN Train';
filenames = dir(fullfile(image_folder, '*.png'));
total_images = numel(filenames);

for n = 1:total_images
    full_name = fullfile(image_folder, filenames(n).name);
    img = imread(full_name);
    img_gray=rgb2gray(img);
    
    cc = medfilt2(img_gray);
    
    T = 155;
    bw = im2bw(cc,T/255);

    %fs = get(0,'ScreenSize');
    %figure('Position',[0 0 fs(3)/2 fs(4)])

    SE = strel('disk',2);
    bw1 = imerode(bw,SE);

    %SE = strel('disk',0);
    bw2 = imdilate(bw1,SE);

    SE2 = strel('disk',4);
    bw3 = imerode(bw2,SE2);

    SE3 = strel('disk',4);

    %SE2 = strel('disk',1);
    bw4 = imdilate(bw3,SE3);

    % SE3 = strel('disk',6);
    % bw5 = imerode(bw4,SE3);
    % figure(9), imshow(bw5);

    % SE3 = strel('disk',1);
    % bw5a = imdilate(bw5,SE3);
    % figure(99), imshow(bw5a);

    bw6 = bwareaopen(bw4,350,8);
    %SE4 = strel('disk',2);
    %bw6 = imerode(bw5,SE4);

    bw6a = imfill(bw6, 'holes');

    SE5 = strel('disk',3);
    bw7 = imdilate(bw6a,SE5);

    if bw7 == 0
        brain1 = 255*uint8(bw7);
        brain_glcm = double(brain1);
        brain_glcm(brain_glcm==0) = NaN;
    else
        cc2 = cc;
        cc2(~bw7)=0;

        cc_resize = imresize(cc2, [256 256]);
        brain1 = cc_resize(cc_resize>0);
        brain_glcm = double(cc_resize);
        brain_glcm(brain_glcm==0) = NaN;
    end

    [N,M,L] = size(brain1);

    His2 =imhist(brain1)/(N*M); %Histogram Ternormalisasi

    Mean(n) = 0;
    for zi = 0:255
        Mean(n) = Mean(n) + (zi*His2(zi+1)); %Mean
    end

    Std2(n) = 0;
    for zi = 0:255
        Std2(n) = Std2(n) + ((zi-Mean(n)).^2*His2(zi+1));
    end
    Std2(n) = sqrt(Std2(n)); %Standard Deviation

    Ent(n) = 0;
    for zi = 0:255
        if His2(zi+1)>0
            Ent(n) = Ent(n) - (His2(zi+1)*log2(His2(zi+1))); %Entropy
        end
    end

    %Ekstraksi Ciri Tekstur Orde Dua (GLCM)
    jarak = 1;
    warning('off','Images:graycomatrix:scaledImageContainsNan');
    GLCM = graycomatrix(brain_glcm,'NumLevels',8, 'GrayLimits',[], 'Offset',[0 jarak; -jarak jarak; -jarak 0; -jarak -jarak]);
    stats = graycoprops(GLCM,{'contrast','homogeneity'});
    warning('on','Images:graycomatrix:scaledImageContainsNan');

    contrast = stats.Contrast;
    contrast0(n) = contrast(1);
    contrast45(n) = contrast(2);
    contrast90(n) = contrast(3);
    contrast135(n) = contrast(4);
    homogeneity = stats.Homogeneity;
    homogeneity0(n) = homogeneity(1);
    homogeneity45(n) = homogeneity(2);
    homogeneity90(n) = homogeneity(3);
    homogeneity135(n) = homogeneity(4);
    
end

input = [Mean;Std2;Ent;contrast0;contrast45;contrast90;contrast135;homogeneity0;homogeneity45;homogeneity90;homogeneity135]
[~,N] = size(input);
t1 = [1;0;0;0]; %Glioma
t2 = [0;1;0;0]; %Meningioma
t3 = [0;0;1;0]; %Metastatic
t4 = [0;0;0;1]; %Normal
target = [t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t2 t2 t2 t3 t3 t3 t3 t3 t3 t3 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4]

rng(27); %yg ini traingd

a = 0;
b = 255;
ra = 0.9;
rb = 0.1;
pa = (((ra-rb) * (input - a)) / (b - a)) + rb

net = newff(pa,target,2,{'logsig','tansig'},'traingd');

bobot_hidden1 = net.IW{1,1}
bias_hidden1 = net.b{1,1}
bobot_keluaran1 = net.LW{2,1}
bias_keluaran1 = net.b{2,1}

net.trainParam.epochs = 25000;
net.performFcn = 'mse';
net.trainParam.lr = 0.5;
% net.trainParam.mc = 0.95;
net.trainParam.goal = 0.02;
% net.trainParam.max_fail = 20000;
net.divideFcn = '';
[net, tr, Y, ee] = train(net,pa,target);

output = sim(net,pa)
trainIndices = vec2ind(output)

hasilbobot_hidden = net.IW{1,1}
hasilbias_hidden = net.b{1,1}
hasilbobot_keluaran = net.LW{2,1}
hasilbias_keluaran = net.b{2,1}
jumlah_iterasi = tr.num_epochs
nilai_keluaran = Y
nilai_error = ee
D = abs(nilai_error).^2;
error_MSE = sum(D(:))/numel(nilai_error)
plotperform(tr)

save jst.mat net
