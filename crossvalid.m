clc;clear;close all;warning off all;
 
image_folder = 'Crossvalid';
filenames = dir(fullfile(image_folder, '*.png'));
total_images = numel(filenames);

for n = 1:total_images
    full_name = fullfile(image_folder, filenames(n).name);
    img = imread(full_name);
    
    %Prerpocessing
    img_gray=rgb2gray(img);
    
    cc = medfilt2(img_gray);
    
    %Tumor Segmentation
    T = 155;
    bw = im2bw(cc,T/255);

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

    His2 =imhist(brain1)/(N*M); %Normalized Histogram

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

    % GLCM Texture Feature Extraction
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

input = [Mean;Std2;Ent;contrast0;contrast45;contrast90;contrast135;homogeneity0;homogeneity45;homogeneity90;homogeneity135];
[~,N] = size(input);
t1 = [1;0;0;0]; %Glioma
t2 = [0;1;0;0]; %Meningioma
t3 = [0;0;1;0]; %Metastatic
t4 = [0;0;0;1]; %Normal
target = [t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t2 t2 t2 t2 t2 t2 t3 t3 t3 t3 t3 t3 t3 t3 t3 t3 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4 t4];
target_ind = vec2ind(target);

a = 0;
b = 255;
ra = 0.9;
rb = 0.1;
pa = (((ra-rb) * (input - a)) / (b - a)) + rb;

k = 10; %Partition and Iteration Value
Indices = crossvalind('Kfold', 50, k); %using K-Fold Cross Validation
counter_acc = 0;
counter_err = 0;
counter_precision = 0;
counter_MSE = 0;
for i = 1:k
    test = (Indices == i);
    trainn = ~test;
    
    trIdn = find(trainn)
    tsIdn = find(test)
    
    rng(27);
    net = newff(pa(:,trIdn),target(:,trIdn),2,{'logsig','tansig'},'traingd');
    
    bobot_hidden1 = net.IW{1,1};
    bias_hidden1 = net.b{1,1};
    bobot_keluaran1 = net.LW{2,1};
    bias_keluaran1 = net.b{2,1};
    
    net.performFcn = 'mse';
    net.trainParam.lr = 0.2;
    net.trainParam.goal = 0.02;
    net.divideFcn = '';
    net.trainParam.epochs = 25000;
    
    [net, tr, Y, ee] = train(net,pa(:,trIdn),target(:,trIdn));
    output2 = sim(net,pa(:,trIdn));
    trainIndices2 = vec2ind(output2);

    hasilbobot_hidden = net.IW{1,1};
    hasilbias_hidden = net.b{1,1};
    hasilbobot_keluaran = net.LW{2,1};
    hasilbias_keluaran = net.b{2,1};
    jumlah_iterasi = tr.num_epochs; %Number of iteration
    nilai_keluaran = Y; %Network output value
    nilai_error = ee; %error value beetwen network and target
    D = abs(nilai_error).^2;
    error_MSE = sum(D(:))/numel(nilai_error)
    counter_MSE = counter_MSE + error_MSE
    
    save jst.mat net
    load jst
    cp = classperf(target_ind(:,tsIdn));
    output1 = sim(net,pa(:,tsIdn));
    testIndices = vec2ind(output1)
    
    plotconfusion(target(:,tsIdn),output1)
    output_trans = testIndices.'
    classperf(cp, output_trans, test(Indices==i))
    cp.CountingMatrix
    cp.CorrectRate * 100
    cp.ErrorRate * 100
    cp.PositivePredictiveValue * 100
    counter_acc = counter_acc + (cp.CorrectRate * 100)
    counter_err = counter_err + (cp.ErrorRate * 100)
    counter_precision = counter_precision + (cp.PositivePredictiveValue * 100)
end

average_accuracy = counter_acc/k
average_err = counter_err/k
% average_precision = counter_precision/k
average_MSE = counter_MSE/k