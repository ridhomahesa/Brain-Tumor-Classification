function varargout = brain_tumor_classify(varargin)
% BRAIN_TUMOR_CLASSIFY MATLAB code for brain_tumor_classify.fig
%      BRAIN_TUMOR_CLASSIFY, by itself, creates a new BRAIN_TUMOR_CLASSIFY or raises the existing
%      singleton*.
%
%      H = BRAIN_TUMOR_CLASSIFY returns the handle to a new BRAIN_TUMOR_CLASSIFY or the handle to
%      the existing singleton*.
%
%      BRAIN_TUMOR_CLASSIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAIN_TUMOR_CLASSIFY.M with the given input arguments.
%
%      BRAIN_TUMOR_CLASSIFY('Property','Value',...) creates a new BRAIN_TUMOR_CLASSIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before brain_tumor_classify_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to brain_tumor_classify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help brain_tumor_classify

% Last Modified by GUIDE v2.5 13-Jul-2017 09:33:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @brain_tumor_classify_OpeningFcn, ...
                   'gui_OutputFcn',  @brain_tumor_classify_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before brain_tumor_classify is made visible.
function brain_tumor_classify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to brain_tumor_classify (see VARARGIN)

% Choose default command line output for brain_tumor_classify
handles.output = hObject;

% Update handles structure
setappdata(0,'datacontainer',hObject);
movegui(hObject,'onscreen')% To display application onscreen
movegui(hObject,'center')  % To display application in the center of screen


% UIWAIT makes brain_tumor_classify wait for user response (see UIRESUME)
% uiwait(handles.figure1);
clear
clc
warning off all
% --- Outputs from this function are returned to the command line.
function varargout = brain_tumor_classify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% INPUT FILE CITRA
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[nama_file,nama_path] = uigetfile({'*.png*'});
if ~isequal(nama_file,0)
    img = imread(fullfile(nama_path,nama_file));
    axes(handles.axes1)
    imshow(img)
    handles.img = img;
    guidata(hObject,handles)
else
    return
end
textLabel = sprintf('%s ',nama_path);
set(handles.edit15, 'String', textLabel);

% PENGUJIAN
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global kelas;

img = handles.img;

% Preprocessing
img_gray=rgb2gray(img);
cc = medfilt2(img_gray);
T = 155;
bw = im2bw(cc,T/255);

% Morphological Operation
SE = strel('disk',2);
bw1 = imerode(bw,SE);

bw2 = imdilate(bw1,SE);

SE2 = strel('disk',4);
bw3 = imerode(bw2,SE2);

SE3 = strel('disk',4);

bw4 = imdilate(bw3,SE3);

bw6 = bwareaopen(bw4,350,8);

bw6a = imfill(bw6, 'holes');
             
SE5 = strel('disk',3);
bw7 = imdilate(bw6a,SE5);

% Jika tidak ada yang tersegmentasi
if bw7 == 0
    cc_resize = 255*uint8(bw7);
    brain1 = cc_resize;
    brain_glcm = double(cc_resize);
    brain_glcm(brain_glcm==0) = NaN;
else
    cc2 = cc;
    cc2(~bw7)=0; % Mengembalikan ke grayscale

    cc_resize = imresize(cc2, [256 256]);
    brain1 = cc_resize(cc_resize>0);
    brain_glcm = double(cc_resize);
    brain_glcm(brain_glcm==0) = NaN;
end

axes(handles.axes11);
imshow(cc_resize);

[N,M,L] = size(brain1);

His2 =imhist(brain1)/(N*M); % Normalized Histogram

Mean = 0;
for zi = 0:255
    Mean = Mean + (zi*His2(zi+1)); %Mean
end
set(handles.edit1,'String',Mean);

Std2 = 0;
for zi = 0:255
    Std2 = Std2 + ((zi-Mean).^2*His2(zi+1));
end
Std2 = sqrt(Std2); %Standard Deviation
set(handles.edit2,'String',Std2);

Ent = 0;
for zi = 0:255
    if His2(zi+1)>0
        Ent = Ent - (His2(zi+1)*log2(His2(zi+1))); %Entropy
    end
end
set(handles.edit3,'String',Ent);

% GLCM Texture Feature Extraction
jarak = 1; %distance beetwen pixel
warning('off','Images:graycomatrix:scaledImageContainsNan');
GLCM = graycomatrix(brain_glcm,'NumLevels',8, 'GrayLimits',[], 'Offset',[0 jarak; -jarak jarak; -jarak 0; -jarak -jarak]);
stats = graycoprops(GLCM,{'contrast','homogeneity'});
warning('on','Images:graycomatrix:scaledImageContainsNan');
contrast = stats.Contrast;
contrast0 = contrast(1);
set(handles.edit19,'String',contrast0);
contrast45 = contrast(2);
set(handles.edit20,'String',contrast45);
contrast90 = contrast(3);
set(handles.edit21,'String',contrast90);
contrast135 = contrast(4);
set(handles.edit22,'String',contrast135);
homogeneity = stats.Homogeneity;
homogeneity0 = homogeneity(1);
set(handles.edit4,'String',homogeneity0);
homogeneity45 = homogeneity(2);
set(handles.edit17,'String',homogeneity45);
homogeneity90 = homogeneity(3);
set(handles.edit5,'String',homogeneity90);
homogeneity135 = homogeneity(4);
set(handles.edit18,'String',homogeneity135);

input = [Mean;Std2;Ent;contrast0;contrast45;contrast90;contrast135;homogeneity0;homogeneity45;homogeneity90;homogeneity135]
a = 0;
b = 255;
ra = 0.9;
rb = 0.1;
norm = (((ra-rb) * (input - a)) / (b - a)) + rb

load jst

output = sim(net,norm)
keluaran = vec2ind(output) % Find largest element in vector

if keluaran == 1
    kelas = 'Glioma'; %Target Class
    set(handles.edit6,'String',kelas);
    set(handles.edit6,'ForegroundColor','blue');
elseif keluaran == 2
    kelas = 'Meningioma';
    set(handles.edit6,'String',kelas);
    set(handles.edit6,'ForegroundColor','red');
elseif keluaran == 3
    kelas = 'Metastatic Bronchogenic Carcinoma';
    set(handles.edit6,'String',kelas);
    set(handles.edit6,'ForegroundColor','green');
elseif keluaran == 4
    kelas = 'Normal';
    set(handles.edit6,'String',kelas);
    set(handles.edit6,'ForegroundColor','black');
else
    kelas = 'Unknown';
    set(handles.edit6,'String',kelas);
    set(handles.edit6,'ForegroundColor','cyan');
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% PELATIHAN
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = uigetdir(pwd, 'Select a folder');
files = dir(fullfile(d, '*.png'));
textLabel = sprintf('%s ',d);
set(handles.edit14, 'String', textLabel);
total_images = numel(files)

for n = 1:total_images
    full_name = fullfile(d, files(n).name);
    
    %Prerpocessing
    img = imread(full_name);
    img_gray=rgb2gray(img);
    
    cc = medfilt2(img_gray);
    
    %Tumor Segmentation
    T = 155;
    bw = im2bw(cc,T/255);

    SE = strel('disk',2);
    bw1 = imerode(bw,SE);

    bw2 = imdilate(bw1,SE);

    SE2 = strel('disk',4);
    bw3 = imerode(bw2,SE2);

    SE3 = strel('disk',4);

    bw4 = imdilate(bw3,SE3);

    bw6 = bwareaopen(bw4,350,8);

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

    His2 =imhist(brain1)/(N*M); % Normalized Histogram

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

    %%Gray Level Coocurence Matrix (GLCM) Texture Feature Extraction
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

rng(27); 
a = 0;
b = 255;
ra = 0.9;
rb = 0.1;
pa = (((ra-rb) * (input - a)) / (b - a)) + rb

%train using BPNN
net = newff(pa,target,str2double(get(handles.edit27,'String')),{'logsig','tansig'},'traingd'); %Desain Jaringan
% net.IW{1,1} = [-0.1920,-0.4876,-1.1113,1.3023,1.1390;...
%     1.7783,-0.5090,-0.6322,0.8257,0.0001;...
%     0.1504,1.5682,0.5776,1.2974,-0.0645;...
%     0.0551,-1.2048,-0.9473,1.0189,1.0549;...
%     -1.0248,1.4706,0.8446,-0.6332,0.4195;...
%     -0.0383,1.4274,-1.3556,0.6851,0.4056;...
%     0.4467,1.0665,0.5537,-1.0882,1.2943;...
%     0.5339,-1.1961,-0.0174,-1.3992,0.9105];
%  
% net.b{1,1} = [2.1220;-1.5157;-0.9094;-0.3031;-0.3031;-0.9094;1.5157;2.1220];
%   
% net.LW{2,1} = [0.1534,-0.6342,-0.5201,0.7730,-0.9427,-0.0202,-0.6641,0.9574];
%  
% net.b{2,1} = [0.4254];

bobot_hidden1 = net.IW{1,1}
bias_hidden1 = net.b{1,1}
bobot_keluaran1 = net.LW{2,1}
bias_keluaran1 = net.b{2,1}

net.trainParam.showWindow = true;
net.trainParam.showCommandLine = true;
net.trainParam.epochs = str2double(get(handles.edit25,'String'));
net.performFcn = 'mse';
net.trainParam.lr = str2double(get(handles.edit26,'String'));
net.trainParam.goal = 0.02;
net.divideFcn = '';
[net, tr, Y, ee] = train(net,pa,target);
output = sim(net,pa)
trainIndices = vec2ind(output)

hasilbobot_hidden = net.IW{1,1}
hasilbias_hidden = net.b{1,1}
hasilbobot_keluaran = net.LW{2,1}
hasilbias_keluaran = net.b{2,1}
jumlah_iterasi = tr.num_epochs %number of iteration
nilai_keluaran = Y %output value
nilai_error = ee %error value
D = abs(nilai_error).^2;
error_MSE = sum(D(:))/numel(nilai_error)
set(handles.edit24,'String',error_MSE);
save jst.mat net

function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
