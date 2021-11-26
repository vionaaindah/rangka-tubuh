function varargout = rangka(varargin)
% RANGKA MATLAB code for rangka.fig
%      RANGKA, by itself, creates a new RANGKA or raises the existing
%      singleton*.
%
%      H = RANGKA returns the handle to a new RANGKA or the handle to
%      the existing singleton*.
%
%      RANGKA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RANGKA.M with the given input arguments.
%
%      RANGKA('Property','Value',...) creates a new RANGKA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rangka_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rangka_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rangka

% Last Modified by GUIDE v2.5 26-Nov-2021 10:38:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rangka_OpeningFcn, ...
                   'gui_OutputFcn',  @rangka_OutputFcn, ...
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


% --- Executes just before rangka is made visible.
function rangka_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rangka (see VARARGIN)

% Choose default command line output for rangka
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rangka wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rangka_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[nama_file,nama_path] = uigetfile(...
    {'*.jpg';'*.bmp';'*.png';'*.tif';},...
    'Pilih Citra');
if ~isequal(nama_file,0)
    handles.data1 = imread(fullfile(nama_path,nama_file));
    guidata(hObject, handles);
    axes(handles.citra);
    imshow(handles.data1);
    title('Citra Asli');
else
    return
end

% --- Executes on button press in result.
function result_Callback(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Img = handles.data1;
[tinggi, lebar] = size(Img);
ambang = 210; % Nilai ini bisa diubah-ubah
biner = zeros(tinggi, lebar);
for baris=1 : tinggi
    for kolom=1 : lebar
        if Img(baris, kolom) >= ambang
            Biner(baris, kolom) = 0;
        else
            Biner(baris, kolom) = 1;
        end
    end
end

img = Biner;
axes(handles.proses1);
handles.data4 = img;
imshow(handles.data4);
title('Citra Biner');
% Template Type
T1 = [0,1,1]';
T2 = T1';
T3 = [ 1,1,0]';
T4 = T3';

[row col plane] = size(img);
imgBin = double(im2bw(img));
ouImg = imgBin;
S = [1 3 5 7];
checKVal = 2;
template = 1;
outBinary  =zeros(row, col);
con = 5;

while (con < 15)
    for i = 2:row-1
        for j = 2:col-1
            window = imgBin(i-1:i+1,j-1:j+1);
            if (template == 1)
                andOp1 = isequal(window(:,2),T1);
                matchTemplate = andOp1;
            end
            if (template == 2)
                andOp1 = isequal(window(2,:),T2);
                matchTemplate = andOp1;
            end
            if (template == 3)
                andOp1 = isequal(window(:,2),T3);
                matchTemplate = andOp1;
            end
            if (template == 4)
                andOp1 = isequal(window(2,:),T4);
                matchTemplate = andOp1;
            end
            %         Connectivity number
            a=2;b=2;
            N0 = window(a,b);
            N1 = window(a,b+1);
            N2 = window(a-1,b+1);
            N3 = window(a-1,b);
            N4 = window(a-1,b-1);
            N5 = window(a,b-1);
            N6 = window(a+1,b-1);
            N7 = window(a+1,b);
            N8 = window(a+1,b+1);
            
            % arr = [(N1-(N3*N1*N2)) (N3-(N5*N3*N4)) (N5-(N7*N5*N6)) (N7-(N1*N7*N8))];
            
            val = [N1 N2 N3 N4 N5 N6 N7 N8];
            arr = [N1<N2 N2<N3 N3<N4 N4<N5 N5<N6 N6<N7 N7<N8 N8<N1];
            Cn = sum(arr);
            EndPoint = abs(N0 - sum(val));

            if imgBin(i,j)==1
                if ((matchTemplate) == 1)
                    if Cn == 1
                        if (EndPoint ~= 0)
                            outBinary(i,j) = 1;
                        end
                    end
                end
            end
        end
    end
    checKVal = sum(sum(outBinary));
    if (checKVal==0)
        con = con+1;
    end

    binVal = find(outBinary==1);
    imgBin(binVal) = 0;
    axes(handles.proses2);
    handles.data3 = outBinary,[];
    imshow(handles.data3);
    title('Proses');
    outBinary  =zeros(row, col);
    template = template+1;

    %     Iteration
    if template == 5
        axes(handles.hasil);
        handles.data5 = imgBin,[];
        imshow(handles.data5);
        title('Bentuk Rangka');
        outBinary  =zeros(row, col);
        template = 1;
    end
end
axes(handles.proses1);
handles.data4 = img;
imshow(handles.data4);
title('Citra Biner');

axes(handles.hasil);
handles.data5 = imgBin,[];
imshow(handles.data5);
title('Bentuk Rangka');

F = handles.data5;
H = ones(4);

[th, lh]=size(H);
[tf, lf]=size(F);
hotx = round(lh/2);
hoty = round(th/2);
 
Xh = [];
Yh = [];
jum_anggota = 0;
% Menentukan koordinat piksel bernilai 1 pada H
for baris = 1 : th
     for kolom = 1 : lh
          if H(baris, kolom) == 1
              jum_anggota = jum_anggota + 1;
              Xh(jum_anggota) = -hotx + kolom;
              Yh(jum_anggota) = -hoty + baris;
          end
     end
end
G = zeros(tf, lf); % Nolkan semua pada hasil dilasi
% Memproses dilasi
for baris = 1 : tf
    for kolom = 1 : lf
        for indeks = 1 : jum_anggota
            if F(baris, kolom) == 1
                xpos = kolom + Xh(indeks);
                ypos = baris + Yh(indeks);
                if (xpos >= 1) && (xpos <= lf) && ...
                        (ypos >= 1) && (ypos <= tf)
                    G(ypos, xpos) = 1;
                end
            end
        end
    end
end
axes(handles.hasil);
handles.data6 = G;
imshow(handles.data6);
title('Bentuk Rangka');

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.citra,'reset');
cla(handles.hasil,'reset');
cla(handles.proses1,'reset');
cla(handles.proses2,'reset');


% --- Executes on button press in biner.
function biner_Callback(hObject, eventdata, handles)
% hObject    handle to biner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Img = handles.data1;
[tinggi, lebar] = size(Img);
ambang = 210; % Nilai ini bisa diubah-ubah
biner = zeros(tinggi, lebar);
for baris=1 : tinggi
    for kolom=1 : lebar
        if Img(baris, kolom) >= ambang
            Biner(baris, kolom) = 0;
        else
            Biner(baris, kolom) = 1;
        end
    end
end
axes(handles.hasil);
handles.data2 = Biner;
imshow(handles.data2);
title('Citra Biner');

% --- Executes on button press in thinning.
function thinning_Callback(hObject, eventdata, handles)
img = handles.data1;
% resize the image
axes(handles.proses1);
handles.data3 = img;
imshow(handles.data3);
title('Input Image');
% Template Type
T1 = [0,1,1]';
T2 = T1';
T3 = [ 1,1,0]';
T4 = T3';

[row col plane] = size(img);
imgBin = double(im2bw(img));
ouImg = imgBin;
S = [1 3 5 7];
checKVal = 2;
template = 1;
outBinary  =zeros(row, col);
con = 5;

while (con < 15)
    for i = 2:row-1
        for j = 2:col-1
            window = imgBin(i-1:i+1,j-1:j+1);
            if (template == 1)
                andOp1 = isequal(window(:,2),T1);
                matchTemplate = andOp1;
            end
            if (template == 2)
                andOp1 = isequal(window(2,:),T2);
                matchTemplate = andOp1;
            end
            if (template == 3)
                andOp1 = isequal(window(:,2),T3);
                matchTemplate = andOp1;
            end
            if (template == 4)
                andOp1 = isequal(window(2,:),T4);
                matchTemplate = andOp1;
            end
            %         Connectivity number
            a=2;b=2;
            N0 = window(a,b);
            N1 = window(a,b+1);
            N2 = window(a-1,b+1);
            N3 = window(a-1,b);
            N4 = window(a-1,b-1);
            N5 = window(a,b-1);
            N6 = window(a+1,b-1);
            N7 = window(a+1,b);
            N8 = window(a+1,b+1);
            
            % arr = [(N1-(N3*N1*N2)) (N3-(N5*N3*N4)) (N5-(N7*N5*N6)) (N7-(N1*N7*N8))];
            
            val = [N1 N2 N3 N4 N5 N6 N7 N8];
            arr = [N1<N2 N2<N3 N3<N4 N4<N5 N5<N6 N6<N7 N7<N8 N8<N1];
            Cn = sum(arr);
            EndPoint = abs(N0 - sum(val));

            if imgBin(i,j)==1
                if ((matchTemplate) == 1)
                    if Cn == 1
                        if (EndPoint ~= 0)
                            outBinary(i,j) = 1;
                        end
                    end
                end
            end
        end
    end
    checKVal = sum(sum(outBinary));
    if (checKVal==0)
        con = con+1;
    end

    binVal = find(outBinary==1);
    imgBin(binVal) = 0;
    axes(handles.proses2);
    handles.data2 = outBinary,[];
    imshow(handles.data2);
    title('Proses');
    outBinary  =zeros(row, col);
    template = template+1;

    %     Iteration
    if template == 5
        axes(handles.hasil);
        handles.data4 = imgBin,[];
        imshow(handles.data4);
        title('Thinning');
        outBinary  =zeros(row, col);
        template = 1;
    end
end
axes(handles.proses1);
handles.data3 = img;
imshow(handles.data3);
title('Input Image');

axes(handles.hasil);
handles.data4 = imgBin,[];
imshow(handles.data4);
title('Thinning');


% --- Executes on button press in dilasi.
function dilasi_Callback(hObject, eventdata, handles)
% hObject    handle to dilasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.data1;
F = im2bw(img, 0.5);
H = ones(4);

[th, lh]=size(H);
[tf, lf]=size(F);
hotx = round(lh/2);
hoty = round(th/2);
 
Xh = [];
Yh = [];
jum_anggota = 0;
% Menentukan koordinat piksel bernilai 1 pada H
for baris = 1 : th
     for kolom = 1 : lh
          if H(baris, kolom) == 1
              jum_anggota = jum_anggota + 1;
              Xh(jum_anggota) = -hotx + kolom;
              Yh(jum_anggota) = -hoty + baris;
          end
     end
end
G = zeros(tf, lf); % Nolkan semua pada hasil dilasi
% Memproses dilasi
for baris = 1 : tf
    for kolom = 1 : lf
        for indeks = 1 : jum_anggota
            if F(baris, kolom) == 1
                xpos = kolom + Xh(indeks);
                ypos = baris + Yh(indeks);
                if (xpos >= 1) && (xpos <= lf) && ...
                        (ypos >= 1) && (ypos <= tf)
                    G(ypos, xpos) = 1;
                end
            end
        end
    end
end
axes(handles.hasil);
handles.data2 = G;
imshow(handles.data2);
title('Hasil Dilasi');
