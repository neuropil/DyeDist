function varargout = DyeDist_dev2(varargin)
% DYEDIST_DEV2 MATLAB code for DyeDist_dev2.fig
%      DYEDIST_DEV2, by itself, creates a new DYEDIST_DEV2 or raises the existing
%      singleton*.
%
%      H = DYEDIST_DEV2 returns the handle to a new DYEDIST_DEV2 or the handle to
%      the existing singleton*.
%
%      DYEDIST_DEV2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYEDIST_DEV2.M with the given input arguments.
%
%      DYEDIST_DEV2('Property','Value',...) creates a new DYEDIST_DEV2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DyeDist_dev2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DyeDist_dev2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DyeDist_dev2

% Last Modified by GUIDE v2.5 07-Mar-2014 11:32:52


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% THINGS TO DO
% ADD NEW DISTRIBUTION ANALYSIS by quadrant (3/6/2014)



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DyeDist_dev2_OpeningFcn, ...
    'gui_OutputFcn',  @DyeDist_dev2_OutputFcn, ...
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
% JAT 11/19/2013 7:07pm

% --- Executes just before DyeDist_dev2 is made visible.
function DyeDist_dev2_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DyeDist_dev2 (see VARARGIN)

% Choose default command line output for DyeDist_dev2
handles.output = hObject;

set(handles.sectList,'Enable','off')
set(handles.imDisplay,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[]);

set(handles.gThresh,'Enable','off');
set(handles.roiButton,'Enable','off');

set(handles.imOpts,'Enable','off');

set(handles.meanT,'Visible','off')
set(handles.sdT,'Visible','off')
set(handles.thresT,'Visible','off')
set(handles.plusSign,'Visible','off')
set(handles.prodSign,'Visible','off')
set(handles.equalSign,'Visible','off')
set(handles.stdT,'Visible','off')

% CHECK IF DATASET ARRAYs are supported
versionCheck = version('-release');
getYear = str2double(versionCheck(1:end-1));
if getYear < 2013
    set(handles.expDATASET,'Enable','off')
end


% Create Column titles
handles.columnNames = {'ROI_Area','Dye_Area','Dye_Ratio','DyeOut_Area','DyeOut_Ratio'};
set(handles.dataTable,'ColumnName',handles.columnNames);

handles.figActive = 0;
handles.hemiCount = 0;

set(handles.infoT,'String','Load Image Files');

handles.SheetCount = 1;
handles.CaseNames = {};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DyeDist_dev2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DyeDist_dev2_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function fileOpts_Callback(~, ~, ~)
% hObject    handle to fileOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function loadFold_Callback(hObject, ~, handles)
% hObject    handle to loadFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.infoT,'String',[]);


handles.WD = uigetdir('C:\');
if ~handles.WD
    return
end

cd(handles.WD)

% Setup initial indicies
dirNames = dir('*.tif');
actualNames = {dirNames.name}';
handles.ImgNames = actualNames;

% Error hanlding %
if isempty(handles.ImgNames)
    warndlg('NO .tif files in folder','NO TIF Files');
    return
end

% User defined number of channels per image
chanNumQ = questdlg('How many channels per image?','Channel Number','1','2','3','3');
chanNum = str2double(chanNumQ);
% Run PARSEIMAGESFILES function
handles.imageINFO = parseImageFiles(handles.ImgNames,chanNum);
% Convert Section Numbers into vector
if isfield(handles.imageINFO,'Image')
    handles.numSections = cell2mat(handles.imageINFO.Image.Nums);
else
    handles.numSections = cell2mat(handles.imageINFO.ROI.Nums);
end
% Create Section Name Tags that will appear in list panel
handles.ListIndex = arrayfun(@(x) ['Section',' ', num2str(x)],handles.numSections, 'UniformOutput', false);
% Load names into list panel
set(handles.sectList,'String',handles.ListIndex);
% User defined RGB designations
chanCond = inputdlg({'Nissl','ROI'}, 'Channel ID', [1 50; 1 50],{'Blue','Red'});

chanCondID = {'Nissl','ROI'};
for fi2 = 1:2
    switch chanCond{fi2,1}
        case 'Red'
            handles.colorId(fi2) = 1;
            handles.colorIdtag{fi2} = chanCondID{fi2};
        case 'Green'
            handles.colorId(fi2) = 2;
            handles.colorIdtag{fi2} = chanCondID{fi2};
        case 'Blue'
            handles.colorId(fi2) = 3;
            handles.colorIdtag{fi2} = chanCondID{fi2};
    end
end

set(handles.infoT,'String',[]);
set(handles.gThresh,'Enable','on');
set(handles.sectList,'Enable','on')
cla(handles.imDisplay);
set(handles.infoT,'String','View sections or Select Start Analysis');

guidata(hObject, handles);



% --- Executes on selection change in sectList.
function sectList_Callback(hObject, ~, handles)
% hObject    handle to sectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sectList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sectList

secChoice = get(hObject,'Value');
handles.sectionNum2use = handles.numSections(secChoice);
handles.section2use = secChoice;

% handles.colorIdtag

if isfield(handles.imageINFO,'Image');
    imageToggle = 0;
    imageFname = handles.imageINFO.Image.FNames{handles.section2use};
    imForsize = imread(imageFname);
    [dim1,dim2,~] = size(imForsize);
else
    imageToggle = 1;
    nisslFname = handles.imageINFO.Nissl.FNames{handles.section2use};
    roiFname = handles.imageINFO.ROI.FNames{handles.section2use};
    imForsize = imread(nisslFname);
    [dim1,dim2] = size(imForsize);
end

image2Disp = uint8(zeros(dim1,dim2,3));
if imageToggle % 2 images each a different channel
    
    nisslMatrix = imread(nisslFname);
    image2Disp(:,:,handles.colorId(1)) = nisslMatrix;
    
    roiMatrix = imread(roiFname);
    image2Disp(:,:,handles.colorId(2)) = roiMatrix;
    
else % Single image 2 channels
    
    imageMatrix = imread(imageFname);
    image2Disp(:,:,handles.colorId(1)) = imageMatrix(:,:,handles.colorId(1));
    image2Disp(:,:,handles.colorId(2)) = imageMatrix(:,:,handles.colorId(2));
    
end

handles.image2show = image2Disp;

cla(handles.imDisplay);
axes(handles.imDisplay);
imshow(handles.image2show);



guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sectList_CreateFcn(hObject, ~, ~)
% hObject    handle to sectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roiButton.
function roiButton_Callback(hObject, ~, handles)
% hObject    handle to roiButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.roiButton,'Enable','off');

caseName = inputdlg('Input Case Name','Case Name',1,{'Case1'});
	
if isempty(caseName)
	    caseName = 'Case1';
end
	
% SET UP TOGGLE FOR BOX AND ALGORITHM
boxQuestion = 1;
while boxQuestion
    
    chooseThresh = questdlg('Measure distribution in Polygon only','Threshold?','Yes','No','Yes');
    
    switch chooseThresh
        case 'Yes' % ONLY USE POLYGON
            boxToggle = 0;
            boxQuestion = 0;
            handles.ExportColnames{handles.SheetCount} = handles.columnNames(1,1:3);
        case 'No'
            boxToggle = 1;
            boxQuestion = 0;
            handles.ExportColnames{handles.SheetCount} = handles.columnNames;
    end
    
end

% set(handles.fileOpts,'Enable','off');

set(handles.infoT,'String',[]);
cla(handles.imDisplay);
set(handles.sectList,'Enable','off');

set(handles.gThresh,'Enable','off');
set(handles.stdT,'Enable','off');

Xcoords = cell(length(handles.numSections),1);
Ycoords = cell(length(handles.numSections),1);
polyArea = zeros(length(handles.numSections),1);
injArea = zeros(length(handles.numSections),1);
oinjArea = zeros(length(handles.numSections),1);

dataTOout = struct;
handles.quadInfo = {};
axes(handles.imDisplay)

%%%%% FIGURE HOW TO KEEP TRACK OF MULTIPLE CASES/SHEETS

Quadrants = cell(1,length(handles.numSections));

for i = 1:length(handles.numSections)
    
    cd(handles.WD)
    
    set(handles.sectList,'Value',i);
    

    
    if isfield(handles.imageINFO,'Image');
        imageToggle = 0;
        imageFname = handles.imageINFO.Image.FNames{i};
        imForsize = imread(imageFname);
        [dim1,dim2,~] = size(imForsize);
    else
        imageToggle = 1;
        nisslFname = handles.imageINFO.Nissl.FNames{i};
        roiFname = handles.imageINFO.ROI.FNames{i};
        imForsize = imread(nisslFname);
        [dim1,dim2] = size(imForsize);
    end
    
    traceImage = uint8(zeros(dim1,dim2,3));
    injThreshImage = uint8(zeros(dim1,dim2,3));
    if imageToggle % 2 images each a different channel
        
        nisslMatrix = imread(nisslFname);
        traceImage(:,:,handles.colorId(1)) = nisslMatrix;
        
        roiMatrix = imread(roiFname);
        injThreshImage(:,:,handles.colorId(2)) = roiMatrix;
        injImage = roiMatrix;
        
    else % Single image 2 channels
        
        imageMatrix = imread(imageFname);
        traceImage(:,:,handles.colorId(1)) = imageMatrix(:,:,handles.colorId(1));
        injThreshImage(:,:,handles.colorId(2)) = imageMatrix(:,:,handles.colorId(2));
        injImage = imageMatrix(:,:,handles.colorId(2));
        
    end

    % Trace Image
    
    cla(handles.imDisplay)
    imshow(traceImage);
    
    set(handles.infoT,'String','Trace Region of Interest');
    [~, Xcoords{i,1}, Ycoords{i,1}] = roipoly(traceImage);    
    polyArea(i) = polyarea(Xcoords{i,1},Ycoords{i,1});
    
    % Create NTS mask
    mNtb_mask = poly2mask(Xcoords{i,1},Ycoords{i,1},dim1,dim2);

    if boxToggle % FIND END
	        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        % Get sample boxes
	        
	        totPixels = dim1*dim2;
	        onePercent = totPixels*0.005;
	        boxDim = round(sqrt(onePercent));
	        halfBD = round(boxDim/2);
	        
	        cla(handles.imDisplay); % REPEAT THIS LINE OUTSIDE TOGGLE
	        %     axes(handles.imDisplay); % REPEAT THIS LINE OUTSIDE TOGGLE
	        imshow(injThreshImage); % REPEAT THIS LINE OUTSIDE TOGGLE
        
	        set(handles.infoT,'String','Click three times on Image');
	        [x_coord, y_coord] = ginput(3);
	        set(handles.infoT,'String',[]);
	        x_coord = round(x_coord);
	        y_coord = round(y_coord);
	        
	        sYcoords = cell(1,3);
	        sXcoords = cell(1,3);
	        for si = 1:length(x_coord)
	            sXcoords{1,si} = round([x_coord(si) - halfBD , x_coord(si) + halfBD, x_coord(si) + halfBD, x_coord(si) - halfBD]);
	            sYcoords{1,si} = round([y_coord(si) - halfBD , y_coord(si) - halfBD, y_coord(si) + halfBD, y_coord(si) + halfBD]);
	        end
	        
	        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        % Extract data and plot boxes
	        
	        pixelsPerbox = cell(1,3);
	        hold on
	        boxPlotsS = cell(1,3);
	        for sti = 1:3
	            samplebox = roipoly(dim1,dim2,sYcoords{1,sti},sXcoords{1,sti});
	            square_mask = poly2mask(sXcoords{1,sti},sYcoords{1,sti},dim1,dim2);
	            pixelsPerbox{1,sti} = injThreshImage(square_mask);
	            [Bi, ~] = bwboundaries(samplebox,'noholes');
	            boxIndices = cell2mat(Bi);
	            boxPlotsS{1,sti} = boxIndices;
	        end
	        
	        % calculate threshold
	        num_inj_pixels = injImage(mNtb_mask);
	        convertPixels2double = single(num_inj_pixels);
	        
	        handles.pixelsBackground = double([pixelsPerbox{1,1};pixelsPerbox{1,2};pixelsPerbox{1,3};convertPixels2double]);
	        
	        handles.pMean = round(mean(handles.pixelsBackground));
            handles.pStd = round(std(handles.pixelsBackground));
            
            pixThresh = handles.pMean + (str2double(get(handles.stdT,'String'))*handles.pStd);
            
            % To ensure that Standard Deviation does not artifically
            % inflate pixel threshold in excess of 255
            if pixThresh >= 255
                pixThresh = pixThresh - (2*handles.pStd);
                
                % Double check of pixel threshold 2/23/2014
                while pixThresh >= 255
                    pixThresh = pixThresh - (handles.pStd/2);
                end
                
            end
            
    else
        cla(handles.imDisplay); 
        imshow(injThreshImage); 

        % calculate threshold
        num_inj_pixels = injImage(mNtb_mask);
        convertPixels2double = single(num_inj_pixels);
        
        handles.pixelsBackground =  double(convertPixels2double);
        
        handles.pMean = round(mean(handles.pixelsBackground));
        handles.pStd = round(std(handles.pixelsBackground));
        
        pixThresh = handles.pMean + (str2double(get(handles.stdT,'String'))*handles.pStd);
        
        modImage = injImage; % save image in new file
        modImage(~mNtb_mask) = 0; % use inverse of mask to get rid of out poly pixels
        threshExceedIndex = modImage > pixThresh;
        % to plot image of injection threshold match
        modImage(~threshExceedIndex) = 0;
        
        injArea(i) = bwarea(threshExceedIndex); % BASELINE 
        
        % CREATE BISECTED POLYGON FUNCTION
        
        [B,~,~,~] = bwboundaries(mNtb_mask);
        
        pixelInfo = regionprops(mNtb_mask,traceImage(:,:,handles.colorId(1)),'Centroid');
        
        % xmin = min(xCoords);
        xmax = max(Xcoords{i,1});
        % ymin = min(yCoords);
        ymax = max(Ycoords{i,1});
        
        % find corner vertices
        xCor = Xcoords{i,1}(2:end);
        yCor = Ycoords{i,1}(2:end);
        
        yCorR = zeros(numel(yCor),1);
        for yi = 1:numel(yCor)
            yCorR(yi) = yCor(yi) + rand;
        end
        
        yCor = yCorR;
        
        xCorR = zeros(numel(xCor),1);
        for xi = 1:numel(xCor)
            xCorR(xi) = xCor(xi) + rand;
        end
        
        xCor = xCorR;
        
        % Top Left corner
        topLeft = (yCor < ymax*0.5 & xCor < xmax*0.5);
        yratio = 0.5;
        xratio = 0.5;
        
        if sum(topLeft) ~= 1
            if sum(topLeft) == 0
                
                while sum(topLeft) ~= 1
                    yratio = yratio + 0.001;
                    yval = ymax*yratio;
                    xratio = xratio + 0.001;
                    xval = ymax*xratio;
                    if sum(yCor < yval) == 0
                        topLeft = (yCor == min(yCor) & xCor < xval);
                    elseif sum(xCor < xval) == 0
                        topLeft = (yCor < yval & xCor == min(xCor));
                    else
                        topLeft = (yCor < yval & xCor < xval);
                    end
                end
                
            elseif sum(topLeft) > 1 % too liberal
                while sum(topLeft) ~= 1
                    yratio = yratio - 0.001;
                    yval = ymax*yratio;
                    xratio = xratio - 0.001;
                    xval = ymax*xratio;
                    if sum(yCor < yval) == 0
                        topLeft = (yCor == min(yCor) & xCor < xval);
                    elseif sum(xCor < xval) == 0
                        topLeft = (yCor < yval & xCor == min(xCor));
                    else
                        topLeft = (yCor < yval & xCor < xval);
                    end
                end
            end
        end
        
        % Bottom Left corner
        botLeft = (yCor > ymax*0.5 & xCor < xmax*0.5);
        yratio = 0.5;
        xratio = 0.5;
        
        if sum(botLeft) ~= 1
            if sum(botLeft) == 0 % too conservative
                while sum(botLeft) ~= 1
                    yratio = yratio - 0.001;
                    yval = ymax*yratio;
                    xratio = xratio + 0.001;
                    xval = ymax*xratio;
                    if sum(yCor > yval) == 0
                        botLeft = (yCor == max(yCor) & xCor < xval);
                    elseif sum(xCor < xval) == 0
                        botLeft = (yCor > yval & xCor == min(xCor));
                    else
                        botLeft = (yCor > yval & xCor < xval);
                    end
                end
            elseif sum(botLeft) > 1 % too liberal
                while sum(botLeft) ~= 1
                    yratio = yratio + 0.001;
                    yval = ymax*yratio;
                    xratio = xratio - 0.001;
                    xval = ymax*xratio;
                    if sum(yCor > yval) == 0
                        botLeft = (yCor == max(yCor) & xCor < xval);
                    elseif sum(xCor < xval) == 0
                        botLeft = (yCor > yval & xCor == min(xCor));
                    else
                        botLeft = (yCor > yval & xCor < xval);
                    end
                end
            end
        end
        
        % Bottom Right corner
        botRight = (yCor > ymax*0.5 & xCor > xmax*0.5);
        yratio = 0.5;
        xratio = 0.5;
        
        if sum(botRight) ~= 1
            if sum(botRight) == 0
                while sum(botRight) ~= 1 % too conservative
                    yratio = yratio - 0.01;
                    xratio = xratio - 0.01;
                    botRight = (yCor > ymax*yratio & xCor > xmax*xratio);
                end
            elseif sum(botRight) > 1
                while sum(botRight) ~= 1 % too liberal
                    yratio = yratio + 0.01;
                    xratio = xratio + 0.01;
                    botRight = (yCor > ymax*yratio & xCor > xmax*xratio);
                end
            end
        end
        
        % Top Right corner
        topRight = (yCor < ymax*0.5 & xCor > xmax*0.5);
        yratio = 0.5;
        xratio = 0.5;
        
        if sum(topRight) ~= 1
            if sum(topRight) == 0
                
                if sum(yCor < ymax*0.5) == 0
                    
                    maxSortx = sort(xCor,'descend');
                    tri = 1;
                    while sum(topRight) ~= 1 % too conservative
                        maxNow = maxSortx(1:tri);
                        xVals = ismember(xCor, maxNow);
                        yratio = yratio + 0.01;
                        topRight = (xVals & yCor < ymax*yratio);
                        tri = 1 + 1;
                    end
                    
                else
                    
                    while sum(topRight) ~= 1 % too conservative
                        yratio = yratio + 0.01;
                        xratio = xratio - 0.01;
                        topRight = (yCor < ymax*yratio & xCor > xmax*xratio);
                    end
                    
                end
            elseif sum(topRight) > 1

                while sum(topRight) ~= 1 % too liberal
                    yratio = yratio - 0.001;
                    xratio = xratio + 0.001;
                    topRight = (yCor < ymax*yratio & xCor > xmax*xratio);
                end
            end
        end
       
        
        % Top Coordinate
        topMidDist = (xCor(topRight) + xCor(topLeft))/2;
        topCoord = min(find(B{1,1}(:,2) == floor(topMidDist) & B{1,1}(:,1) < ymax*0.6));
        
        % Bottom Coordinate
        botMidDist = (xCor(botRight) + xCor(botLeft))/2;
        botCoord = max(find(B{1,1}(:,2) == floor(botMidDist) & B{1,1}(:,1) > ymax*0.6));
        
        % Left Coordinate
        leftMidDist = (yCor(topLeft) + yCor(botLeft))/2;
        leftCoord = min(find(B{1,1}(:,1) == floor(leftMidDist) & B{1,1}(:,2) < xmax*0.6));
        
        % Right Coordinate
        rightMidDist = (yCor(topRight) + yCor(botRight))/2;
        rightCoord = max(find(B{1,1}(:,1) == ceil(rightMidDist) & B{1,1}(:,2) > xmax*0.6));
        
        % Quadrant 1
        
        % Quadrant 1 will be the TOP LEFT QUADRANT
        
        % Set the three anchor vertices (centriod, top , left)
        q1_coord1 = [pixelInfo.Centroid(1) ; pixelInfo.Centroid(2)];
        q1_coord2 = [(B{1,1}(leftCoord,2));(B{1,1}(leftCoord,1))];
        q1_coord3 = [(B{1,1}(topCoord,2));(B{1,1}(topCoord,1))];
        
        % Create X and Y coordinates out of vertices
        quadrant_1_Xcoords = [q1_coord1(1);q1_coord2(1);q1_coord3(1)];
        quadrant_1_Ycoords = [q1_coord1(2);q1_coord2(2);q1_coord3(2)];
        
        q1_test = ((yCor < quadrant_1_Ycoords(1) & yCor < quadrant_1_Ycoords(2)) &...
            (xCor < quadrant_1_Xcoords(1) & xCor < quadrant_1_Xcoords(3)));
        
        % Get indices for found coordinates
        Quad1_coord_indices = find(q1_test == 1);
        Quad1Flip = flipud(Quad1_coord_indices);
        
        % Insert found vertices into quadrant coordinates
        quad1Xtemp = quadrant_1_Xcoords(1:2);
        quad1Ytemp = quadrant_1_Ycoords(1:2);
        quad1Xout = quad1Xtemp;
        quad1Yout = quad1Ytemp;
        for q1i2 = 1:length(Quad1Flip);
            addXquad1 = xCor(Quad1Flip(q1i2));
            addYquad1 = yCor(Quad1Flip(q1i2));
            quad1Xout = [quad1Xout ; addXquad1];
            quad1Yout = [quad1Yout ; addYquad1];
        end
        
        quad1Xfinal = [quad1Xout ; quadrant_1_Xcoords(3)];
        quad1Yfinal = [quad1Yout ; quadrant_1_Ycoords(3)];
        
        quadCoords.TL = [quad1Xfinal quad1Yfinal];
        
        handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,1} = quadCoords.TL;
%         h = impoly(injImage, quad_1);
        
        % Quadrant 1 Mask
%         Quadrants.TL = createMask(h);
        
        % Quadrant 2
        
        % Quadrant 2 will be the BOTTOM LEFT QUADRANT
        
        q2_coord1 = [pixelInfo.Centroid(1) ; pixelInfo.Centroid(2)];
        q2_coord2 = [(B{1,1}(botCoord,2));(B{1,1}(botCoord,1))];
        q2_coord3 = [(B{1,1}(leftCoord,2));(B{1,1}(leftCoord,1))];
        
        quadrant_2_Xcoords = [q2_coord1(1);q2_coord2(1);q2_coord3(1)];
        quadrant_2_Ycoords = [q2_coord1(2);q2_coord2(2);q2_coord3(2)];
        
        q2_test = ((yCor > quadrant_2_Ycoords(1) & yCor > quadrant_2_Ycoords(3)) &...
            (xCor < quadrant_2_Xcoords(1) & xCor < quadrant_2_Xcoords(2)));
        
        Quad2_coord_indices = find(q2_test == 1);
        
        % Insert found vertices into quadrant coordinates
        quad2Xtemp = quadrant_2_Xcoords(1:2);
        quad2Ytemp = quadrant_2_Ycoords(1:2);
        quad2Xout = quad2Xtemp;
        quad2Yout = quad2Ytemp;
        for q1i2 = 1:length(Quad2_coord_indices);
            addXquad2 = xCor(Quad2_coord_indices(q1i2));
            addYquad2 = yCor(Quad2_coord_indices(q1i2));
            quad2Xout = [quad2Xout ; addXquad2];
            quad2Yout = [quad2Yout ; addYquad2];
        end
        
        quad2Xfinal = [quad2Xout ; quadrant_2_Xcoords(3)];
        quad2Yfinal = [quad2Yout ; quadrant_2_Ycoords(3)];
        
        quadCoords.BL = [quad2Xfinal quad2Yfinal];
        
        % Create polygon handles for quadrant 2
        
%         h2 = impoly(gca, quad_2);
        
        % Quadrant 2 Mask
%         Quadrants.BL = createMask(h2);
        handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,2} = quadCoords.BL;
        
        % Quadrant 3
        
        % Quadrant 3 will be the BOTTOM RIGHT QUADRANT
        
        q3_coord1 = [(B{1,1}(rightCoord,2));(B{1,1}(rightCoord,1))];
        q3_coord2 = [pixelInfo.Centroid(1) ; pixelInfo.Centroid(2)];
        q3_coord3 = [(B{1,1}(botCoord,2));(B{1,1}(botCoord,1))];
        
        quadrant_3_Xcoords = [q3_coord1(1);q3_coord2(1);q3_coord3(1)];
        quadrant_3_Ycoords = [q3_coord1(2);q3_coord2(2);q3_coord3(2)];
        
        % Find vertices
        q3_test = ((yCor > quadrant_3_Ycoords(1) & yCor > quadrant_3_Ycoords(2)) &...
            (xCor > quadrant_3_Xcoords(2) & xCor > quadrant_3_Xcoords(3)));
        
        Quad3_coord_indices = find(q3_test == 1);
        Quad3Flip = flipud(Quad3_coord_indices);
        % Insert found vertices into quadrant coordinates
        
        quad3Xout = quadrant_3_Xcoords;
        quad3Yout = quadrant_3_Ycoords;
        for q1i2 = 1:length(Quad3Flip);
            addXquad3 = xCor(Quad3Flip(q1i2));
            addYquad3 = yCor(Quad3Flip(q1i2));
            quad3Xout = [quad3Xout ; addXquad3];
            quad3Yout = [quad3Yout ; addYquad3];
        end
        
        quad3Xfinal = quad3Xout;
        quad3Yfinal = quad3Yout;
        
        quadCoords.BR = [quad3Xfinal quad3Yfinal];
        
%         h3 = impoly(gca, quad_3);
        
        % Quadrant 3 Mask
%         Quadrants.BR = createMask(h3);
        handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,3} = quadCoords.BR;
        % Quadrant 4
        
        % Quadrant 4 will be the TOP RIGHT QUADRANT
        
        q4_coord1 = [(B{1,1}(rightCoord,2));(B{1,1}(rightCoord,1))];
        q4_coord2 = [(B{1,1}(topCoord,2));(B{1,1}(topCoord,1))];
        q4_coord3 = [pixelInfo.Centroid(1) ; pixelInfo.Centroid(2)];
        
        quadrant_4_Xcoords = [q4_coord1(1);q4_coord2(1);q4_coord3(1)];
        quadrant_4_Ycoords = [q4_coord1(2);q4_coord2(2);q4_coord3(2)];
        
        q4_test = ((yCor < quadrant_4_Ycoords(1) & yCor < quadrant_4_Ycoords(3)) &...
            (xCor > quadrant_4_Xcoords(2) & xCor > quadrant_4_Xcoords(3)));
        
        Quad4_coord_indices = find(q4_test == 1);
        Quad4Flip = flipud(Quad4_coord_indices);
        % Insert found vertices into quadrant coordinates
        quad4Xtemp = quadrant_4_Xcoords(1);
        quad4Ytemp = quadrant_4_Ycoords(1);
        quad4Xout = quad4Xtemp;
        quad4Yout = quad4Ytemp;
        for q1i2 = 1:length(Quad4Flip);
            addXquad4 = xCor(Quad4Flip(q1i2));
            addYquad4 = yCor(Quad4Flip(q1i2));
            quad4Xout = [quad4Xout ; addXquad4];
            quad4Yout = [quad4Yout ; addYquad4];
        end

        quad4Xfinal = [quad4Xout ; quadrant_4_Xcoords(2:3)];
        quad4Yfinal = [quad4Yout ; quadrant_4_Ycoords(2:3)];
        
        quadCoords.BL = [quad4Xfinal quad4Yfinal];
        
%         h4 = impoly(gca, quad_4);
        
        % Quadrant 4 Mask
%         Quadrants.BL = createMask(h4);
        handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,4} = quadCoords.BL;

        % CALCULATE NUM PIXELS ABOVE THRESHHOLD IN EACH QUADRANT AND DIVIDE
        % BY TOTAL
        
        % injArea(i) AREA of whole polygon that exceeds threshold
        
        injImage(~mNtb_mask) = 0; % use inverse of mask to get rid of out poly pixels
        threshExceedIndex = injImage > pixThresh;
        % to plot image of injection threshold match
        injImage(~threshExceedIndex) = 0;
        
        injArea(i) = bwarea(threshExceedIndex); % BASELINE 
        
        quadInfo = struct;
        pixelsPerquad = cell(1,3);
        hold on
        quadPlotsS = cell(1,3);
        for sti = 1:4
            samplequad = roipoly(dim1,dim2,handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,sti}(:,2),...
                handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,sti}(:,1));
            quad_mask = poly2mask(handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,sti}(:,1),...
                handles.quadInfo.(strcat('sheet',num2str(handles.SheetCount))){1,sti}(:,2),dim1,dim2);
            pixelsPerquad{1,sti} = injImage(quad_mask);
            [Qi, ~] = bwboundaries(samplequad,'noholes');
            quadIndices = cell2mat(Qi);
            quadPlotsS{1,sti} = quadIndices;
            
            quadInfo.(strcat('quad',num2str(sti))).quadMask = injImage;
            quadInfo.(strcat('quad',num2str(sti))).quadMask(~quad_mask) = 0;
            quadInfo.(strcat('quad',num2str(sti))).threshExceedIndex = quadInfo.(strcat('quad',num2str(sti))).quadMask > pixThresh;
            quadInfo.(strcat('quad',num2str(sti))).quadMask(quadInfo.(strcat('quad',num2str(sti))).threshExceedIndex) = 0;
            quadInfo.(strcat('quad',num2str(sti))).quadArea = bwarea(quadInfo.(strcat('quad',num2str(sti))).threshExceedIndex);
            quadInfo.(strcat('quad',num2str(sti))).quadRatio = quadInfo.(strcat('quad',num2str(sti))).quadArea / injArea(i);
            
            plot(quadPlotsS{1,sti}(:,1),quadPlotsS{1,sti}(:,2),'y');
        end
        
        Quadrants{1,i} = quadInfo;
        
    end
    % To ensure that Standard Deviation does not artifically
    % inflate pixel threshold in excess of 255
    if pixThresh >= 255
        pixThresh = pixThresh - (2*handles.pStd);
        
        % Double check of pixel threshold 2/23/2014
        while pixThresh >= 255
            pixThresh = pixThresh - (handles.pStd/2);
        end
        
    end
    
    set(handles.meanT,'Visible','on')
    set(handles.sdT,'Visible','on')
    set(handles.thresT,'Visible','on')
    set(handles.plusSign,'Visible','on')
    set(handles.prodSign,'Visible','on')
    set(handles.equalSign,'Visible','on')
    set(handles.stdT,'Visible','on')
    
    set(handles.meanT,'String',num2str(handles.pMean));
    set(handles.sdT,'String',num2str(handles.pStd));
    set(handles.thresT,'String',num2str(pixThresh));
        
    if boxToggle
        % Get box outside polygon
        maxPpoint = min(Ycoords{i,1});
        minPpoint = max(Ycoords{i,1});
        leftPpoint = min(Xcoords{i,1});
        rightPpoint = max(Xcoords{i,1});
        midHorz = ((rightPpoint - leftPpoint)/2) + leftPpoint;
        midVert = ((maxPpoint - minPpoint)/2) + minPpoint;
        
        leftBpoint =  leftPpoint - ((midHorz - leftPpoint)/4);
        rightBpoint = rightPpoint + ((rightPpoint - midHorz)/4);
        topBpoint = maxPpoint - ((midVert - maxPpoint)/4);
        bottomBpoint = minPpoint + ((minPpoint - midVert)/4);
        
        BoxXcoords = round([leftBpoint , rightBpoint, rightBpoint, leftBpoint]);
        BoxYcoords = round([topBpoint , topBpoint, bottomBpoint, bottomBpoint]);
        
        Out_square_mask = poly2mask(BoxXcoords,BoxYcoords,dim1,dim2);
        
        outerbox = roipoly(dim1,dim2,BoxYcoords,BoxXcoords);
        
        [oBi, ~] = bwboundaries(outerbox,'noholes');
        OutboxIndices = cell2mat(oBi);
        
        % To get outer ring mask
        outerRingMask = Out_square_mask & ~mNtb_mask;
        
        modImage = injImage; % save image in new file
        modImage(~mNtb_mask) = 0; % use inverse of mask to get rid of out poly pixels
        threshExceedIndex = modImage > pixThresh;
        
        % Outside injection area
        
        modImage2 = injImage;
        modImage2(~outerRingMask) = 0;
        OthreshExceedIndex = modImage2 > pixThresh;
        
        % to plot image of injection threshold match
        modImage(~threshExceedIndex) = 0;
        modImage2(~OthreshExceedIndex) = 0;
        
        % calculate area
        injArea(i) = bwarea(threshExceedIndex);
        
        % ADD CODE to INCREASE BOX PROGRAMMTICALLY
        
        minRow = min(OutboxIndices(:,2));
        maxRow = max(OutboxIndices(:,2));
        minCol = min(OutboxIndices(:,1));
        maxCol = max(OutboxIndices(:,1));
        
        line1 = modImage2(minRow,minCol:maxCol);
        line2 = modImage2(minRow:maxRow,maxCol);
        line3 = modImage2(maxRow,minCol:maxCol);
        line4 = modImage2(minRow:maxRow,minCol);
        
        newbox = 0;
        fline = struct;
        for line = 1:4
            switch line
                case 1
                    percBor = sum(line1 > pixThresh)/numel(line1);
                    
                    nMinRow1 = minRow;
                    while percBor > 0.1 && nMinRow1 ~= dim2;
                        nMinRow1 = nMinRow1 - 1;
                        nline1 = injImage(nMinRow1,minCol:maxCol);
                        
                        percBor = sum(nline1 > pixThresh)/numel(nline1);
                        newbox = 1;
                    end
                    
                    fline.line1.y = nMinRow1;
                    fline.line2.y = nMinRow1;
                    
                case 2
                    percBor = sum(line2 > pixThresh)/numel(line2);
                    
                    nMaxCol2 = maxCol;
                    while percBor > 0.1 && nMaxCol2 ~= dim1;
                        nMaxCol2 = nMaxCol2 + 1;
                        nline2 = injImage(minRow:maxRow,nMaxCol2);
                        
                        percBor = sum(nline2 > pixThresh)/numel(nline2);
                        newbox = 1;
                    end
                    
                    fline.line2.x = nMaxCol2;
                    fline.line3.x = nMaxCol2;
                    
                case 3
                    percBor = sum(line3 > pixThresh)/numel(line3);
                    
                    nMaxRow3 = maxRow;
                    while percBor > 0.1 && nMaxRow3 ~= dim2;
                        nMaxRow3 = nMaxRow3 + 1;
                        nline3 = injImage(nMaxRow3,minCol:maxCol);
                        
                        percBor = sum(nline3 > pixThresh)/numel(nline3);
                        newbox = 1;
                    end
                    
                    fline.line3.y = nMaxRow3;
                    fline.line4.y = nMaxRow3;
                    
                    
                case 4
                    percBor = sum(line4 > pixThresh)/numel(line4);
                    
                    nMinCol4 = minCol;
                    while percBor > 0.1 && nMinCol4 ~= dim1;
                        nMinCol4 = nMinCol4 - 1;
                        nline4 = injImage(minRow:maxRow,nMinCol4);
                        
                        percBor = sum(nline4 > pixThresh)/numel(nline4);
                        newbox = 1;
                    end
                    
                    fline.line1.x = nMinCol4;
                    fline.line4.x = nMinCol4;
                    
            end
        end

        if newbox
            nBoxXc = zeros(1,4);
            nBoxYc = zeros(1,4);
            for nb = 1:4
                nBoxXc(nb) = fline.(strcat('line',num2str(nb))).x;
                nBoxYc(nb) = fline.(strcat('line',num2str(nb))).y;
            end
            
            %RECALCULATE NEW BOX
            
            NEW_Out_square_mask = poly2mask(nBoxXc,nBoxYc,dim1,dim2);
            
            NEW_outerbox = roipoly(dim1,dim2,nBoxYc,nBoxXc);
            
            [NoBi, ~] = bwboundaries(NEW_outerbox,'noholes');
            NewOutboxIndices = cell2mat(NoBi);
            
            % To get outer ring mask
            NEWouterRingMask = NEW_Out_square_mask & ~mNtb_mask;
            
            % Outside injection area
            modImage3 = injImage;
            modImage3(~NEWouterRingMask) = 0;
            NEWOthreshExceedIndex = modImage3 > pixThresh;
            
            % to plot image of injection threshold match
            modImage(~threshExceedIndex) = 0;
            modImage3(~NEWOthreshExceedIndex) = 0;
            
            % calculate out injection area
            oinjArea(i) = bwarea(NEWOthreshExceedIndex);
        end
        
        blankImage2 = uint8(zeros(dim1,dim2,3));
        blankImage2(:,:,1) = modImage;
        
        if newbox
            blankImage2(:,:,3) = modImage3;
        else
            blankImage2(:,:,3) = modImage2;
        end
        
        cla(handles.imDisplay)
        imshow(blankImage2);
        hold on
        plot(Xcoords{i,1}, Ycoords{i,1},'-y');
        
        if newbox
            plot(NewOutboxIndices(:,1),NewOutboxIndices(:,2),'r')
        else
            plot(OutboxIndices(:,1),OutboxIndices(:,2),'r')
        end
        
        hold on
        for ip = 1:3
            plot(boxPlotsS{1,ip}(:,1),boxPlotsS{1,ip}(:,2),'y')
        end
        
        set(handles.infoT,'String','Press enter for next image');
        handles.figActive = 1;
        guidata(hObject, handles);
        
        % Save file image or not
        dataFrame = getframe(handles.imDisplay); % LINE BROKEN
        [handles.im2save,~] = frame2im(dataFrame);
        guidata(hObject, handles);
        
        pause
        
        set(handles.infoT,'String',[]);

        % Area of polygon
        dataTOout.polyArea(i,1) = round(polyArea(i));
        % Area of injection in polygon
        dataTOout.injArea(i,1) = round(injArea(i));
        % Ratio of injection in polygon
        dataTOout.polyInjRatio(i,1) = round(((injArea(i)/polyArea(i))*1000))/1000;
        % Area of injection outside polygon
        dataTOout.outinjArea(i,1) = round(oinjArea(i));
        % Ratio of injection outstide polygon
        tempOArea = oinjArea(i,1)/(oinjArea(i)+injArea(i));
        dataTOout.injOutRatio(i,1) = round((tempOArea*1000))/1000;
        
        
    else

        
%         blankImage2 = uint8(zeros(dim1,dim2,3));
%         blankImage2(:,:,1) = modImage;
%         
%         cla(handles.imDisplay)
%         imshow(blankImage2);
        hold on
        plot(Xcoords{i,1}, Ycoords{i,1},'-r');
        
        set(handles.infoT,'String','Press enter for next image');
        handles.figActive = 1;
        guidata(hObject, handles);
        
        % Save file image or not
        dataFrame = getframe(handles.imDisplay); % LINE BROKEN
        [handles.im2save,~] = frame2im(dataFrame);
        guidata(hObject, handles);
        
        pause
        
        set(handles.infoT,'String',[]);
        
        % Area of polygon
        dataTOout.polyArea(i,1) = round(polyArea(i));
        % Area of injection in polygon
        dataTOout.injArea(i,1) = round(injArea(i));
        % Ratio of injection in polygon
        dataTOout.polyInjRatio(i,1) = round(((injArea(i)/polyArea(i))*1000))/1000;
        % Area of injection outside polygon
        dataTOout.outinjArea(i,1) = NaN;
        % Ratio of injection outstide polygon
        dataTOout.injOutRatio(i,1) = NaN;
    end
    tempDout = [dataTOout.polyArea , dataTOout.injArea , dataTOout.polyInjRatio,...
        dataTOout.outinjArea , dataTOout.injOutRatio];
    
    set(handles.dataTable,'Data',num2cell(tempDout))
    
    % INCORPORATE ADDITIOINAL SHEETS
end

handles.AllDATA.(strcat('sheet',num2str(handles.SheetCount))) = get(handles.dataTable,'Data');
handles.AllDATA2.(strcat('sheet',num2str(handles.SheetCount))) = Quadrants;
handles.CaseNames{handles.SheetCount} = caseName;
 	
set(handles.dataTable,'Data','');
handles.SheetCount = handles.SheetCount + 1;
% DETERMINE IF USER WANTS TO ANALYZE ANOTHER DATA SET
msgbox('To Start another image set select LOAD FOLDER', 'INFO','help');
 	
cla(handles.imDisplay)
% if handles.hemiCount == 1;
%     cla(handles.imDisplay)
%     set(handles.infoT,'String','Select opposite Hemisphere');
% elseif handles.hemiCount == 2;
%     set(handles.infoT,'String','Analysis Complete!');
%     cla(handles.imDisplay)
%     set(handles.gThresh,'Enable','off')
%     set(handles.roiButton,'Enable','off')
%     set(handles.hemiPanel,'Visible','off')
%     set(handles.fileOpts,'Enable','on');
%     set(handles.expXl,'Enable','on');
%     
% end





guidata(hObject, handles);

% --- Executes on button press in gThresh.
function gThresh_Callback(hObject, ~, handles)
% hObject    handle to gThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% CHANGE TEXT MESSAGE
set(handles.infoT,'String','Select Draw ROI');
% ENABLE ROI DRAWING BUTTON
set(handles.roiButton,'Enable','on');
% DISABLE THRESHOLD BUTTON
set(handles.gThresh,'Enable','off');
% DISABLE SECTIONLIST
set(handles.sectList,'Enable','off')

guidata(hObject, handles);




function stdT_Callback(hObject, eventdata, handles)
% hObject    handle to stdT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stdT as text
%        str2double(get(hObject,'String')) returns contents of stdT as a double

% std2use = str2double(get(hObject,'String'));
%
% handles.pThresh = handles.pMean + (std2use*handles.pStd);
%
% set(handles.thresT,'String',num2str(handles.pThresh));
%
% guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function stdT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function imOpts_Callback(hObject, eventdata, handles)
% hObject    handle to imOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function expImg_Callback(hObject, eventdata, handles)
% hObject    handle to expImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.figActive == 1;
    handles.saveLoc = uigetdir('C:\');
    cd(handles.saveLoc)
    
    caseName = inputdlg('Save file ID','FILE name',[1 30],{'Image'});
    
    saveName = strcat('AnnalysisImage_',date,'_',char(caseName),'.tiff');
    
    imwrite(handles.im2save,saveName);
end


% --------------------------------------------------------------------
function clrel_Callback(hObject, eventdata, handles)
% hObject    handle to clrel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exitquest = questdlg('Are you sure you want to RESTART the program?','Restart?',...
    'Yes','No','Yes');

switch exitquest
    case 'Yes'
        delete(handles.figure1);
        DyeDist_dev2;
    case 'No'
        return
end


% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exitquest = questdlg('Are you sure you want to EXIT the program?','Exit?',...
    'Yes','No','Yes');

switch exitquest
    case 'Yes'
        delete(handles.figure1);
    case 'No'
        return
end

% --------------------------------------------------------------------
function expOpts_Callback(hObject, eventdata, handles)
% hObject    handle to expOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function expXl_Callback(hObject, eventdata, handles)
% hObject    handle to expXl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

expquest = questdlg('Are you sure you want to EXPORT TO EXCEL?','EXCEL?',...
    'Yes','No','Yes');

switch expquest
    case 'Yes'
        handles.saveLoc = uigetdir('C:\');
        cd(handles.saveLoc)
        
        versionCheck = version('-release');
        getYear = str2double(versionCheck(1:end-1));
        if getYear < 2013
            DS_toggle = 0;
        else
            DS_toggle = 1;
        end
        
        numofSheets = numel(fieldnames(handles.AllDATA));
        
        orientCh = questdlg('Sagittal or Coronal sections?','Orientation','Sagittal','Coronal','Sagittal');
        
        switch orientCh
            case 'Sagittal'
                orientKey = [{'TL';'BL';'BR';'TR'} , {'AD';'AV';'PV';'PD'}];
            case 'Coronal'
                orientKey = [{'TL';'BL';'BR';'TR'} , {'LD';'LV';'MV';'MD'}];
        end
        
        % For generic data
        for nsi = 1:numofSheets
            
            tempCols = handles.ExportColnames{nsi};
            tempSheet = handles.CaseNames{nsi};
            tempData = handles.AllDATA.(strcat('sheet',num2str(nsi)));
            
            % clean up data
            if isnan(tempData{1,4})
                tempData = tempData(:,1:3);
            end
            
            if DS_toggle
                outData = cell2dataset(tempData,'VarNames',tempCols);
            else
                % Do something
            end
            saveName = char(strcat('ImData_',date,'_',tempSheet,'.xlsx'));
            export(outData,'XLSfile',saveName)
        end
        
        % For quadrant analysis % ADD EXPORT STUFF
        for nsi = 1:numofSheets
            
            tempCols = handles.ExportColnames{nsi};
            tempSheet = handles.CaseNames{nsi};
            tempData = handles.AllDATA.(strcat('sheet',num2str(nsi)));
            
            % clean up data
            if isnan(tempData{1,4})
                tempData = tempData(:,1:3);
            end
            
            if DS_toggle
                outData = cell2dataset(tempData,'VarNames',tempCols);
            else
                % Do something
            end
            saveName = char(strcat('ImData_',date,'_',tempSheet,'.xlsx'));
            export(outData,'XLSfile',saveName)
        end
        
    case 'No'
        return
end

% --------------------------------------------------------------------
function expML_Callback(hObject, eventdata, handles)
% hObject    handle to expML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function expSTRUCT_Callback(hObject, eventdata, handles)
% hObject    handle to expSTRUCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

expquest = questdlg('Are you sure you want to EXPORT TO STRUCT?','STRUCT?',...
    'Yes','No','Yes');

switch expquest
    case 'Yes'
        handles.saveLoc = uigetdir('C:\');
        cd(handles.saveLoc)
        
        numofSheets = numel(fieldnames(handles.AllDATA));
        OutData = struct;
        for nsi = 1:numofSheets
            
            tempCols = handles.ExportColnames{nsi};
            tempSheet = handles.CaseNames{nsi};
            tempData = handles.AllDATA.(strcat('sheet',num2str(nsi)));
            
            % clean up data
            if isnan(tempData{1,4})
                tempData = tempData(:,1:3);
            end
            
            
            for di = 1:length(tempCols)
                OutData.(tempSheet{1}).(tempCols{di}) = tempData(:,di);
            end
        end
        
        saveName = char(strcat('ImDataST_',date,'.mat'));
        save(saveName,'OutData');
        
    case 'No'
        return
end


% --------------------------------------------------------------------
function expDATASET_Callback(hObject, eventdata, handles)
% hObject    handle to expDATASET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

expquest = questdlg('Are you sure you want to EXPORT TO DATASET?','DATASET?',...
    'Yes','No','Yes');

switch expquest
    case 'Yes'
        handles.saveLoc = uigetdir('C:\');
        cd(handles.saveLoc)
        
        numofSheets = numel(fieldnames(handles.AllDATA));
        
        outDS = dataset;
        
        for nsi = 1:numofSheets
            
            tempCols = handles.ExportColnames{nsi};
            tempSheet = handles.CaseNames{nsi};
            tempData = handles.AllDATA.(strcat('sheet',num2str(nsi)));
            tempSections = (1:1:size(tempData,1))';
            tempCases = repmat(handles.CaseNames{nsi},numel(tempSections),1);
            
            % clean up data
            if isnan(tempData{1,4})
                tempData = tempData(:,1:3);
            end

            tempCols = ['Cases' , 'Sections' , tempCols];
            tempData = [num2cell(tempCases) , num2cell(tempSections) , tempData];
            
            tempDS = cell2dataset(tempData,'VarNames',tempCols);
            outDS = [outDS ; tempDS];

        end        
            saveName = char(strcat('ImDataDS_',date,'.mat'));
            save(saveName,'outDS');
        
    case 'No'
        return
end
