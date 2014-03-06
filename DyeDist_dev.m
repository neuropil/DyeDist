function varargout = DyeDist_dev(varargin)
% DYEDIST_DEV MATLAB code for DyeDist_dev.fig
%      DYEDIST_DEV, by itself, creates a new DYEDIST_DEV or raises the existing
%      singleton*.
%
%      H = DYEDIST_DEV returns the handle to a new DYEDIST_DEV or the handle to
%      the existing singleton*.
%
%      DYEDIST_DEV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYEDIST_DEV.M with the given input arguments.
%
%      DYEDIST_DEV('Property','Value',...) creates a new DYEDIST_DEV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DyeDist_dev_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DyeDist_dev_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DyeDist_dev

% Last Modified by GUIDE v2.5 04-Mar-2014 14:57:15


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% THINGS TO DO
% ADD NEW DISTRIBUTION ANALYSIS by quadrant (3/6/2014)
% ADD ABILITY TO NAME SHEETS FOR DIFFERENT DATA SETS 
% FINIALIZE EXPORT FEATURES


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DyeDist_dev_OpeningFcn, ...
    'gui_OutputFcn',  @DyeDist_dev_OutputFcn, ...
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

% --- Executes just before DyeDist_dev is made visible.
function DyeDist_dev_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DyeDist_dev (see VARARGIN)

% Choose default command line output for DyeDist_dev
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

% UIWAIT makes DyeDist_dev wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DyeDist_dev_OutputFcn(~, ~, handles)
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
handles.ImgNames = actualNames(cellfun(@(x) strcmp('s',x(1)),actualNames));

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

axes(handles.imDisplay)

%%%%% FIGURE HOW TO KEEP TRACK OF MULTIPLE CASES/SHEETS

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
        cla(handles.imDisplay); % REPEAT THIS LINE OUTSIDE TOGGLE
        %     axes(handles.imDisplay); % REPEAT THIS LINE OUTSIDE TOGGLE
        imshow(injThreshImage); % REPEAT THIS LINE OUTSIDE TOGGLE

        % calculate threshold
        num_inj_pixels = injImage(mNtb_mask);
        convertPixels2double = single(num_inj_pixels);
        
        handles.pixelsBackground =  double(convertPixels2double);
        
        handles.pMean = round(mean(handles.pixelsBackground));
        handles.pStd = round(std(handles.pixelsBackground));
        
        pixThresh = handles.pMean + (str2double(get(handles.stdT,'String'))*handles.pStd);
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
        modImage = injImage; % save image in new file
        modImage(~mNtb_mask) = 0; % use inverse of mask to get rid of out poly pixels
        threshExceedIndex = modImage > pixThresh;
        % to plot image of injection threshold match
        modImage(~threshExceedIndex) = 0;
        
        injArea(i) = bwarea(threshExceedIndex);
        
        blankImage2 = uint8(zeros(dim1,dim2,3));
        blankImage2(:,:,1) = modImage;
        
        cla(handles.imDisplay)
        imshow(blankImage2);
        hold on
        plot(Xcoords{i,1}, Ycoords{i,1},'-y');
        
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
        DyeDist_dev;
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
caseName = inputdlg('Save file ID','FILE name',[1 30],{'data'});

for nsi = 1:numofSheets
    
    tempCols = handles.ExportColnames{nsi};
    tempSheet = handles.CaseNames{nsi};
    tempData = handles.AllDATA.(strcat('sheet',num2str(nsi)));
    
    if DS_toggle
        outData = cell2dataset(tempData,'VarNames',tempCols);
    else
        % Do something
    end
    
    export(outData,'XLSfile',saveName,'sheet',tempSheet)
end


saveName = strcat('AnnalysisOutput_',date,'_',char(caseName),'.xlsx');

export(expData,'XLSfile',saveName);

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


% --------------------------------------------------------------------
function expCELLARRAY_Callback(hObject, eventdata, handles)
% hObject    handle to expCELLARRAY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function expDATASET_Callback(hObject, eventdata, handles)
% hObject    handle to expDATASET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
