function varargout = DyeDist(varargin)
% DYEDIST MATLAB code for DyeDist.fig
%      DYEDIST, by itself, creates a new DYEDIST or raises the existing
%      singleton*.
%
%      H = DYEDIST returns the handle to a new DYEDIST or the handle to
%      the existing singleton*.
%
%      DYEDIST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYEDIST.M with the given input arguments.
%
%      DYEDIST('Property','Value',...) creates a new DYEDIST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DyeDist_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DyeDist_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DyeDist

% Last Modified by GUIDE v2.5 27-Feb-2014 15:03:14


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THINGS TO DO
% Add message for DRAW POLYGON
% Add export to Excel



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DyeDist_OpeningFcn, ...
    'gui_OutputFcn',  @DyeDist_OutputFcn, ...
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

% --- Executes just before DyeDist is made visible.
function DyeDist_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DyeDist (see VARARGIN)

% Choose default command line output for DyeDist
handles.output = hObject;

set(handles.sectList,'Enable','off')
set(handles.hemiPanel,'Visible','off')
set(handles.imDisplay,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[]);

set(handles.gThresh,'Enable','off');
set(handles.roiButton,'Enable','off');

set(handles.meanT,'Visible','off')
set(handles.sdT,'Visible','off')
set(handles.thresT,'Visible','off')
set(handles.plusSign,'Visible','off')
set(handles.prodSign,'Visible','off')
set(handles.equalSign,'Visible','off')
set(handles.stdT,'Visible','off')

set(handles.expXl,'Enable','off')

handles.figActive = 0;
handles.hemiCount = 0;

set(handles.infoT,'String','Load Image Files');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DyeDist wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DyeDist_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in leftB.
function leftB_Callback(hObject, ~, handles)
% hObject    handle to leftB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of leftB

set(handles.infoT,'String',[]);
set(handles.gThresh,'Enable','on');
set(handles.sectList,'Enable','on')
cla(handles.imDisplay);
lcheck = get(hObject,'Value');

if lcheck
    handles.HemiS = 'L';
    set(handles.rightB,'Value',0)
else
    handles.HemiS = 'R';
    set(handles.rightB,'Value',1)
end

set(handles.infoT,'String','View sections or Select Start Analysis');

guidata(hObject, handles);



% --- Executes on button press in rightB.
function rightB_Callback(hObject, ~, handles)
% hObject    handle to rightB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rightB

set(handles.infoT,'String',[]);
set(handles.gThresh,'Enable','on');
set(handles.sectList,'Enable','on')
cla(handles.imDisplay);
rcheck = get(hObject,'Value');

if rcheck
    handles.HemiS = 'R';
    set(handles.leftB,'Value',0)
else
    handles.HemiS = 'L';
    set(handles.leftB,'Value',1)
end

set(handles.infoT,'String','View sections or Select Start Analysis');

guidata(hObject, handles);

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

sectNumExt = cellfun(@(x) str2double(regexp(x,'\d*','match')),handles.ImgNames);

[handles.SectionIndex,sortOrder] = sort(sectNumExt);

handles.numSections = unique(handles.SectionIndex);

handles.ImgNames = handles.ImgNames(sortOrder);

handles.hemisIndex = cellfun(@(x) regexp(x,'L|R','match'),handles.ImgNames);

handles.ListIndex = arrayfun(@(x) strcat('Section_',num2str(x)),handles.numSections, 'UniformOutput', false);

set(handles.hemiPanel,'Visible','on')
set(handles.sectList,'String',handles.ListIndex);

handles.channels = inputdlg('How many channels?','Channels',[1 30],{'1'});

switch char(handles.channels)
    case '1'
        options.Interpreter = 'tex';
        options.Default = 'R';
        chan1 = questdlg('Which channel?','Channels','R','G','B',options);
        switch chan1
            case 'R'
                handles.colorId = 1;
            case 'G'
                handles.colorId = 2;
            case 'B'
                handles.colorId = 3;
        end
    case '2'
        chan2 = inputdlg({'Flur1_Name','Flur2_Name'}, 'Channel ID', [1 40; 1 40],{'R, Injection','B, dapi'});
        
        flurS = [strsplit(chan2{1},','); strsplit(chan2{2},',')];
        for fi2 = 1:2
            switch flurS{fi2,1}
                case 'R'
                    handles.colorId(fi2) = 1;
                    handles.colorIdtag{fi2} = strtrim(flurS{fi2,2});
                case 'G'
                    handles.colorId(fi2) = 2;
                    handles.colorIdtag{fi2} = strtrim(flurS{fi2,2});
                case 'B'
                    handles.colorId(fi2) = 3;
                    handles.colorIdtag{fi2} = strtrim(flurS{fi2,2});
            end
        end
    case '3'
        chan3 = inputdlg({'R','G','B'}, 'Channel ID', [1 40; 1 40; 1 40],{'Injection', 'dapi', 'GFP'});
        handles.colorIdtag = chan3;
    otherwise
        warndlg('TOO many channels: input 1-3');
        return
        
end

set(handles.infoT,'String','Choose Hemisphere');



guidata(hObject, handles);



% --- Executes on selection change in sectList.
function sectList_Callback(hObject, ~, handles)
% hObject    handle to sectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sectList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sectList

secChoice = get(hObject,'Value');
handles.section2use = handles.numSections(secChoice);
secIndex = handles.SectionIndex == handles.section2use;

sec2useIndex = (secIndex == 1) & (strcmp(handles.HemiS,handles.hemisIndex));

if sum(sec2useIndex) == 0;
    
    if strcmp(handles.HemiS,'L')
        warnSide = 'Left';
    else
        warnSide = 'Right';
    end
    
    warnMessage =  sprintf('No %s image for this section',warnSide);
    
    warndlg(warnMessage,'Hemisphere Warning');
    
    return
    
else
    
    sections = handles.ImgNames(sec2useIndex);
    
    imForsize = imread(sections{1});
    
    [dim1,dim2] = size(imForsize);
    
    blankImage = uint8(zeros(dim1,dim2,3));
    switch char(handles.channels)
        case '1'
            for i = 1:3
                if i == handles.colorId;
                    section2incl = imread(sections{1});
                    blankImage(:,:,i) = section2incl;
                else
                    continue
                end
            end
            
        case '2'
            if isempty(strfind(sections{1},handles.colorIdtag{1}))
                % then section 1 is associated with colorId tag 2
                if length(sections) == 2
                    section1 = imread(sections{1});
                    blankImage(:,:,handles.colorId(2)) = section1;
                    
                    section2 = imread(sections{2});
                    blankImage(:,:,handles.colorId(1)) = section2;
                else
                    section1 = imread(sections{1});
                    blankImage(:,:,handles.colorId(2)) = section1;
                end
                
            else
                if length(sections) == 2
                    section1 = imread(sections{1});
                    blankImage(:,:,handles.colorId(1)) = section1;
                    
                    section2 = imread(sections{2});
                    blankImage(:,:,handles.colorId(2)) = section2;
                else
                    section1 = imread(sections{1});
                    blankImage(:,:,handles.colorId(1)) = section1;
                end
            end
            
        case '3'
            
    end
    
    handles.image2show = blankImage;
    
    cla(handles.imDisplay);
    axes(handles.imDisplay);
    imshow(handles.image2show);
    
end


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

set(handles.fileOpts,'Enable','off');

set(handles.infoT,'String',[]);
cla(handles.imDisplay);
set(handles.sectList,'Enable','off');

set(handles.gThresh,'Enable','off');
set(handles.stdT,'Enable','off');

conditionQ = questdlg('Identify ROI tag','Condition',handles.colorIdtag{1},handles.colorIdtag{2},handles.colorIdtag{2});

Xcoords = cell(length(handles.numSections),1);
Ycoords = cell(length(handles.numSections),1);
polyArea = zeros(length(handles.numSections),1);
injArea = zeros(length(handles.numSections),1);
oinjArea = zeros(length(handles.numSections),1);

dataTOout = struct;

axes(handles.imDisplay)

% Check column titles
curTableN = get(handles.dataTable,'ColumnName');
if isempty(curTableN{1})
    if strcmp(handles.HemiS,'L')
        columnNames = {'LRoi_Area','L_Inj_Area','L_Inj_Ratio','L_InjO_Area',...
            'L_InjO_Ratio'};
        set(handles.dataTable,'ColumnName',columnNames);
    elseif strcmp(handles.HemiS,'R')
        columnNames = {'RRoi_Area','R_Inj_Area','R_Inj_Ratio','R_InjO_Area',...
            'R_InjO_Ratio'};
        set(handles.dataTable,'ColumnName',columnNames);
    end
else
    if strcmp(handles.HemiS,'L') && ~strcmp(curTableN{1}(1),'L')
        columnNames = [curTableN' , {'LRoi_Area','L_Inj_Area','L_Inj_Ratio','L_InjO_Area',...
            'L_InjO_Ratio'}];
        set(handles.dataTable,'ColumnName',columnNames);
    elseif strcmp(handles.HemiS,'R') && ~strcmp(curTableN{1}(1),'R')
        columnNames = [curTableN', {'RRoi_Area','R_Inj_Area','R_Inj_Ratio','R_InjO_Area',...
            'R_InjO_Ratio'}];
        set(handles.dataTable,'ColumnName',columnNames);
    end
end

% Check table data
currentTdata = get(handles.dataTable,'Data');

if strcmp(currentTdata{1},'');
    previousTable = 0;
else
    [PreRows,~] = size(currentTdata);
    preTable = currentTdata;
    previousTable = 1;
end


for i = 1:length(handles.numSections)
    
    cd(handles.WD)
    
    handles.section2use = handles.numSections(i);
    secIndex = handles.SectionIndex == handles.section2use;
    sec2useIndex = (secIndex == 1) & (strcmp(handles.HemiS,handles.hemisIndex));
    sections = handles.ImgNames(sec2useIndex);
    
    % get condition section and put in green channel
    
    conIndex = ~isempty(cellfun(@(x) strfind(x,conditionQ), sections, 'UniformOutput',false));
    numSecs = 1:2;
    
    if isempty(sections)
        continue
        
    else
        
        if length(sections) > 1;
            oppIndex = numSecs(~ismember(1:2,conIndex));
            injection = sections{oppIndex};
            injImage = imread(injection);
        else
            injArea(i) = NaN;
            oinjArea(i) = NaN;
            dataTOout.datapolyInjRatio(i) = NaN;
            dataTOout.outinjArea(i) = NaN;
            dataTOout.injOutRatio(i) = NaN;
        end
        
        section2plot = sections{conIndex};
        sec2draw = imread(section2plot);
        [dim1,dim2] = size(sec2draw);
        % Create 3 images using blank image background
        
        blankImage = uint8(zeros(dim1,dim2,3));
        % Trace Image
        traceImage = blankImage;
        traceImage(:,:,2) = sec2draw;
        
        cla(handles.imDisplay)
        imshow(traceImage);
        
        set(handles.infoT,'String','Trace Region of Interest');
        [~, Xcoords{i,1}, Ycoords{i,1}] = roipoly(traceImage);
        
        polyArea(i) = polyarea(Xcoords{i,1},Ycoords{i,1});
        
        % Create NTS mask
        mNtb_mask = poly2mask(Xcoords{i,1},Ycoords{i,1},dim1,dim2);
        
        if length(sections) > 1;
            
            % Create injection image for thresholding
            
            injThreshImage = blankImage;
            injThreshImage(:,:,1) = injImage;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get sample boxes
            
            totPixels = dim1*dim2;
            onePercent = totPixels*0.005;
            boxDim = round(sqrt(onePercent));
            halfBD = round(boxDim/2);
            
            cla(handles.imDisplay);
            axes(handles.imDisplay);
            imshow(injThreshImage);
            
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
            
            if pixThresh < 15
                pixTtoggle = 0;
            else
                pixTtoggle = 1;
            end
            
            if pixTtoggle
                
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
            end
            
            if pixTtoggle
                
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
                dataFrame = getframe(handles.imDisplay);
                [handles.im2save,~] = frame2im(dataFrame);
                guidata(hObject, handles);
                
                pause
                
                set(handles.infoT,'String',[]);
                
            end
            
            if pixTtoggle
                
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
                
                tempDout = [dataTOout.polyArea , dataTOout.injArea , dataTOout.polyInjRatio,...
                    dataTOout.outinjArea , dataTOout.injOutRatio];
                
            else
                
                dataTOout.polyArea(i,1) = round(polyArea(i));
                dataTOout.injArea(i,1) = NaN;
                dataTOout.polyInjRatio(i,1) = NaN;
                dataTOout.outinjArea(i,1) = NaN;
                dataTOout.injOutRatio(i,1) = NaN;
                
                tempDout = [dataTOout.polyArea , dataTOout.injArea , dataTOout.polyInjRatio,...
                    dataTOout.outinjArea , dataTOout.injOutRatio];
                
            end
            
            
            if ~previousTable
                
                polyAreaDT = [dataTOout.polyArea , dataTOout.injArea , dataTOout.polyInjRatio,...
                    dataTOout.outinjArea , dataTOout.injOutRatio];
                polyAreaDT = num2cell(polyAreaDT);
                
            else
                % find out length of previous table add blanks to the lesser one
                % and bind together
                
                tempAreaDT = [dataTOout.polyArea , dataTOout.injArea , dataTOout.polyInjRatio,...
                    dataTOout.outinjArea , dataTOout.injOutRatio];
                tempAreaDT = num2cell(tempAreaDT);
                
                [CurRows ,~] = size(tempDout);
                
                if CurRows < PreRows % current rows less than previous rows
                    diffRow = PreRows - CurRows;
                    filler = repmat({''},diffRow,5);
                    combineOut = [tempAreaDT ; filler];
                else
                    diffRow = CurRows - PreRows;
                    filler = repmat({''},diffRow,5);
                    combineOut = [tempAreaDT ; filler];
                end
                
                polyAreaDT =  [preTable , combineOut];
                
            end
            
            set(handles.dataTable,'Data',polyAreaDT);
            
        else
            polyAreaDT = [polyArea , injArea];
            set(handles.dataTable,'Data',polyAreaDT);
            
        end
        
    end
    
end


if strcmp(handles.HemiS,'L')
    set(handles.rightB,'Enable','on');
elseif strcmp(handles.HemiS,'R')
    set(handles.leftB,'Enable','on');
end

handles.hemiCount = handles.hemiCount + 1;

if handles.hemiCount == 1;
    cla(handles.imDisplay)
    set(handles.infoT,'String','Select opposite Hemisphere');
elseif handles.hemiCount == 2;
    set(handles.infoT,'String','Analysis Complete!');
    cla(handles.imDisplay)
    set(handles.gThresh,'Enable','off')
    set(handles.roiButton,'Enable','off')
    set(handles.hemiPanel,'Visible','off')
    set(handles.fileOpts,'Enable','on');
    set(handles.expXl,'Enable','on');
    
end



guidata(hObject, handles);

% --- Executes on button press in gThresh.
function gThresh_Callback(hObject, ~, handles)
% hObject    handle to gThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% secChoice = get(handles.sectList,'Value');
%
% handles.section2use = handles.numSections(secChoice);
% secIndex = handles.SectionIndex == handles.section2use;
%
% sec2useIndex = (secIndex == 1) & (strcmp(handles.HemiS,handles.hemisIndex));
% sections = handles.ImgNames(sec2useIndex);
%
% imForsize = imread(sections{1});
%
% [dim1,dim2] = size(imForsize);
%
% blankImage = uint8(zeros(dim1,dim2,3));
% switch char(handles.channels)
%     case '1'
%         for i = 1:3
%             if i == handles.colorId;
%                 section2incl = imread(sections{1});
%                 blankImage(:,:,i) = section2incl;
%             else
%                 continue
%             end
%         end
%
%     case '2'
%         if isempty(strfind(sections{1},handles.colorIdtag{1}))
%             % then section 1 is associated with colorId tag 2
%             if length(sections) == 2
%                 section1 = imread(sections{1});
%                 blankImage(:,:,handles.colorId(2)) = section1;
%
%                 section2 = imread(sections{2});
%                 blankImage(:,:,handles.colorId(1)) = section2;
%             else
%                 section1 = imread(sections{1});
%                 blankImage(:,:,handles.colorId(2)) = section1;
%             end
%
%         else
%             if length(sections) == 2
%                 section1 = imread(sections{1});
%                 blankImage(:,:,handles.colorId(1)) = section1;
%
%                 section2 = imread(sections{2});
%                 blankImage(:,:,handles.colorId(2)) = section2;
%             else
%                 section1 = imread(sections{1});
%                 blankImage(:,:,handles.colorId(1)) = section1;
%             end
%         end
%
%     case '3'
%
% end
%
% handles.image2show = blankImage;
% handles.imageforThresh = blankImage(:,:,1);

% figure;
% imshow(handles.imageforThresh)

% totPixels = dim1*dim2;
% onePercent = totPixels*0.005;
% boxDim = round(sqrt(onePercent));
% halfBD = round(boxDim/2);
%
% cla(handles.imDisplay);
% axes(handles.imDisplay);
% imshow(handles.image2show);
%
% [x_coord, y_coord] = ginput(3);
% x_coord = round(x_coord);
% y_coord = round(y_coord);
%
% Ycoords = cell(1,3);
% Xcoords = cell(1,3);
% for i = 1:length(x_coord)
%     Xcoords{1,i} = round([x_coord(i) - halfBD , x_coord(i) + halfBD, x_coord(i) + halfBD, x_coord(i) - halfBD]);
%     Ycoords{1,i} = round([y_coord(i) - halfBD , y_coord(i) - halfBD, y_coord(i) + halfBD, y_coord(i) + halfBD]);
% end
%
% pixelsPerbox = cell(1,3);
% hold on
% for i = 1:3
%     samplebox = roipoly(dim1,dim2,Ycoords{1,i},Xcoords{1,i});
%
%     square_mask = poly2mask(Xcoords{1,i},Ycoords{1,i},dim1,dim2);
%
%     pixelsPerbox{1,i} = handles.imageforThresh(square_mask);
%
%     [Bi, ~] = bwboundaries(samplebox,'noholes');
%
%     boxIndices = cell2mat(Bi);
%
%     plot(boxIndices(:,1),boxIndices(:,2),'y')
%
% end
%
% handles.pixelsBackground = double([pixelsPerbox{1,1};pixelsPerbox{1,2};pixelsPerbox{1,3}]);
%
% handles.pMean = round(mean(handles.pixelsBackground));
% handles.pStd = round(std(handles.pixelsBackground));
%
% handles.pThresh = handles.pMean + (str2double(get(handles.stdT,'String'))*handles.pStd);
%
% set(handles.meanT,'Visible','on')
% set(handles.sdT,'Visible','on')
% set(handles.thresT,'Visible','on')
% set(handles.plusSign,'Visible','on')
% set(handles.prodSign,'Visible','on')
% set(handles.equalSign,'Visible','on')
% set(handles.stdT,'Visible','on')
%
% set(handles.meanT,'String',num2str(handles.pMean));
% set(handles.sdT,'String',num2str(handles.pStd));
% set(handles.thresT,'String',num2str(handles.pThresh));
set(handles.infoT,'String','Select Draw ROI');

set(handles.roiButton,'Enable','on');

set(handles.gThresh,'Enable','off');

set(handles.sectList,'Enable','off')

set(handles.leftB,'Enable','off');
set(handles.rightB,'Enable','off');

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


delete(handles.figure1);
AnnalysisGUI_v02;


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

outData = get(handles.dataTable,'Data');
outDataCol = get(handles.dataTable,'ColumnName');outDataCol = outDataCol';

outData(cellfun(@(x) isempty(x), outData)) = {nan};

convert = cell2mat(outData);

expData = dataset(convert(:,1),convert(:,2),convert(:,3),convert(:,4),...
    convert(:,5),convert(:,6),convert(:,7),convert(:,8),convert(:,9),...
    convert(:,10));
expData.Properties.VarNames = outDataCol;

caseName = inputdlg('Save file ID','FILE name',[1 30],{'data'});

saveName = strcat('AnnalysisOutput_',date,'_',char(caseName),'.xlsx');

export(expData,'XLSfile',saveName);
