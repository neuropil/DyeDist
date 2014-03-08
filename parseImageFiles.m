function [imageID] = parseImageFiles(imFileNames,chanNum)
%PARSEIMAGEFILES Summary of this function goes here
%   Detailed explanation goes here

% Ignore mismatches; throw them out and throw out warning message 



% Check if number of inputs are correct
if nargin < 2 || nargin > 2
    error('Must input two arguments');
end
% Check if Channel number input is numeric and 3 or less
if ~isnumeric(chanNum)
    error('Second input argument must be numeric')
elseif chanNum > 3
    error('Second input cannot be higher than 3')
end

% If one channel
% Look for number and condition separated by '_'
if chanNum == 1
    
    nisslCond = '';
    roiCond = '';
    
    imIter = 1;
    nToggle = 1; rToggle = 1;
    while (nToggle || rToggle || imIter < length(imFileNames)) 
        
        str = imFileNames{imIter};
        [~,~,e] = fileparts(str);
        foundNames = strsplit(str(1:strfind(str,e)-1),'_');
        
        if ismember(nisslCond,foundNames) || ismember(roiCond,foundNames)
            imIter = imIter + 1;
           continue 
        end
        
        cond1Sl = listdlg('PromptString','Select Condition Name','SelectionMode','Single',...
            'ListString',foundNames);
        
        condSel = foundNames{cond1Sl};
        
        cQuest = sprintf('Is %s the Nissl or ROI condition?',condSel);
        
        cond1Q = questdlg(cQuest, 'Condition', 'Nissl', 'ROI', 'Nissl');
        
        switch cond1Q
            case 'Nissl'
                nisslCond = condSel;
                nToggle = 0;
            case 'ROI'
                roiCond = condSel;
                rToggle = 0;
        end

        imIter = imIter + 1;
    end
    
    % Get list of both condition files
    conTypes = {'Nissl', 'ROI'};
    conds = {nisslCond , roiCond};
    nisslFiles = cell(length(imFileNames),1);
    roiFiles = cell(length(imFileNames),1);
    for ci = 1:2
        condCheck = conds{ci};
        imCount = 1;
        for imI = 1:length(imFileNames)
            if strcmp(regexp(imFileNames{imI},condCheck,'match'),condCheck);
                switch conTypes{ci}
                    case 'Nissl'
                        nisslFiles{imCount} = imFileNames{imI};
                        imCount = imCount + 1;
                    case 'ROI'
                        roiFiles{imCount} = imFileNames{imI};
                        imCount = imCount + 1;
                end
            else
                continue
            end
        end
    end
    % Clean up
    nisslFiles = nisslFiles(cellfun(@(x) ~isempty(x), nisslFiles));
    roiFiles = roiFiles(cellfun(@(x) ~isempty(x), roiFiles));

    % Extract Numbers
    combinedFiles = {nisslFiles , roiFiles};
    nisslFnums = cell(length(imFileNames),1);
    roiFnums = cell(length(imFileNames),1);
    for conI = 1:2
        
        tempFnames = combinedFiles{conI};
        for fileI = 1:length(tempFnames)
            switch conTypes{conI}
                case 'Nissl'
                     nisslFnums{fileI} = str2double(regexp(tempFnames{fileI},'[0-9]','match'));
                case 'ROI'
                     roiFnums{fileI} = str2double(regexp(tempFnames{fileI},'[0-9]','match'));
            end
        end
    end
    
    % Clean up
    nisslFnums = nisslFnums(cellfun(@(x) ~isempty(x), nisslFnums));
    roiFnums = roiFnums(cellfun(@(x) ~isempty(x), roiFnums));
    
    % Sort and save into struct for output
    [~,nisslOrder] = sort(cell2mat(nisslFnums));
    [~,imageOrder] = sort(cell2mat(roiFnums));
    
    imageID.Nissl.Nums = nisslFnums(nisslOrder);
    imageID.Nissl.FNames = nisslFiles(nisslOrder);
    imageID.ROI.Nums = roiFnums(imageOrder);
    imageID.ROI.FNames = roiFiles(imageOrder);
  
% If more than one channels
% Look for number
elseif chanNum > 1
    
    imFnums = cellfun(@(x) str2double(regexp(x,'[0-9]+','match')), imFileNames);
    
    imFnums = num2cell(imFnums);
    
    % Sort and save into struct for output
    [~,imageOrder] = sort(cell2mat(imFnums));
    
    imageID.Image.Nums = imFnums(imageOrder);
    imageID.Image.FNames = imFileNames(imageOrder);
    
end






end

