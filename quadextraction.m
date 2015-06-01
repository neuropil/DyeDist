function [output] = quadextraction(polyType, plotflag)
% QUADEXTRACTION MATLAB code for quadextraction.m
%      QUADEXTRACTION: takes specific dataset output from SClayers GUI and 
%      plots and exports orientation data.
%
%      H = quadextraction(polyType) returns struct variable with selected
%                                    orientations.
%      INPUTS: 
%             plotType: 
%                      'whole' : if user wants to examine whole polygon
%                      analysis.
%                      REQUIRES mfile with 'Whole' in file name
%                      Any string input other than 'whole' will initiate
%                      quadrant analysis
%                      REQUIRES mfile with 'Quadrant in file name
%
%            plotflag:
%                     0 or 1 to output plot
%                     NO INPUT will initiate DEFAULT: 0
%
%
%      EXAMPLES:
%                [out] = quadextraction('whole')
%                [out] = quadextraction('quad')
%                [out] = quadextraction('whole',1)
%
%
% John A. Thompson 2014
% email questions or bugs to john.arthur.thompson@gmail.com



if nargin == 0
    warndlg('Not enough input arguments')
    return
elseif nargin == 1
    plotflag = 0;
    w = warndlg('Plotting is turned off','A Warning Dialog');
    waitfor(w);
end


if strcmp(polyType,'whole')
    
    [fileName, Location] = uigetfile;
    
    cd(Location)
    
    load(fileName)
    
    workSpaceCheck = who;
    
    if ismember('outDS',workSpaceCheck) == 0;
       warndlg('For whole polygon analysis choose WholePoly');
       return
    end
    
    allML = outDS.Dye_Ratio;
    lateral = allML(1:round(numel(allML)/2));
    medial = allML(round(numel(allML)/2):end);
    mediallateral = [mean(lateral) , mean(medial)];
    medlatSD = [std(lateral), std(medial)];
    
    if plotflag
        
        bar(mediallateral,0.5)
        set(gca, 'XTickLabel',{'lateral','medial'})
        hold on
        errorbar(mediallateral,medlatSD,'rx')
        ylabel('Fraction of SC area with Chat Expression')
        
    end
    
    output.medial = medial;
    output.lateral = lateral;
    
else
    
    [fileName, Location] = uigetfile;
    
    cd(Location)
    
    load(fileName)
    
    workSpaceCheck = who;
    
    if ismember('quadOutData',workSpaceCheck) == 0;
       warndlg('For whole polygon analysis choose Quadrant');
       return
    end

    sectionS = unique(quadOutData.SectionID);
    
    % RESORT
    
    numbers = cellfun(@(x) str2double(x), regexp(sectionS,'[0-9]{1,2}','match'));
    
    [~, newSort] = sort(numbers);
    
    reOrderList = sectionS(newSort,:);
    
    a = zeros(length(sectionS),1);
    p = zeros(length(sectionS),1);
    v = zeros(length(sectionS),1);
    d = zeros(length(sectionS),1);
    
    letters = {'A', 'P', 'V', 'D'};
    for si = 1:length(sectionS)
        secIndex = strcmp(reOrderList{si},quadOutData.SectionID);
        
        tempDS = quadOutData(secIndex,:);
        
        
        for li = 1:length(letters)
            switch letters{li}
                
                case 'A'
                    aindex = zeros(1,2);
                    idCount = 1;
                    for ci = 1:length(tempDS)
                        if strfind(tempDS.QuadID{ci},'A')
                            aindex(idCount) = tempDS.AreaRatio{ci};
                            idCount = idCount + 1;
                        end
                    end
                    
                    a(si) = mean(aindex);
                    
                case 'P'
                    pindex = zeros(1,2);
                    idCount = 1;
                    for ci = 1:length(tempDS)
                        if strfind(tempDS.QuadID{ci},'P')
                            pindex(idCount) = tempDS.AreaRatio{ci};
                            idCount = idCount + 1;
                        end
                    end
                    
                    p(si) = mean(pindex);
                    
                case 'V'
                    vindex = zeros(1,2);
                    idCount = 1;
                    for ci = 1:length(tempDS)
                        if strfind(tempDS.QuadID{ci},'V')
                            vindex(idCount) = tempDS.AreaRatio{ci};
                            idCount = idCount + 1;
                        end
                    end
                    
                    v(si) = mean(vindex);
                    
                case 'D'
                    dindex = zeros(1,2);
                    idCount = 1;
                    for ci = 1:length(tempDS)
                        if strfind(tempDS.QuadID{ci},'D')
                            dindex(idCount) = tempDS.AreaRatio{ci};
                            idCount = idCount + 1;
                        end
                    end
                    
                    d(si) = mean(dindex);
            end
        end
    end
    
    
    %%
    
    anterposter = [mean(a) , mean(p)];
    stdantpost = [std(a) , std(p)];
    
    dorsVentral = [mean(d) , mean(v)];
    stddorven = [std(d) , std(v)];
    
    if plotflag
        
        figure;
        bar(anterposter,0.5);
        set(gca, 'XTickLabel',{'anterior','posterior'})
        hold on
        errorbar(anterposter,stdantpost,'rx')
        ylabel('Fraction of SC area with Chat Expression')
        figure;
        bar(dorsVentral,0.5);
        set(gca, 'XTickLabel',{'dorsal','ventral'})
        hold on
        errorbar(dorsVentral,stddorven,'rx')
        ylabel('Fraction of SC area with Chat Expression')
        
    end
    
    output.anterior = a;
    output.posterior = p;
    output.dorsal = d;
    output.ventral = v;
    
end

