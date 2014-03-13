% Get files

cd C:\Users\John\Desktop\BethSCMatlabOutput

data2L = 12;

% Medial to Lateral

allML = outDS.Dye_Ratio;
lateral = allML(1:round(numel(allML)/2));
medial = allML(round(numel(allML)/2):end);
mediallateral = [mean(lateral) , mean(medial)];
medlatSD = [std(lateral), std(medial)];
bar(mediallateral,0.5)
set(gca, 'XTickLabel',{'lateral','medial'})
hold on
errorbar(mediallateral,medlatSD,'rx')
ylabel('Fraction of SC area with Chat Expression')

%%

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
    secIndex = strcmp(sectionS{si},quadOutData.SectionID);
    
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


anterposter = [mean(a) , mean(p)];
stdantpost = [std(a) , std(p)];
bar(anterposter,0.5);
set(gca, 'XTickLabel',{'anterior','posterior'})
hold on
errorbar(anterposter,stdantpost,'rx')
ylabel('Fraction of SC area with Chat Expression')

dorsVentral = [mean(d) , mean(v)];
stddorven = [std(d) , std(v)];
bar(dorsVentral,0.5);
set(gca, 'XTickLabel',{'dorsal','ventral'})
hold on
errorbar(dorsVentral,stddorven,'rx')
ylabel('Fraction of SC area with Chat Expression')

