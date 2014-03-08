function [Quadrants] = QSectPoly_v01(xCoords, yCoords, mask, imageChannel)

Quadrants = struct;

[B,~,~,~] = bwboundaries(mask);

pixelInfo = regionprops(mask,imageChannel,'Centroid');

% xmin = min(xCoords);
xmax = max(xCoords);
% ymin = min(yCoords);
ymax = max(yCoords);

% find corner vertices
xCor = xCoords(2:end);
yCor = yCoords(2:end);

yCorR = zeros(numel(yCor),1);
for i = 1:numel(yCor)
   yCorR(i) = yCor(i) + rand; 
end

yCor = yCorR;

xCorR = zeros(numel(xCor),1);
for i = 1:numel(xCor)
   xCorR(i) = xCor(i) + rand; 
end

xCor = xCorR;

%% Top Left corner
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

%% Bottom Left corner
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

%% Bottom Right corner
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
        
%% Top Right corner
topRight = (yCor < ymax*0.5 & xCor > xmax*0.5); 
yratio = 0.5;
xratio = 0.5;

if sum(topRight) ~= 1
    if sum(topRight) == 0
        
        if sum(yCor < ymax*0.5) == 0
            
            maxSortx = sort(xCor,'descend');
            i = 1;
            while sum(topRight) ~= 1 % too conservative
                maxNow = maxSortx(1:i);
                xVals = ismember(xCor, maxNow);
                yratio = yratio + 0.01;
                topRight = (xVals & yCor < ymax*yratio);
                i = 1 + 1;
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
            yratio = yratio - 0.01;
            xratio = xratio + 0.01;
            topRight = (yCor < ymax*yratio & xCor > xmax*xratio);
        end
    end
end

%% THIS IS CORRECT

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

%% Quadrant 1

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

quad_1 = [quad1Xfinal quad1Yfinal];

% Create polygon handles for quadrant 1

h = impoly(imageChannel, quad_1);

%% Quadrant 1 Mask
Quadrants.TL = createMask(h);

%% Quadrant 2

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

quad_2 = [quad2Xfinal quad2Yfinal];

% Create polygon handles for quadrant 2

h2 = impoly(gca, quad_2);

%% Quadrant 2 Mask
Quadrants.BL = createMask(h2);

%% Quadrant 3

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

quad_3 = [quad3Xfinal quad3Yfinal];

h3 = impoly(gca, quad_3);

%% Quadrant 3 Mask
Quadrants.BR = createMask(h3);

%% Quadrant 4

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

quad_4 = [quad4Xfinal quad4Yfinal];

h4 = impoly(gca, quad_4);

%% Quadrant 4 Mask
Quadrants.BL = createMask(h4);




