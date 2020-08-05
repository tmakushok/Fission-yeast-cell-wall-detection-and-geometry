%% Main steps of the script:
%% 1) Detecting the perimeter of each cell and finding it's axis and tips
%% 2) Calculating length and width of each cell
clear;      
%% Parameters
% Path to the files to analyze   
ImageBaseName = '_InputImages\';   
% Path and filename for output of the results
FileHistLen = '_OutputImagesAndVariables\HistLength.fig';
FileHistWid = '_OutputImagesAndVariables\HistWidth.fig';
FileResPix = '_OutputImagesAndVariables\ResultInPixels.mat';
FileResMicr = '_OutputImagesAndVariables\ResultInMicrons.mat';
FileResPixTx = '_OutputImagesAndVariables\ResultInPixels.txt';
FileResMicrTx = '_OutputImagesAndVariables\ResultInMicrons.txt';
FileWeirdShape = '_OutputImagesAndVariables\NbWeirdlyShapedCells.txt';
% Minimal area of a cell and cell contour, in pixels
MinArea = 1000;             
MinAreaCont = 100; 
% Defining the size limits for objects to be recognized as being S.pombe cells
MinCellWidth = 45;
MaxCellWidth = 90; 
MinCellLength = 90;
MaxCellLength = 175;
% Parameter used for filling in the space between the edges of the
% cells determined automaticly (in pixels)
MaxEdgeDistance = 100;
% Length of the lines along which to take intensity profiles
% to obtain maximum values for cell length
LineLength = (MaxCellLength + MinCellLength)/2;
% Length of the lines along which to take intensity profiles
% to obtain maximum values for cell width
LineLengthW = (MaxCellWidth + MinCellWidth)/2;
% The value added to the automaticly detected half cell width to
% get the initial position of line center for lines along which we take
% intensity profiles to obtain maximum values for width
HWPlus = 10;
% Maximal distance between centers of two detected areas when they
% are still considered as the same cell detected twice
DistDeDoubling = 45;
% Number of bins for the histogram
NbBins = 50;
% Pixel size
PxSize = 0.0628;
% Max number of bad cells that can be selected manually
GinputNb = 30;     
% Figure numbers for visualization
ImageFigNb = 100;
AllBinCellsFigNb = 200;
%% 
FigNb = 1;
TotalIntens = [];
AverageIntens = [];
Result = [];
TotInt = 0;
AverInt = 0;
ImFiles = dir(ImageBaseName);
for i_File = 3:length(ImFiles)              % Loop on the image files to analyse  
    % i_File starts at 3 because first two elements of AverageImFiles are
    % files '.' and '..'
    close all;        
    FileName = ImFiles(i_File).name;  
    FilePath = [ImageBaseName FileName];    
    InitImage = imread(FilePath);
    [m, n] = size(InitImage);       
    FullImWithOneCell = [];     % To accumulate segmented cells in one image
    figure, imshow(InitImage);
    imagesc(InitImage);    
    %% Determining the edges of the cells using 'edge'    
    [CannyMed13_2, ThresCanny] = edge(InitImage, 'canny');
    figure, imshow(CannyMed13_2);  
    %% Taking off very small contours        
    LabelsCont = bwlabel(CannyMed13_2);
    StatsCont = regionprops(LabelsCont, 'Area', 'FilledImage', 'BoundingBox', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
    FullImWithOneCell = zeros(m, n);        % To have an image of the full size, but with only one filled cell on it	
	% Selecting objects that are big enough
    [Condition] = find([StatsCont.Area] > MinAreaCont);
    BigCont = ismember(LabelsCont, Condition);    
%     figure, imshow(BigCont);           
    CellEnds = []; 
    CellWidthEnds = [];    
    AddedToRes = 0;        % To keep track of the number of elements added in 'Result' for this image
    CellsPixels = cell(0,0);
    [OldLenRes, a] = size(Result);
    for i = 1:length(StatsCont)            
        % Choosing all edges longer than a min and not close to image border              
        if (StatsCont(i).Area > MinAreaCont) & (StatsCont(i).BoundingBox(1) > 2) & (StatsCont(i).BoundingBox(2) > 2)
            %% Filling with 1s the space between edge lines horisontally and vertically   
            %% Horizontal fill
            [m_OneCont, n_OneCont] = size(StatsCont(i).FilledImage);
            InputEdges = StatsCont(i).FilledImage;
%             figure, imshow(InputEdges);            
            % Putting some 0s all around the edge of the image (otherwise, edge analysed is too close to the border of the image)
            AddCanvas = 2.5*MaxEdgeDistance;
            InputEdges = [zeros(AddCanvas, n_OneCont + 2*AddCanvas); 
                zeros(m_OneCont, AddCanvas), InputEdges, zeros(m_OneCont, AddCanvas); 
                zeros(AddCanvas, n_OneCont + 2*AddCanvas)];
            [m_OneCont, n_OneCont] = size(InputEdges);             
%             figure, imshow(InputEdges);             
            InsideFlag = 0;
            j_Sobel = 1;               
            BinaryImage = zeros(m_OneCont, n_OneCont);
            for i_Sobel = 1:m_OneCont           % loop on line numbers
                while j_Sobel < (n_OneCont - MaxEdgeDistance) % + ???      % number of the column inside of the line
                    if (InputEdges(i_Sobel, j_Sobel) == 1)
                        BinaryImage(i_Sobel, j_Sobel) = 1;
                        for j_Search = 1:MaxEdgeDistance        % Looking where is the next proximal edge
                            if (InputEdges(i_Sobel, j_Sobel + j_Search) == 1)
                                InsideFlag = 1;
                                break;
                            end
                        end

                        if (InsideFlag == 1)            % Filling internal part of the cell with 1s
                           for j_Output = 1:j_Search
                               BinaryImage(i_Sobel, j_Sobel + j_Output) = 1;
                           end  
                           InsideFlag = 0;
                        end
                        j_Sobel = j_Sobel + j_Search;
                    else 
                        j_Sobel = j_Sobel + 1;
                    end                
                end
                j_Sobel = 1;
            end       
%             figure, imshow(BinaryImage);
            %% Vertical lines fill
            InsideFlag = 0;
            i_Sobel = 1;                       
            for j_Sobel = 1:n_OneCont           % loop on column numbers
                while i_Sobel < (m_OneCont - MaxEdgeDistance) % + ???      % number of the line inside of the column
                    if (BinaryImage(i_Sobel, j_Sobel) == 1)
                        for i_Search = 1:MaxEdgeDistance        % Looking where is the next proximal edge
                            if (BinaryImage(i_Sobel + i_Search, j_Sobel) == 1)
                                InsideFlag = 1;
                                break;
                            end
                        end

                        if (InsideFlag == 1)            % Filling internal part of the cell with 1s
                           for i_Output = 0:i_Search
                               BinaryImage(i_Sobel + i_Output, j_Sobel) = 1;
                           end  
                           InsideFlag = 0;
                        end
                        i_Sobel = i_Sobel + i_Search;
                    else 
                        i_Sobel = i_Sobel + 1;
                    end                
                end
                i_Sobel = 1;
            end
%             figure, imshow(BinaryImage);
            %% Extracting coordinates of the upper left corner and width of the
            %% bounding box that contained the image with one single edge
            lc_x = StatsCont(i).BoundingBox(1) - 0.5; 
            lc_y = StatsCont(i).BoundingBox(2) - 0.5; 
            wid_x = StatsCont(i).BoundingBox(3); 
            wid_y = StatsCont(i).BoundingBox(4); 
            % Extracting the filled cell from BinaryImage inside the box 
            % corresponding to the initial BoundingBox
            CroppedIm = BinaryImage((AddCanvas + 1):(AddCanvas + wid_y), (AddCanvas + 1):(AddCanvas + wid_x));
            %% Adding next 'layer' to the segmented image: next cell on
            %% the black BkGd added to the ones accumulated previously            
            OneCellInNature = [zeros(lc_y - 1, n); 
                zeros(wid_y, lc_x - 1), CroppedIm, zeros(wid_y, n - (lc_x - 1) - wid_x); 
                zeros(m - (lc_y - 1) - wid_y, n)];            
            FullImWithOneCell = FullImWithOneCell + OneCellInNature;            
            figure(AllBinCellsFigNb); 
            imshow(FullImWithOneCell); 
            imagesc(FullImWithOneCell);
            %% Measuring parameters of the filled cells (before it was done for the cell contour only)   
            Labels = bwlabel(OneCellInNature);
            Stats = regionprops(Labels, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'PixelList');           
            %% Discarding cells that have an area that is too small
            if Stats.Area < MinArea
                continue
            end
            %% Measuring cell length and checking its conformity to normal yeast values
            % Calculation of the coordinates of cell ends in format [x1, y1; x2, y2]                        
            x0 = Stats.Centroid(1);
            y0 = Stats.Centroid(2);
            Angle = Stats.Orientation;            
            %% Finding one of cell tips: going along the major axis of the ellipse
            %% towards the center until finding the line with a point that is white on BigRegions
            BorderFlag = 0;
            HalfLength = (Stats.MajorAxisLength / 2) + HWPlus;
            LengthFlag = 0;            
            BigRegions = OneCellInNature;            
            while (HalfLength > 0)
				x1 = x0 + HalfLength * cosd(Angle);
                y1 = y0 - HalfLength * sind(Angle); 
                % Condition for the center of the line being inside the image 
                if (uint16(x1) > 0) && (uint16(x1) < n + 1) && (uint16(y1) > 0) && (uint16(y1) < m + 1)
                    % At first, calculations are in the normal XY axis, but with the center in the
                    % center of the cell
                    % Calculating coordinates of the center of this perpendicular line                    
                    LineCenter = [HalfLength * cosd(Angle), HalfLength * sind(Angle)];  
                    % Calculating coordinates of the ends of this perpendicular line
                    % It lies on the line parallel to cell axis, so angle is Angle
                    % The line has a fixed length: LineLength
                    Line1_x = [LineCenter(1) - (LineLengthW/2)*sind(Angle), LineCenter(1) + (LineLengthW/2)*sind(Angle)];
                    Line1_y = [LineCenter(2) + (LineLengthW/2)*cosd(Angle), LineCenter(2) - (LineLengthW/2)*cosd(Angle)];  
                    % Passage back to normal Oxy axis
                    CellCenter = Stats.Centroid;                    
                    Line1_x = Line1_x + CellCenter(1);
                    Line1_y = -Line1_y + CellCenter(2);               
                    % Visualizing the position of the line                    
%                     figure, imshow(OneCellInNature); 
%                     line(Line1_x, Line1_y);  
%                     line(CellCenter(1), CellCenter(2), 'Color', [0 .8 0], 'Marker', 'o');  
%                      line(LineCenter(1) + CellCenter(1), -LineCenter(2) +
%                      CellCenter(2), 'Color', [.8 0 0], 'Marker', 'o');  
                    % Not to have problems (NaN) when calculating intensity
                    % profile, take as coordinates of the line ends
                    % the min/max between the calculated values and the
                    % borders of the image
                    for i_NaN = 1:2
                       if (uint16(Line1_x(i_NaN)) < 1)
                           Line1_x(i_NaN) = 1;
                       end
                       if (uint16(Line1_x(i_NaN)) > n)
                           Line1_x(i_NaN) = n;
                       end                       
                       if (uint16(Line1_y(i_NaN)) < 1)
                           Line1_y(i_NaN) = 1;
                       end
                       if (uint16(Line1_y(i_NaN)) > m)
                           Line1_y(i_NaN) = m;
                       end                              
                    end
                    % Calculating intensity profiles for the line 
                    PlotProfile = (improfile(OneCellInNature, Line1_x, Line1_y))';
                    LengthFlag = max(PlotProfile);   % Flag showing that the 'width limit' of the cell was attained 
%                     figure, plot(PlotProfile, 's', 'MarkerSize', 3); grid on;                     
%                     title('Intensity profile');
%                     xlabel('Position, in pixels');
%                     ylabel('Intensity');                                                                                
                    if (LengthFlag == 1)   
                        ElemEnds = [x1, y1];
                        break
                    else
                        HalfLength = HalfLength - 1;
                    end                
                else
                    BorderFlag = 1;     % The cell is at the border of the image
                    break                   
                end
            end
            if (BorderFlag == 1)    % If the cell was at the border, not continue with the analysis
                continue
            end
            %% Finding the second true cell tip
            BorderFlag = 0;
            HalfLength = (Stats.MajorAxisLength / 2) + HWPlus;
            LengthFlag = 0;
            while (HalfLength > 0)
				x2 = x0 - HalfLength * cosd(Angle);
                y2 = y0 + HalfLength * sind(Angle);                         
                if (uint16(x2) > 0) && (uint16(x2) < n + 1) && (uint16(y2) > 0) && (uint16(y2) < m + 1)
                    % At first, calculations in the normal XY axis, but with the center in the
                    % center of the cell                    
                    LineCenter = [-HalfLength * cosd(Angle), -HalfLength * sind(Angle)];  
                    % Calculating coordinates of the ends of this perpendicular line
                    % It lies on the line parallel to cell axis, so angle is Angle
                    % The line has a fixed length: LineLength
                    Line1_x = [LineCenter(1) - (LineLengthW/2)*sind(Angle), LineCenter(1) + (LineLengthW/2)*sind(Angle)];
                    Line1_y = [LineCenter(2) + (LineLengthW/2)*cosd(Angle), LineCenter(2) - (LineLengthW/2)*cosd(Angle)];  
                    % Back to normal Oxy axis                                    
                    Line1_x = Line1_x + CellCenter(1);
                    Line1_y = -Line1_y + CellCenter(2);               
                    % Checking the position of the line
%                     figure(FigNb); FigNb = FigNb + 1;
%                     imshow(OneCellInNature); 
%                     line(Line1_x, Line1_y);  
%                     line(CellCenter(1), CellCenter(2), 'Color', [0 .8 0], 'Marker', 'o');  
%                     line(LineCenter(1) + CellCenter(1), -LineCenter(2) + CellCenter(2), 'Color', [.8 0 0], 'Marker', 'o');  
                    % Not to have problems (NaN) when calculating intensity
                    % profile, we take as coordinates of the line ends
                    % the min/max between the calculated values and the
                    % borders of the image
                    for i_NaN = 1:2
                       if (uint16(Line1_x(i_NaN)) < 1)
                           Line1_x(i_NaN) = 1;
                       end
                       if (uint16(Line1_x(i_NaN)) > n)
                           Line1_x(i_NaN) = n;
                       end                       
                       if (uint16(Line1_y(i_NaN)) < 1)
                           Line1_y(i_NaN) = 1;
                       end
                       if (uint16(Line1_y(i_NaN)) > m)
                           Line1_y(i_NaN) = m;
                       end                              
                    end
                    % Calculating intensity profiles for the line 
                    PlotProfile = (improfile(OneCellInNature, Line1_x, Line1_y))';     
                    LengthFlag = max(PlotProfile);   % Flag showing that we attained the 'width limit' of the cell
%                     figure, plot(PlotProfile, 's', 'MarkerSize', 3); grid on;                     
%                     title('Intensity profile');
%                     xlabel('Position, in pixels');
%                     ylabel('Intensity');                           
                    if (LengthFlag == 1)   
                        CellEnds = [CellEnds; [ElemEnds, x2, y2, CellCenter(1), CellCenter(2)]];                        
                        break
                    else
                        HalfLength = HalfLength - 1;
                    end
                else
                    BorderFlag = 1;     % The cell is at the border of the image
                    break                   
                end                    
            end                        
            if BorderFlag == 1      % The cell is at the border of the image               
                continue
            end
            
            CellLength = sqrt((x2-x1)^2 + (y2-y1)^2);
            [CE_Lin, a] = size(CellEnds);
            if (CellLength < MinCellLength) | (CellLength > MaxCellLength)                               
                CellEnds(CE_Lin, :) = [];            % Taking this cell from the matrix with cell ends of good cells
                continue 
            end    
            %% Measuring cell width and checking its conformity to normal values           
            % CellWidthEnds format for one StatsCont entry is [x1, y1, x2, y2, x_CellCenter, y_CellCenter])                                         
            %% Get maximum width end 1: Looking at 'intensity profiles' perpendicular
            %% to the detected width
            BorderFlag = 0;
            HalfWidth = (Stats.MinorAxisLength / 2) + HWPlus;
            WidthFlag = 0;                        
            while (HalfWidth > 0)
				x1 = x0 - HalfWidth * sind(Angle);
                y1 = y0 - HalfWidth * cosd(Angle);                  
                if (uint16(x1) > 0) & (uint16(x1) < n + 1) & (uint16(y1) > 0) & (uint16(y1) < m + 1) % To make sure coordinates are inside the image
                    % At first, calculations in the normal XY axis, but with the center in the
                    % center of the cell
                    LineCenter = [HalfWidth * cosd(Angle + 90), HalfWidth * sind(Angle + 90)];  
                    % Calculating coordinates of the ends of this perpendicular line
                    % It lies on the line parallel to cell axis, so angle is Angle
                    % The line has a fixed length: LineLength
                    Line1_x = [LineCenter(1) - (LineLength/2)*sind(Angle + 90), LineCenter(1) + (LineLength/2)*sind(Angle + 90)];
                    Line1_y = [LineCenter(2) + (LineLength/2)*cosd(Angle + 90), LineCenter(2) - (LineLength/2)*cosd(Angle + 90)];  
                    % Passage back to normal Oxy axis                  
                    Line1_x = Line1_x + CellCenter(1);
                    Line1_y = -Line1_y + CellCenter(2);               
                    % Checking the position of the line
%                     figure(FigNb); FigNb = FigNb + 1;
%                     imshow(OneCellInNature); 
%                     line(Line1_x, Line1_y);  
%                     line(CellCenter(1), CellCenter(2), 'Color', [0 .8 0], 'Marker', 'o');  
%                     line(LineCenter(1) + CellCenter(1), -LineCenter(2) + CellCenter(2), 'Color', [.8 0 0], 'Marker', 'o');  
                    % Not to have problems (NaN) when calculating intensity
                    % profile, we take as coordinates of the line ends
                    % the min/max between the calculated values and the
                    % borders of the image
                    for i_NaN = 1:2
                       if (uint16(Line1_x(i_NaN)) < 1)
                           Line1_x(i_NaN) = 1;
                       end
                       if (uint16(Line1_x(i_NaN)) > n)
                           Line1_x(i_NaN) = n;
                       end                       
                       if (uint16(Line1_y(i_NaN)) < 1)
                           Line1_y(i_NaN) = 1;
                       end
                       if (uint16(Line1_y(i_NaN)) > m)
                           Line1_y(i_NaN) = m;
                       end                              
                    end
                    % Calculating intensity profiles for the line 
                    PlotProfile = (improfile(OneCellInNature, Line1_x, Line1_y))';     
                    WidthFlag = max(PlotProfile);   % Flag showing that we attained the 'width limit' of the cell
%                     figure, plot(PlotProfile, 's', 'MarkerSize', 3); grid on;                     
%                     title('Intensity profile');
%                     xlabel('Position, in pixels');
%                     ylabel('Intensity');                     
                    if (WidthFlag == 1)    % End of cell width is reached                       
                        ElemWidthEnds = [x1, y1];                        
                        break
                    else
                        HalfWidth = HalfWidth - 1;
                    end                
                else
                    BorderFlag = 1;     % The cell is at the border of the image
                    break                   
                end
            end
            if BorderFlag == 1      % The cell is at the border of the image               
                continue
            end            
            %% Get maximum width end 2: Looking at 'intensity profiles' perpendicular
            %% to the detected width
            HalfWidth = (Stats.MinorAxisLength / 2) + HWPlus;
            BorderFlag = 0;
            WidthFlag = 0;
            while (HalfWidth > 0)
				x2 = x0 + HalfWidth * sind(Angle);
                y2 = y0 + HalfWidth * cosd(Angle);               
                if (uint16(x2) > 0) & (uint16(x2) < n + 1) & (uint16(y2) > 0) & (uint16(y2) < m + 1)
                    % At first, calculations in the cartesian XY axis, but with the center in the
                    % center of the cell
                    LineCenter = [-HalfWidth * cosd(Angle + 90), -HalfWidth * sind(Angle + 90)];
                    % Calculating coordinates of the ends of this perpendicular line
                    % It lies on the line parallel to cell axis, so angle is Angle
                    % The line has a fixed length: LineLength
                    Line1_x = [LineCenter(1) - (LineLength/2)*sind(Angle + 90), LineCenter(1) + (LineLength/2)*sind(Angle + 90)];
                    Line1_y = [LineCenter(2) + (LineLength/2)*cosd(Angle + 90), LineCenter(2) - (LineLength/2)*cosd(Angle + 90)];  
                    % Passage back to normal Oxy axis
                    %CellCenter = Stats.Centroid;                    
                    Line1_x = Line1_x + CellCenter(1);
                    Line1_y = -Line1_y + CellCenter(2);               
                    % Checking the position of the line
%                     figure, imshow(OneCellInNature); 
%                     line(Line1_x, Line1_y);  
%                     line(CellCenter(1), CellCenter(2), 'Color', [0 .8 0], 'Marker', 'o');  
%                     line(LineCenter(1) + CellCenter(1), -LineCenter(2) + CellCenter(2), 'Color', [.8 0 0], 'Marker', 'o');  
                    % Not to have problems (NaN) when calculating intensity
                    % profile, we take as coordinates of the line ends
                    % the min/max between the calculated values and the
                    % borders of the image
                    for i_NaN = 1:2
                       if (uint16(Line1_x(i_NaN)) < 1)
                           Line1_x(i_NaN) = 1;
                       end
                       if (uint16(Line1_x(i_NaN)) > n)
                           Line1_x(i_NaN) = n;
                       end                       
                       if (uint16(Line1_y(i_NaN)) < 1)
                           Line1_y(i_NaN) = 1;
                       end
                       if (uint16(Line1_y(i_NaN)) > m)
                           Line1_y(i_NaN) = m;
                       end                              
                    end
                    % Calculating intensity profiles for the line 
                    PlotProfile = (improfile(OneCellInNature, Line1_x, Line1_y))';     
                    WidthFlag = max(PlotProfile);   % Flag showing that the 'width limit' of the cell was attained 
%                     figure, plot(PlotProfile, 's', 'MarkerSize', 3); grid on;                     
%                     title('Intensity profile');
%                     xlabel('Position, in pixels');
%                     ylabel('Intensity');                                         
                    if (WidthFlag == 1)                          
                        CellWidthEnds = [CellWidthEnds; [ElemWidthEnds, x2, y2, CellCenter(1), CellCenter(2)]];
                        break
                    else
                        HalfWidth = HalfWidth - 1;
                    end
                else
                    BorderFlag = 1;     % The cell is at the border of the image
                    break                   
                end                    
            end                        
            if BorderFlag == 1      % The cell is at the border of the image
                continue
            end            
            CellWidth = sqrt((x2-x1)^2 + (y2-y1)^2);
            [CWE_Lin, a] = size(CellWidthEnds);
            [CE_Lin, a] = size(CellEnds);
            if (CellWidth < MinCellWidth) || (CellWidth > MaxCellWidth)                 
                CellWidthEnds(CWE_Lin, :) = [];      % Taking this cell from stored 'half-good' cells
                CellEnds(CE_Lin, :) = [];            
                continue
            end       
            %% Visualisation of geometrical parameters of the selected good cells 
%             % Visualisation of the cell length
%             figure(ImageFigNb);
%             line([CellEnds(CE_Lin, 1), CellEnds(CE_Lin, 3)], [CellEnds(CE_Lin, 2), CellEnds(CE_Lin, 4)]);            
%             % Visualisation of the cell width            
%             line([CellWidthEnds(CWE_Lin, 1), CellWidthEnds(CWE_Lin, 3)], [CellWidthEnds(CWE_Lin, 2), CellWidthEnds(CWE_Lin, 4)]);            
%             % Visualisation of the cell center        
%             line(StatsCont(i).Centroid(1), StatsCont(i).Centroid(2), 'Color', [.8 0 0], 'Marker', 'o');  
            %% Solving the problem of having two detected areas overlapping 
            %% (steming out from inner and outer edges of a cell) with slightly different widths 
            AddedFlag = 0;            
            [LenWE, a] = size(CellWidthEnds); 
            switch LenWE
                case 0
                    continue            
                case 1                                       
                    Result = [Result; CellLength, CellWidth, CellCenter];
                    AddedToRes = AddedToRes + 1;
                    CellsPixels = {Stats.PixelList};                     
                    continue
            end
            for i_UnDoubl = (LenWE - 1):-1:1  % back count as more chances to have an overlap with one of the 'closely previous' areas             
                % Calculation of the distance between current cell and the
                % cell number i_DeDoubling
                CentersDist = sqrt((CellWidthEnds(i_UnDoubl, 5) - CellCenter(1))^2 + (CellWidthEnds(i_UnDoubl, 6) - CellCenter(2))^2);                
                if CentersDist <= double(DistDeDoubling)    % Centers of the two cells are very close
                    CellWidthOld = sqrt((CellWidthEnds(i_UnDoubl, 3) - CellWidthEnds(i_UnDoubl, 1))^2 + (CellWidthEnds(i_UnDoubl, 4) - CellWidthEnds(i_UnDoubl, 2))^2);
                    if CellWidth < CellWidthOld
                        % Determining the index of the line in Results that
                        % should be removed. For this, find in 'Result' 
                        % the line with the same cell center position as in
                        % CellWidthEnds(i_UnDoubl).
                        % 'OldLenRes' is to take away the contribution in 'Result' of
                        % lines coming from previous images
                        [ResLen, a] = size(Result);
                        i_ToRemove_CP = find(Result((OldLenRes + 1):ResLen, 3) == CellWidthEnds(i_UnDoubl, 5) & Result((OldLenRes + 1):ResLen, 4) == CellWidthEnds(i_UnDoubl, 6));       
                        i_ToRemove = i_ToRemove_CP + OldLenRes;
                        Result(i_ToRemove, :) = [];                            
                        % Remove those lines also in geometry matrixes                        
                        CellWidthEnds(i_UnDoubl, :) = [];
                        CellEnds(i_UnDoubl, :) = [];
                        % Adding of the good information to 'Result'
                        Result = [Result; CellLength, CellWidth, CellCenter]; 
                        % Taking away from the pixels set
                        CellsPixels(i_ToRemove_CP) = [];
                        % Adding to the pixels set
                        CellsPixels = [CellsPixels; {Stats.PixelList}];                                                
                        AddedFlag = 1;
                        break
                    else                        
                        AddedFlag = 1;   
                        % Taking away the last cell's information from geometry arrays
                        [LinG, a] = size(CellEnds);
                        CellWidthEnds(LinG, :) = [];
                        CellEnds(LinG, :) = []; 
                        break
                    end
                end
            end
            if AddedFlag == 0   % There was no cells with cell centers so close that it could have been an overlap
                Result = [Result; CellLength, CellWidth, CellCenter];                                
                CellsPixels = [CellsPixels; {Stats.PixelList}];		
                AddedToRes = AddedToRes + 1;
            end             
        end             
    end     
    %% Visualisation of cell objects after the un-overlapping procedure
    FinalImage = zeros(m, n);
    % Filling with 1s the places where 'final good' cells are detected
    for i_vis = 1:length(CellsPixels)
        for i_OnePix = 1:length(CellsPixels{i_vis})     % Better to use 'size', but here it doesn't matter
            FinalImage(CellsPixels{i_vis}(i_OnePix, 2), CellsPixels{i_vis}(i_OnePix, 1)) = 1;
        end
    end
    % Representing borders of the segmented image overlaid with initial image 
    BWoutline = bwperim(FinalImage);   % To find perimeter pixels in binary image
    Segout = InitImage; 
    OutlineGreyLevel = max(InitImage(:));
    Segout(BWoutline) = OutlineGreyLevel;     
    figure, imshow(Segout, []);
    % Putting geometrical parameters on this overlaid image
    [Len_Res, a] = size(Result);     
    for i_vis = 1:(Len_Res - OldLenRes)     %:-1:(Len_Res - AddedToRes - OldLenRes + 1)
        % Visualisation of the cell length 
        line([CellEnds(i_vis, 1), CellEnds(i_vis, 3)], [CellEnds(i_vis, 2), CellEnds(i_vis, 4)]);            
        % Visualisation of the cell width            
        line([CellWidthEnds(i_vis, 1), CellWidthEnds(i_vis, 3)], [CellWidthEnds(i_vis, 2), CellWidthEnds(i_vis, 4)]);                    
    end
end
%% Save 'Result' array, with data in pixels format
save(FileResPix, 'Result');      
save(FileResPixTx, 'Result', '-ASCII');      
%% Conversion from pixels to microns. Saving the results.
for i_R = 1:Len_Res
    Result(i_R) = Result(i_R) * PxSize;
end
save(FileResMicr, 'Result');      
save(FileResMicrTx, 'Result', '-ASCII');  
%% Histogram visualisation of the maximum length and width of all the cells
%% on all images in the folder analysed and saving of the histograms
figure,
[n_hist, xout] = hist(Result(:, 1), NbBins);
bar(xout, n_hist);
title('Histogram of cells length'); 
xlabel('Cell length, in microns');
saveas(gcf, FileHistLen);

figure,
[n_hist, xout] = hist(Result(:, 2), NbBins);
bar(xout, n_hist);
title('Histogram of cells width');
xlabel('Cell width, in microns');
saveas(gcf, FileHistWid);
