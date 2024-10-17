% count_cells_stereoscope_Ngo2024.m
%
% Read 2D stereoscope images, detect cells by thresholding, 
% culling small and large objects. (Objects = immune cells.)
% Then high-pass filter and again threshold, etc.; combine (OR) the 
% two binary images.
%
% ------------------------------------------------------------------------------------------------
% For data and analysis described in Ngo et al. 2024
%     Copied "as used" with code dependent on particular file structures.
%     The function findObjects() takes an image as an input, and is
%     generally applicable, regardless of file structure.
% ------------------------------------------------------------------------------------------------
% 
% Threshold value automatically set to be = median + (median-min).
% 
% Direct counting of objects, and inference of number of objects from total 
% fluorescence intensity / intensity of median object (which assumes
% median object is a single cell)
% Also returns total median-subtracted intensity, median object intensity,
% and median object area
% Writes to CSV file: overall statistics for each image, and also for the
% entire set a CSV file with all the cell intensity information.
%
% Inputs:
%   fishNoArray : array of fish numbers to examine. (Required)
%       Enter just one number to examine one image dataset (Recommend: 
%             printOutput = true; showEdgeImage = true)
%       Enter an array to examine multiple images, probably
%             irf8-sgRNA injected: fishNoArray = 1:67;
%             uninjected WT fish:  fishNoArray = 68:87;
%             (Recommend: printOutput = false; showEdgeImage = false)
%   fileName : Default 'C3.tif' for mCherry/Macrophages
%   pathNameBase : parent folder containing \F{n}... fish folders
%       default 'C:\Users\Raghu\Documents\Temporary\Macrophages Stereomicroscope Irf8 Julia 21Mar2023\';
%       likely other location: 'E:\Macrophages Stereomicroscope Irf8 Julia 21Mar2023\';
%   scale : um/px; default 0.8125 um/px (8x Zoom)
%   radius_range_um : min and max cell radii to allow for objects, microns
%             Default 3, 20 um
%   sigma_HP : high pass filter size, px. Should be roughly equal to cell
%             diameter; default ceil(2*radius_range_um(1)/scale)
%   printOutput : print object finding output for each image; default 
%       true if fishNoArray is a single element; false otherwise
%   showEdgeImage : if true, superimpose edges of binary object image on
%                   original image; default true if fishNoArray
%                   is a single element; false otherwise
%   outputCSVfilename = filename for CSV output; will add csv.
%                   empty for none. (default)
%
% Outputs:
%     cellData: structure with number of elements = length(fishNoArray); Fields:
%         Nobjects : number of detected objects
%         Nobjects_inferred : number of objects inferred
%         totalIntensity : total (median-subtracted) intensity
%         medianObjectIntensity : total intensity of median object 
%         medianObjectArea: Area (um^2) of median object
%         thresh: intensity threshold used for this image
%
% Raghuveer Parthasarathy
% Dec. 22, 2023
% January 24, 2024: make more modular
% Last modified January 24, 2024

function cellData = count_cells_stereoscope_Ngo2024(fishNoArray, fileName, ...
    pathNameBase, scale, radius_range_um, sigma_HP, outputCSVfilename, printOutput, ...
    showEdgeImage)

if ~exist('fileName', 'var') || isempty(fileName)
    fileName = 'C3.tif'; % mCherry image
end
if ~exist('pathNameBase', 'var') || isempty(pathNameBase)
    pathNameBase = 'C:\Users\Raghu\Documents\Temporary\Macrophages Stereomicroscope Irf8 Julia 21Mar2023\';
end
if ~exist('scale', 'var') || isempty(scale)
    scale = 0.8125; % um/px (8x Zoom)
end
if ~exist('radius_range_um', 'var') || isempty(radius_range_um)
    radius_range_um = [3 20]; % um
end
if ~exist('sigma_HP', 'var') || isempty(sigma_HP)
    sigma_HP = ceil(2*radius_range_um(1)/scale);
end
if ~exist('outputCSVfilename', 'var') 
    outputCSVfilename = [];
end
if ~exist('printOutput', 'var') || isempty(printOutput)
    printOutput = length(fishNoArray)==1;
end
if ~exist('showEdgeImage', 'var') || isempty(showEdgeImage)
    showEdgeImage = length(fishNoArray)==1;
end

%% Examine each image

radius_range_px = radius_range_um/scale; % min, max to allow, pixels

nFish = length(fishNoArray);

cellData = repmat(struct('Nobjects', [], ...
    'totalIntensity', [], 'medianObjectIntensity', [], ...
    'medianObjectArea', []), [nFish 1]);

% All images
all_CellIntensities = [];
for j=1:nFish
    if mod(j,10)==0
        fprintf('Image %d ... ', j)
    end
    pathName = strcat(pathNameBase, 'F', num2str(fishNoArray(j)));
    im = imread(strcat(pathName, '\', fileName));
    [cellData(j).Nobjects, cellData(j).totalIntensity, ...
        cellData(j).medianObjectIntensity, ...
        medianObjectAreaPx, thisImage_TotalIntensities, im_bin_combined] = ...
        findObjects(im, radius_range_px, sigma_HP, printOutput, showEdgeImage);
    cellData(j).medianObjectArea = medianObjectAreaPx*scale*scale;
    all_CellIntensities = [all_CellIntensities thisImage_TotalIntensities];
end
fprintf('\n')

% histogram of all cell intensities
% figure; histogram(all_CellIntensities, 1000)
% xlabel('Cell Intensities')

% CSV output
if ~isempty(outputCSVfilename)
    outputMatrix = [fishNoArray', [cellData.totalIntensity]', [cellData.Nobjects]', ...
        [cellData.medianObjectIntensity]', [cellData.medianObjectArea]'];
    T = array2table(outputMatrix);
    T.Properties.VariableNames(1:5) = {'fishNo','totalIntensity', 'Nobjects', ...
        'medianObjectIntensity', 'MedianObjectArea_um2'};
    writetable(T,strcat(outputCSVfilename, '.csv'));

    % all cell intensities
    T2 = array2table(all_CellIntensities');
    T2.Properties.VariableNames = {'CellIntensities'};
    writetable(T2,strcat(outputCSVfilename, '_allCellIntensities', '.csv'));
end

end

function [Nobjects, totalIntensity, medianObjectIntensity, ...
          medianObjectArea, all_TotalIntensities, im_bin_combined] = ...
         findObjects(im, radius_range_px, sigma_HP, printOutput, showEdgeImage)
    % Direct counting of objects, and inferring number of objects from total 
    % fluorescence intensity / intensity of median object (which assumes
    % median object is a single cell)
    %
    % Processing: 
    %     - threshold as median + (median-min)
    %     - Closure; disk of radius 2 (hard-coded)
    %     - high-pass filtering, and then thresholding again
    %     - combine both above-threshold etc. binary images
    %
    % Inputs
    %     im : 2D image
    %     showEdgeImage: if true, superimpose edges of binary object image on
    %          original image
    %
    %     See above header for details about other inputs
    %
    % Outputs:
    %     Nobjects : number of objects (above threhold etc. connected
    %         regions) found
    %     totalIntensity : total median-subtracted intensity of the image
    %     medianObjectIntensity : median object intensity (image median
    %         subtracted)
    %     medianObjectArea : median object area (px)
    %     all_TotalIntensities : all the object (cell) total intensities
    %     im_bin_combined : binary image, combined (OR) from thresholded and
    %                   HP-filtered + thresholded images

    class_im = class(im);
    maxVal = intmax(class(im)); % max of this class, to use later.
    % figure('name', 'Original image'); imshow(im, [])

    % Total image intensity, median-subtracted
    median_intensity_im = double(median(im, "all"));
    min_intensity_im = double(min(im, [], 'all'));
    totalIntensity = sum(double(im)-median_intensity_im, 'all');

    % Threshold.
    % Median is essentially equal to background; signal occupies a small
    % fraction of pixels.
    % Threshold value = median + (median-min)
    thresh = 2*median_intensity_im - min_intensity_im;
    im_bin = im > thresh;
    % fprintf('Threshold: %.1f\n', thresh);
    % figure('name', 'above threshold'); imshow(im_bin, [])
    
    % Opening (smooth and separate objects)
    ste = strel('disk', 2);
    im_bin = imopen(im_bin, ste);
    % figure('name', 'above threshold; opened'); imshow(im_bin, [])
    
    % Remove small objects
    im_bin = bwareaopen(im_bin, round(pi*radius_range_px(1)*radius_range_px(1)));
    % figure('name', 'above threshold; remove small'); imshow(im_bin, [])
    
    % Remove large objects. Loop, so can fill in binary image also
    stats1 = regionprops(im_bin, 'area', 'PixelIdxList');
    delIndex = [];
    for j=1:length(stats1)
        if stats1(j).Area > (pi*radius_range_px(2)*radius_range_px(2))
            im_bin(stats1(j).PixelIdxList) = 0;
            delIndex = [delIndex j];
        end
    end
    % figure('name', 'above threshold; remove small and large'); imshow(im_bin, [])

    % High pass filter
    im_lp = imgaussfilt(im, sigma_HP);
    im_hp = double(im)-double(im_lp);
    % figure('name', 'high pass'); imshow(im_hp, [])
    % figure; histogram(im_hp(:), 1000)
    % set(gca, 'yscale', 'log')
    im_hp(im_hp<0) = 0.0;
    thresh_hp = median(im_hp, 'all') + 5*std(im_hp, [], 'all');
    im_hp_bin = im_hp>thresh_hp;
    % figure('name', 'high pass, threshold, original'); imshow(im_hp_bin, [])
    % Opening; same ste
    im_hp_bin = imopen(im_hp_bin, ste);
    % Remove small objects
    im_hp_bin = bwareaopen(im_hp_bin, round(pi*radius_range_px(1)*radius_range_px(1)));
    % figure('name', 'high pass, threshold, remove small'); imshow(im_hp_bin, [])
    
    % overall binary: combine first thresholded image and HP-filter-derived
    % image
    im_bin_combined = im_bin | im_hp_bin;
    % figure('name', 'overall binary'); imshow(im_bin_2, [])

    % Find objects
    % Use original image, median-subtracted, for intensity of objects
    im_for_intensity = double(im) - median_intensity_im;
    % background-subtracted intensities
    stats = regionprops(im_bin_combined, im_for_intensity, 'area', ...
        "Centroid", "MeanIntensity", "PixelIdxList");
    Nobjects = length(stats);

    % Total intensity of each object; also save to a combined array
    for j=1:Nobjects
        % make double to avoid possible overflow
        stats(j).objectTotalIntensity = stats(j).Area * double(stats(j).MeanIntensity);
    end
    all_TotalIntensities = [stats.objectTotalIntensity];

    % Infer number of objects: Total image fluorescence / median object
    % intensity, where median object intensity = median object area *
    % mean(mean object intensity)
    % Note that total image fluorescence has been median-subtracted
    if Nobjects>0
        medianObjectArea = median([stats.Area]);
        meanObjectIntensity = mean([stats.MeanIntensity]);
        medianObjectIntensity = medianObjectArea*meanObjectIntensity;
    else
        medianObjectArea = 0;
        medianObjectIntensity = 0;
    end
        
    if printOutput
        fprintf('\n')
        fprintf('Total median-subtracted intensity: %.3e\n', totalIntensity)
        fprintf('Number of objects found: %d\n', Nobjects)
        fprintf('median Object Intensity %.2e\n', medianObjectIntensity)
        fprintf('median Object Area (px) %.2e\n', medianObjectArea)
    end
    if showEdgeImage
        % Make an image with the final binary images as edges on the original
        % images
        im_bin_dilate = imdilate(im_bin_combined, ste);
        im_bin_edge = im_bin_dilate & ~im_bin_combined;
        % imwrite(uint8(255.0*im_bin_edge), 'im_bin_edge.png')
        max_im = double(max(im(:)));
        im_rescale = double(im)*double(maxVal)/max_im; % rescale to max of class
        im_colorEdge = repmat(im_rescale, [1 1 3]);
        im_colorEdge(:,:,1) = max(im_rescale, double(maxVal)*double(im_bin_edge));
        im_colorEdge(:,:,2) = max(im_rescale, double(maxVal)*double(im_bin_edge));
        figure('name', 'original with objects'); imshow(cast(im_colorEdge, class_im), [])
    end

end
