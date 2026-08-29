function adult_RHT_analysis_single_section_outliers_removed()
% ADULT_RHT_ANALYSIS_SINGLE_SECTION_OUTLIERS_REMOVED - Adult Opn4/RHT axon
% comparative analysis (Speer Lab), single optical section + outlier-excluded
% variant.
%
% Same pipeline as adult_RHT_analysis_single_section.m, but specific
% channel/animal combinations are excluded from the analysis (see
% config.exclude below):
%   - Green channel omitted for WT animals 927, 928, 929, 930, 931, 932
%   - Both channels omitted for KO animal 937
% Excluded data is dropped everywhere it would otherwise be used: the Step 1
% normalization save-out and montages, the per-genotype normalization limits,
% the 1D profiles, alignment, statistics, and the average/shifted/heatmap
% outputs. A green-excluded WT animal still contributes its RED channel.
%
% Operates on single optical sections (2 um Z depth) instead of MIPs. Reads the
% rotated/cropped RGB single-section TIFFs (1130x634), splits out the red &
% green RHT-tracer channels (blue/DAPI is dropped), and runs the full WT-vs-KO
% 1D profile comparison.
%
% Genotype is read from the filename: 'c57' or 'wtmate' = WT, 'ko' = KO.
% Animal number = the leading digits of the filename (independent replicate).
% Filenames end in the section number (..._RGB_Z###); the section number is
% variable and does not affect genotype/animal parsing.
%
% Complete integrated analysis pipeline
% 
% Performs complete analysis in one script:
%   STEP 1: Image normalization (linear stretch)
%   STEP 2: 1D histogram computation (bin 80 + bin 20)
%   STEP 3: Alignment and statistical analysis
%
% Features:
%   - Single configuration section
%   - Optional step skipping
%   - Dual-bin strategy (bin 80 for alignment, bin 20 for plotting)
%   - Mann-Whitney U tests on per-animal summary statistics
%   - Post-hoc power analysis (asterisks when p<0.05 AND power>0.8)
%   - Enhanced legends with Mean, SD, and N
%   - Heatmaps of WT vs KO intensity differences
%   - Average images from shifted data
%   - Multiple output formats (.png, .fig, .eps)
%   - Shifted original images
%
% Author: Integrated pipeline (Steps 1-3 combined)
% Date: 2024 (Modified)

%% ===================== MASTER CONFIGURATION =====================
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║     SCN AXON ANALYSIS PIPELINE - MODIFIED VERSION         ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

config = struct();

% ===== PATHS =====
% Input: rotated/cropped RGB single-section TIFFs (1130x634), one per animal.
% Single optical sections (2 um Z depth), NOT max-intensity projections. Red &
% green are the RHT tracers; blue (DAPI) is dropped on load.
config.imgDir = "C:\Users\user\RHT_project\Analysis\Single_section_TIFFs_rot_crop";
config.outRoot = "C:\Users\user\RHT_project\Analysis\Final_1D_analysis_single_section";

% ===== STEP CONTROL =====
config.runStep1 = true;   % Image normalization
config.runStep2 = true;   % Histogram computation
config.runStep3 = true;   % Alignment and statistics

% ===== IMAGE PARAMETERS =====
config.pixelSize = 0.6214809;  % microns per pixel (this dataset's voxel X/Y)
config.scaleBarMicrons = 50;   % Scale bar length in microns

% ===== BINNING PARAMETERS (DUAL-BIN STRATEGY) =====
config.binSize_alignment = 80;  % Bin size for robust alignment (valley detection)
config.binSize_analysis = 20;   % Bin size for smooth plotting

% ===== ALIGNMENT PARAMETERS =====
config.maxShift = 100;              % Maximum shift allowed in pixels
config.outlierThreshold = 2.0;      % Standard deviations for outlier detection
config.smoothWindow = 1;            % Smoothing for valley detection (1 = no smoothing)
config.valleySearchMin = 480;       % Minimum X position for valley search (pixels)
config.valleySearchMax = 670;       % Maximum X position for valley search (pixels)
% Images are 1130 x 634; midline intensity dip sits ~500-650 px. Window above
% brackets that with a small margin. (Valley search is on the X axis, width
% 1130, unchanged from the MIP data; the height change does not affect it.)

% ===== STATISTICAL PARAMETERS =====
config.alpha = 0.05;                % Significance level for Mann-Whitney tests
config.powerThreshold = 0.8;        % Minimum power threshold for asterisks

% ===== LATERALITY MAPPING =====
% Each channel reports the projection of one eye (dual-eye injection of spectrally
% distinct tracers). The SCN lobe ipsilateral to the injected eye lies on a different
% half of the image for each eye, so the mapping from image half to laterality is
% CHANNEL-SPECIFIC and cannot be a single fixed label.
%   'negative' = x < 0 (left of the aligned midline)
%   'positive' = x > 0 (right of the aligned midline)
config.ipsiHalf.Red   = 'negative';  % red channel: ipsilateral SCN is x < 0
config.ipsiHalf.Green = 'positive';  % green channel: ipsilateral SCN is x > 0

% ===== OUTLIER EXCLUSIONS =====
% Per-channel lists of animal numbers to omit from the analysis. An excluded
% (channel, animal) pair is dropped everywhere: normalization limits, 1D
% profiles, alignment, statistics, and average/shifted/heatmap images. An
% animal excluded from only one channel still contributes the other channel.
% To exclude BOTH channels for an animal, list it under both Red and Green.
% Leave a list empty ([]) to exclude nothing for that channel.
config.exclude.Red   = [937];                         % KO 937: both channels omitted
config.exclude.Green = [927 928 929 930 931 932 937]; % WT green outliers + KO 937

% ===== OUTPUT PARAMETERS =====
config.saveShiftedImages = true;    % Save shifted original images
config.saveHeatmaps = true;         % Save heatmaps of WT vs KO differences
config.saveAverageImages = true;    % Save average images from shifted data
config.imageFormat = 'tif';         % Format for shifted images

% ===== FIGURE PARAMETERS =====
config.showFigures = true;          % Display figures during creation
config.savePNG = true;              % Save PNG versions
config.saveFIG = true;              % Save .fig files
config.saveEPS = true;              % Save .eps files
config.closeFIG = true;             % Close .fig files after saving
config.figureResolution = 600;      % DPI for saved figures (PNG and EPS)

% ===== OUTPUT DIRECTORIES =====
config.linearStretchDir = fullfile(config.outRoot, '01_Linear_Stretch');
config.histogramDir = fullfile(config.outRoot, '02_Histograms');
config.alignmentDir = fullfile(config.outRoot, '03_Alignment_with_Stats');
config.shiftedImagesDir = fullfile(config.outRoot, '04_Shifted_Images');
config.heatmapsDir = fullfile(config.outRoot, '05_Heatmaps');
config.averageImagesDir = fullfile(config.outRoot, '06_Average_Images');
% Parallel outputs using the maximum-intensity projection ACROSS ALIGNED
% ANIMALS (instead of the mean across animals). Same displays and heatmaps,
% built from the per-genotype/channel MIP of the aligned stack.
config.heatmapsMipDir = fullfile(config.outRoot, '05b_Heatmaps_MIP');
config.averageImagesMipDir = fullfile(config.outRoot, '06b_Average_Images_MIP');
% Parallel outputs using the STANDARD DEVIATION across aligned animals
% (like ImageJ Z-project -> Standard Deviation). Highlights where signal is
% variable across animals rather than where it is bright.
config.heatmapsStdDir = fullfile(config.outRoot, '05d_Heatmaps_STD');
config.averageImagesStdDir = fullfile(config.outRoot, '06d_Average_Images_STD');

% ===== INTERNAL PATHS (auto-generated) =====
config.normalizedDataFile = fullfile(config.outRoot, 'normalized_data.mat');
config.histogramDataFile = fullfile(config.outRoot, 'histogram_data.mat');
config.alignedDataFile = fullfile(config.outRoot, 'aligned_data_with_stats.mat');

%% ===================== RUN PIPELINE =====================
fprintf('Configuration loaded:\n');
fprintf('  Input directory: %s\n', config.imgDir);
fprintf('  Output directory: %s\n', config.outRoot);
fprintf('  Pixel size: %.4f µm/pixel\n', config.pixelSize);
fprintf('  Alignment bin size: %d pixels\n', config.binSize_alignment);
fprintf('  Analysis/plotting bin size: %d pixels\n', config.binSize_analysis);
fprintf('\n');

% Create main output directory
if ~exist(config.outRoot, 'dir')
    mkdir(config.outRoot);
end

% Run steps
if config.runStep1
    fprintf('═══════════════════════════════════════════════════════════\n');
    fprintf(' STEP 1/3: IMAGE NORMALIZATION\n');
    fprintf('═══════════════════════════════════════════════════════════\n');
    step1_ImageNormalization(config);
end

if config.runStep2
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════\n');
    fprintf(' STEP 2/3: HISTOGRAM COMPUTATION (DUAL-BIN)\n');
    fprintf('═══════════════════════════════════════════════════════════\n');
    step2_HistogramComputation(config);
end

if config.runStep3
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════\n');
    fprintf(' STEP 3/3: ALIGNMENT & STATISTICAL ANALYSIS\n');
    fprintf('═══════════════════════════════════════════════════════════\n');
    step3_AlignmentAndStatistics(config);
end

% Final summary
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              PIPELINE COMPLETE!                            ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Output locations:\n');
if config.runStep1
    fprintf('  Step 1: %s\n', config.linearStretchDir);
end
if config.runStep2
    fprintf('  Step 2: %s\n', config.histogramDir);
end
if config.runStep3
    fprintf('  Step 3: %s\n', config.alignmentDir);
    fprintf('  Shifted images: %s\n', config.shiftedImagesDir);
    if config.saveHeatmaps
        fprintf('  Heatmaps: %s\n', config.heatmapsDir);
    end
    if config.saveAverageImages
        fprintf('  Average images: %s\n', config.averageImagesDir);
    end
    fprintf('  Documentation: %s\n', fullfile(config.alignmentDir, 'Documentation'));
    fprintf('  Final data: %s\n', config.alignedDataFile);
end
fprintf('\n');

end

%% ═════════════════════════════════════════════════════════════
%% STEP 1: IMAGE NORMALIZATION
%% ═════════════════════════════════════════════════════════════

function step1_ImageNormalization(config)
% Load and normalize images using linear histogram stretching

fprintf('[1/4] Setting up directories...\n');
setupStep1Directories(config);

fprintf('[2/4] Loading RGB single-section images and splitting channels...\n');
channelData = loadExtractedChannels(config);

fprintf('[3/4] Performing linear histogram stretch (1%%-99%% saturation)...\n');
linearData = performLinearStretch(channelData, config);

fprintf('[4/4] Creating and saving montages...\n');
createStep1Montages(linearData, config);

fprintf('Saving normalized data...\n');
saveNormalizedData(channelData, linearData, config);

fprintf('✓ Step 1 Complete!\n');
end

function setupStep1Directories(config)
if ~exist(config.linearStretchDir, 'dir'), mkdir(config.linearStretchDir); end

subdirs = {
    fullfile(config.linearStretchDir, 'WT_Red')
    fullfile(config.linearStretchDir, 'WT_Green')
    fullfile(config.linearStretchDir, 'KO_Red')
    fullfile(config.linearStretchDir, 'KO_Green')
    fullfile(config.linearStretchDir, 'Montages')
};

for i = 1:length(subdirs)
    if ~exist(subdirs{i}, 'dir')
        mkdir(subdirs{i});
    end
end
end

function channelData = loadExtractedChannels(config)
fprintf('  Loading from: %s\n', config.imgDir);

channelData = struct();
channelData.Red = struct();
channelData.Green = struct();

% Input is one RGB single-section TIFF per animal (red & green = RHT tracers,
% blue = DAPI, dropped here). Filenames end in the section number, e.g.
% 928b_c57_t1_RGB_Z003. Glob the RGB single-section files.
rgb_files = dir(fullfile(config.imgDir, '*_RGB_Z*.tif'));
fprintf('  Found %d RGB single-section files\n', length(rgb_files));

allRed = [];
allGreen = [];
imageInfo = [];
WT_animals = [];
KO_animals = [];

for i = 1:length(rgb_files)
    [~, basename, ~] = fileparts(rgb_files(i).name);

    % Determine genotype from filename tags:
    %   WT = 'c57' or 'wtmate'   |   KO = 'ko' (case-insensitive)
    % Check WT tags first so a 'wtmate' file is never mis-caught by 'ko'.
    if contains(basename, 'c57', 'IgnoreCase', true) || ...
       contains(basename, 'wtmate', 'IgnoreCase', true)
        genotype = 'WT';
    elseif contains(basename, 'ko', 'IgnoreCase', true)
        genotype = 'KO';
    else
        warning('Could not determine genotype for %s, skipping', basename);
        continue;
    end

    % Animal number = leading digits of the filename (e.g. 931, 933b -> 933,
    % 5186b -> 5186). Each distinct number is an independent biological
    % replicate. A trailing letter after the digits is ignored.
    tokens = regexp(basename, '^(\d+)', 'tokens', 'once');
    if ~isempty(tokens)
        animalNum = str2double(tokens{1});
    else
        warning('Could not extract animal number from %s, skipping', basename);
        continue;
    end

    % Load the RGB single-section TIFF and split channels
    rgb_path = fullfile(config.imgDir, rgb_files(i).name);
    rgb = imread(rgb_path);

    if ndims(rgb) ~= 3 || size(rgb, 3) < 3
        warning('%s is not a 3-channel RGB image, skipping', basename);
        continue;
    end

    redChannel   = double(rgb(:, :, 1));   % red   = RHT tracer
    greenChannel = double(rgb(:, :, 2));   % green = RHT tracer
    % blue (:,:,3) = DAPI, intentionally dropped

    allRed = cat(3, allRed, redChannel);
    allGreen = cat(3, allGreen, greenChannel);
    imageInfo = [imageInfo; struct('genotype', genotype, 'animal', animalNum, 'filename', basename)];

    if strcmp(genotype, 'WT')
        WT_animals = [WT_animals, animalNum];
    else
        KO_animals = [KO_animals, animalNum];
    end

    fprintf('    Loaded: %s (%s animal %d)\n', basename, genotype, animalNum);
end

channelData.Red.images = allRed;
channelData.Green.images = allGreen;
channelData.imageInfo = imageInfo;
channelData.WT_animals = unique(WT_animals);
channelData.KO_animals = unique(KO_animals);

fprintf('  Total: %d images (%d WT, %d KO)\n', size(allRed, 3), ...
        length(channelData.WT_animals), length(channelData.KO_animals));
end

function linearData = performLinearStretch(channelData, config)
linearData = struct();
linearData.imageInfo = channelData.imageInfo;

channels = {'Red', 'Green'};

for ch = 1:length(channels)
    channelName = channels{ch};
    images = channelData.(channelName).images;
    
    fprintf('  Processing %s channel...\n', channelName);
    
    WT_idx = find(strcmp({channelData.imageInfo.genotype}, 'WT'));
    KO_idx = find(strcmp({channelData.imageInfo.genotype}, 'KO'));

    % Indices that actually inform the normalization reference: excluded
    % (channel, animal) pairs must not shift the stretch limits for the rest.
    % The stretch is still APPLIED to every image of the genotype so the
    % stretched_all array stays index-aligned with imageInfo; excluded images
    % are dropped later at the profile/stats stage.
    WT_keep = applyExclusions(WT_idx, channelName, channelData.imageInfo, config);
    KO_keep = applyExclusions(KO_idx, channelName, channelData.imageInfo, config);

    % WT normalization (reference computed from kept images only)
    WT_images = images(:, :, WT_idx);
    [WT_stretched, WT_lims, WT_gmin, WT_gmax] = ...
        stretchBlock(WT_images, images(:, :, WT_keep));

    % KO normalization (reference computed from kept images only)
    KO_images = images(:, :, KO_idx);
    [KO_stretched, KO_lims, KO_gmin, KO_gmax] = ...
        stretchBlock(KO_images, images(:, :, KO_keep));

    if numel(WT_keep) < numel(WT_idx) || numel(KO_keep) < numel(KO_idx)
        fprintf(['    %s: normalization reference excludes %d WT + %d KO ' ...
                 'outlier image(s)\n'], channelName, ...
                 numel(WT_idx) - numel(WT_keep), numel(KO_idx) - numel(KO_keep));
    end
    % Reconstruct
    stretched_all = zeros(size(images));
    stretched_all(:, :, WT_idx) = WT_stretched;
    stretched_all(:, :, KO_idx) = KO_stretched;
    
    linearData.(channelName).images = stretched_all;
    
    % Save individual images (excluded channel/animal pairs are skipped)
    saveNormalizedImages(stretched_all, channelData.imageInfo, channelName, ...
                        config.linearStretchDir, config);
end
end

function [stretched, lims, gmin, gmax] = stretchBlock(allImages, refImages)
% Linear-stretch a block of images to 0-255. The normalization REFERENCE
% (the mat2gray min/max and the 1st/99th percentile limits) is computed from
% refImages only, so excluded outliers do not perturb it, but the same fixed
% mapping is applied to every image in allImages. This reproduces the original
% mat2gray -> prctile(1,99) -> imadjust*255 behaviour when refImages ==
% allImages.
if isempty(refImages)
    % Degenerate: nothing to build a reference from; fall back to all images.
    refImages = allImages;
end

% mat2gray reference range from the kept images
gmin = min(refImages(:));
gmax = max(refImages(:));
if gmax <= gmin
    gmax = gmin + 1;   % avoid divide-by-zero on a flat block
end

% Normalize the reference the same way mat2gray would, then take percentiles
ref_norm = (double(refImages) - gmin) / (gmax - gmin);
ref_norm = min(max(ref_norm, 0), 1);
lims = [prctile(ref_norm(:), 1), prctile(ref_norm(:), 99)];

% Apply the SAME normalization + stretch to every image in the block
stretched = zeros(size(allImages));
for i = 1:size(allImages, 3)
    img = (double(allImages(:, :, i)) - gmin) / (gmax - gmin);
    img = min(max(img, 0), 1);
    stretched(:, :, i) = imadjust(img, lims, [0 1]) * 255;
end
end

function saveNormalizedImages(images, imageInfo, channelName, baseDir, config)
% Write one stretched TIFF per image, but skip any (channel, animal) pair that
% is excluded via config.exclude, so the normalization folder on disk contains
% only the animals that actually enter the analysis for this channel.
excludeAnimals = [];
if nargin >= 5 && isfield(config, 'exclude') && isfield(config.exclude, channelName)
    excludeAnimals = config.exclude.(channelName);
end

nSaved = 0;
nSkipped = 0;
for i = 1:size(images, 3)
    genotype = imageInfo(i).genotype;
    animal = imageInfo(i).animal;

    if ismember(animal, excludeAnimals)
        fprintf('      excluded from %s save-out: %s animal %d\n', ...
                channelName, genotype, animal);
        nSkipped = nSkipped + 1;
        continue;
    end

    outDir = char(fullfile(baseDir, [genotype '_' channelName]));
    outFile = char(fullfile(outDir, sprintf('%s_%s_animal%d_LinearStretch.tif', ...
                       genotype, channelName, animal)));
    
    imwrite(uint8(images(:, :, i)), outFile);
    nSaved = nSaved + 1;
end

fprintf('    Saved %d %s images (%d excluded)\n', nSaved, channelName, nSkipped);
end

function createStep1Montages(linearData, config)
channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};

for ch = 1:length(channels)
    channelName = channels{ch};
    
    for g = 1:length(genotypes)
        genotype = genotypes{g};
        
        idx = find(strcmp({linearData.imageInfo.genotype}, genotype));
        % Drop excluded (channel, animal) pairs so the montage matches the
        % saved TIFFs and the downstream analysis for this channel.
        idx = applyExclusions(idx, channelName, linearData.imageInfo, config);
        animalNumbers = [linearData.imageInfo(idx).animal];
        linear_images = linearData.(channelName).images(:, :, idx);
        
        labeled_images = addAnimalLabelsToImages(linear_images, animalNumbers);
        
        fig = figure('Position', [100 100 1200 800]);
        montage(labeled_images, 'Size', [NaN 5], 'BorderSize', 2);
        title(sprintf('%s %s - Linear Stretch', genotype, channelName), 'FontSize', 14);
        
        drawnow;
        
        outFileBase = char(fullfile(config.linearStretchDir, 'Montages', ...
                              sprintf('Montage_%s_%s_LinearStretch', genotype, channelName)));
        savefig(fig, [outFileBase '.fig']);
        print(fig, [outFileBase '.png'], '-dpng', sprintf('-r%d', config.figureResolution));
        try
            exportgraphics(fig, [outFileBase '.eps'], 'ContentType', 'vector');
        catch ME
            warning('EPS export failed for %s (PNG/FIG still saved OK): %s', outFileBase, ME.message);
        end
        
        try
            close(fig);
        catch
            % Figure handle already invalid (can happen over remote desktop) - nothing to close
        end
    end
end

fprintf('  Montages created\n');
end

function labeled_images = addAnimalLabelsToImages(images, animalNumbers)
[imgHeight, imgWidth, nImages] = size(images);
labeled_images = zeros(imgHeight, imgWidth, 3, nImages, 'uint8');

for i = 1:nImages
    img_gray = uint8(images(:, :, i));
    img_rgb = cat(3, img_gray, img_gray, img_gray);
    
    label_text = sprintf('Animal %d', animalNumbers(i));
    position = [15, 15];
    img_labeled = insertText(img_rgb, position, label_text, ...
                            'FontSize', 48, 'Font', 'Arial Bold', ...
                            'BoxColor', 'black', 'BoxOpacity', 0.9, ...
                            'TextColor', 'yellow', 'AnchorPoint', 'LeftTop');
    
    labeled_images(:, :, :, i) = img_labeled;
end
end

function saveNormalizedData(channelData, linearData, config)
data.channelData = channelData;
data.linearData = linearData;
data.config = config;

save(char(config.normalizedDataFile), 'data', '-v7.3');
fprintf('  Saved: %s\n', config.normalizedDataFile);
end

%% ═════════════════════════════════════════════════════════════
%% STEP 2: HISTOGRAM COMPUTATION
%% ═════════════════════════════════════════════════════════════

function step2_HistogramComputation(config)
% Compute 1D histograms at both bin sizes

fprintf('[1/3] Setting up directories...\n');
setupStep2Directories(config);

fprintf('[2/3] Loading normalized data...\n');
data = load(char(config.normalizedDataFile));
imageInfo = data.data.channelData.imageInfo;
normalizedData = data.data.linearData;

fprintf('[3/3] Computing dual-bin profiles (bin %d + bin %d)...\n', ...
        config.binSize_alignment, config.binSize_analysis);
profileData = computeDualBinProfiles(normalizedData, imageInfo, config);

createStep2Montages(profileData, config);

saveHistogramData(profileData, imageInfo, config);

fprintf('✓ Step 2 Complete!\n');
end

function setupStep2Directories(config)
if ~exist(config.histogramDir, 'dir'), mkdir(config.histogramDir); end

subdirs = {
    fullfile(config.histogramDir, 'Linear', 'Montages')
};

for i = 1:length(subdirs)
    if ~exist(subdirs{i}, 'dir')
        mkdir(subdirs{i});
    end
end
end

function profileData = computeDualBinProfiles(normalizedData, imageInfo, config)
channels = {'Red', 'Green'};
profileData = struct();

for ch = 1:length(channels)
    channelName = channels{ch};
    images = normalizedData.(channelName).images;
    
    fprintf('  Processing %s channel...\n', channelName);
    
    % Compute at both bin sizes
    [profiles_align, x_coords_align] = computeBinnedProfilesX(images, config.binSize_alignment);
    [profiles_analysis, x_coords_analysis] = computeBinnedProfilesX(images, config.binSize_analysis);
    
    WT_idx = find(strcmp({imageInfo.genotype}, 'WT'));
    KO_idx = find(strcmp({imageInfo.genotype}, 'KO'));

    % Drop excluded (channel, animal) pairs. This is the stage that feeds all
    % downstream profiles, alignment, and statistics, so filtering here removes
    % the outliers from every WT-vs-KO comparison for this channel.
    WT_idx = applyExclusions(WT_idx, channelName, imageInfo, config);
    KO_idx = applyExclusions(KO_idx, channelName, imageInfo, config);
    profileData.(channelName).WT.profiles_align = profiles_align(:, WT_idx);
    profileData.(channelName).WT.x_coords_align = x_coords_align;
    profileData.(channelName).WT.animals = [imageInfo(WT_idx).animal];
    
    profileData.(channelName).KO.profiles_align = profiles_align(:, KO_idx);
    profileData.(channelName).KO.x_coords_align = x_coords_align;
    profileData.(channelName).KO.animals = [imageInfo(KO_idx).animal];
    
    % Store analysis data (bin 20)
    profileData.(channelName).WT.profiles_analysis = profiles_analysis(:, WT_idx);
    profileData.(channelName).WT.x_coords_analysis = x_coords_analysis;
    
    profileData.(channelName).KO.profiles_analysis = profiles_analysis(:, KO_idx);
    profileData.(channelName).KO.x_coords_analysis = x_coords_analysis;
    
    fprintf('    Bin %d: %d bins, Bin %d: %d bins\n', ...
            config.binSize_alignment, length(x_coords_align), ...
            config.binSize_analysis, length(x_coords_analysis));
end
end

function [profiles_X, x_coords] = computeBinnedProfilesX(images, binSize)
[height, width, numImages] = size(images);

numBins = floor(width / binSize);
profiles_X = zeros(numBins, numImages);
x_coords = zeros(numBins, 1);

for b = 1:numBins
    startCol = (b-1) * binSize + 1;
    endCol = min(b * binSize, width);
    x_coords(b) = (startCol + endCol) / 2;
    
    for img = 1:numImages
        binRegion = images(:, startCol:endCol, img);
        profiles_X(b, img) = mean(binRegion(:));
    end
end
end

function createStep2Montages(profileData, config)
channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};

for ch = 1:length(channels)
    channelName = channels{ch};
    
    for g = 1:length(genotypes)
        genotype = genotypes{g};
        
        profiles = profileData.(channelName).(genotype).profiles_analysis;
        x_coords = profileData.(channelName).(genotype).x_coords_analysis;
        animals = profileData.(channelName).(genotype).animals;
        
        fig = figure('Position', [100 100 1000 600]);
        hold on;
        
        if strcmp(genotype, 'WT')
            colors = lines(50);
        else
            colors = jet(50);
        end
        
        for i = 1:size(profiles, 2)
            colorIdx = mod(i-1, size(colors, 1)) + 1;
            plot(x_coords, profiles(:, i), 'LineWidth', 1.5, ...
                 'Color', colors(colorIdx, :), ...
                 'DisplayName', sprintf('Animal %d', animals(i)));
        end
        
        mean_profile = mean(profiles, 2);
        std_profile = std(profiles, 0, 2);
        
        plot(x_coords, mean_profile, 'k-', 'LineWidth', 3, 'DisplayName', 'Mean');
        plot(x_coords, mean_profile + std_profile, 'k--', 'LineWidth', 2, 'DisplayName', '+1 SD');
        plot(x_coords, mean_profile - std_profile, 'k--', 'LineWidth', 2, 'DisplayName', '-1 SD');
        
        xlabel('X Position (pixels)', 'FontSize', 12);
        ylabel('Mean Intensity', 'FontSize', 12);
        title(sprintf('%s %s - Bin %d (n=%d)', genotype, channelName, ...
                     config.binSize_analysis, size(profiles, 2)), 'FontSize', 14);
        legend('Location', 'eastoutside');
        grid on;
        
        outFileBase = char(fullfile(config.histogramDir, 'Linear', 'Montages', ...
                              sprintf('Montage_%s_%s', genotype, channelName)));
        savefig(fig, [outFileBase '.fig']);
        print(fig, [outFileBase '.png'], '-dpng', sprintf('-r%d', config.figureResolution));
        try
            exportgraphics(fig, [outFileBase '.eps'], 'ContentType', 'vector');
        catch ME
            warning('EPS export failed for %s (PNG/FIG still saved OK): %s', outFileBase, ME.message);
        end
        
        try
            close(fig);
        catch
            % Figure handle already invalid (can happen over remote desktop) - nothing to close
        end
    end
end

fprintf('  Montages created\n');
end

function saveHistogramData(profileData, imageInfo, config)
data.allProfileData.Linear = profileData;
data.imageInfo = imageInfo;
data.binSize_alignment = config.binSize_alignment;
data.binSize_analysis = config.binSize_analysis;
data.normalizationMethods = {'Linear'};

save(char(config.histogramDataFile), 'data', '-v7.3');
fprintf('  Saved: %s\n', config.histogramDataFile);
end

%% ═════════════════════════════════════════════════════════════
%% STEP 3: ALIGNMENT AND STATISTICAL ANALYSIS
%% ═════════════════════════════════════════════════════════════

function step3_AlignmentAndStatistics(config)
% Align using bin 80, analyze at bin 20

fprintf('[1/10] Setting up directories...\n');
setupStep3Directories(config);

fprintf('[2/10] Loading histogram data...\n');
loaded = load(char(config.histogramDataFile));
allProfileData = loaded.data.allProfileData;
imageInfo = loaded.data.imageInfo;

fprintf('[3/10] Aligning profiles (bin %d → bin %d)...\n', ...
        config.binSize_alignment, config.binSize_analysis);
alignedData = alignAllProfiles(allProfileData, config);

fprintf('[4/10] Performing Mann-Whitney tests on per-animal summary statistics...\n');
statsData = performMannWhitneyTests(alignedData, config);

fprintf('[5/10] Creating comparison plots...\n');
createBeforeAfterPlots(alignedData, config);
createWTvsKOComparisonPlots(alignedData, statsData, config);

fprintf('[6/10] Creating QC plots...\n');
createValleyDiagnosticPlots(alignedData, allProfileData, config);
createAlignmentPlots(alignedData, config);

if config.saveShiftedImages
    fprintf('[7/10] Saving shifted images...\n');
    saveShiftedImages(alignedData, imageInfo, config);
end

if config.saveAverageImages
    fprintf('[8/10] Creating average images (mean across animals)...\n');
    createAverageImages(alignedData, imageInfo, config, 'mean');
    fprintf('[8/10] Creating average images (MIP across animals)...\n');
    createAverageImages(alignedData, imageInfo, config, 'max');
    fprintf('[8/10] Creating average images (STD across animals)...\n');
    createAverageImages(alignedData, imageInfo, config, 'std');
end

if config.saveHeatmaps
    fprintf('[9/10] Creating heatmaps (mean across animals)...\n');
    createHeatmaps(alignedData, config, 'mean');
    fprintf('[9/10] Creating heatmaps (MIP across animals)...\n');
    createHeatmaps(alignedData, config, 'max');
    fprintf('[9/10] Creating heatmaps (STD across animals)...\n');
    createHeatmaps(alignedData, config, 'std');
end

saveAlignedDataWithStats(alignedData, statsData, imageInfo, config);

fprintf('[10/10] Generating documentation files...\n');
generateDocumentationFiles(alignedData, statsData, config);

fprintf('✓ Step 3 Complete!\n');
end

function setupStep3Directories(config)
if ~exist(config.alignmentDir, 'dir'), mkdir(config.alignmentDir); end

subdirs = {'Aligned_Profiles', 'Outliers', 'QC_Plots', 'Before_After_Comparison', 'Statistical_Analysis'};
for i = 1:length(subdirs)
    dirPath = fullfile(config.alignmentDir, subdirs{i});
    if ~exist(dirPath, 'dir'), mkdir(dirPath); end
end

if config.saveShiftedImages
    if ~exist(config.shiftedImagesDir, 'dir'), mkdir(config.shiftedImagesDir); end
    
    channels = {'Red', 'Green'};
    genotypes = {'WT', 'KO'};
    for ch = 1:length(channels)
        for g = 1:length(genotypes)
            dirPath = fullfile(config.shiftedImagesDir, channels{ch}, genotypes{g});
            if ~exist(dirPath, 'dir'), mkdir(dirPath); end
        end
    end
end

if config.saveHeatmaps
    if ~exist(config.heatmapsDir, 'dir'), mkdir(config.heatmapsDir); end
end

if config.saveAverageImages
    if ~exist(config.averageImagesDir, 'dir'), mkdir(config.averageImagesDir); end
end
end

function alignedData = alignAllProfiles(allProfileData, config)
% Align using bin 80, apply to bin 20

channels = {'Red', 'Green'};
alignedData = struct();
normMethod = 'Linear';

for ch = 1:length(channels)
    channelName = channels{ch};
    
    % Get bin 80 data for alignment
    WT_profiles_align = allProfileData.(normMethod).(channelName).WT.profiles_align;
    WT_x_coords_align = allProfileData.(normMethod).(channelName).WT.x_coords_align;
    WT_animals = allProfileData.(normMethod).(channelName).WT.animals;
    
    KO_profiles_align = allProfileData.(normMethod).(channelName).KO.profiles_align;
    KO_x_coords_align = allProfileData.(normMethod).(channelName).KO.x_coords_align;
    KO_animals = allProfileData.(normMethod).(channelName).KO.animals;
    
    % Get bin 20 data for analysis
    WT_profiles_analysis = allProfileData.(normMethod).(channelName).WT.profiles_analysis;
    WT_x_coords_analysis = allProfileData.(normMethod).(channelName).WT.x_coords_analysis;
    
    KO_profiles_analysis = allProfileData.(normMethod).(channelName).KO.profiles_analysis;
    KO_x_coords_analysis = allProfileData.(normMethod).(channelName).KO.x_coords_analysis;
    
    % Ensure correct orientation
    WT_profiles_align = ensureColumnProfiles(WT_profiles_align, length(WT_x_coords_align));
    KO_profiles_align = ensureColumnProfiles(KO_profiles_align, length(KO_x_coords_align));
    WT_profiles_analysis = ensureColumnProfiles(WT_profiles_analysis, length(WT_x_coords_analysis));
    KO_profiles_analysis = ensureColumnProfiles(KO_profiles_analysis, length(KO_x_coords_analysis));
    
    % Combine for valley detection
    combined_profiles_align = [WT_profiles_align, KO_profiles_align];
    num_WT = size(WT_profiles_align, 2);
    num_KO = size(KO_profiles_align, 2);
    
    x_coords_align = WT_x_coords_align(:);
    
    % Detect valleys in bin 80 data
    [valley_positions, valley_intensities] = detectValleys(combined_profiles_align, x_coords_align, config);
    
    mean_valley_position = mean(valley_positions);
    mean_valley_intensity = mean(valley_intensities);
    mean_valley_pixel = x_coords_align(round(mean_valley_position));
    
    % Calculate shifts in bin 80 units
    shifts_bin80 = zeros(1, num_WT + num_KO);
    for i = 1:(num_WT + num_KO)
        shifts_bin80(i) = round(mean_valley_position - valley_positions(i));
    end
    
    % Convert to pixels
    shifts_pixels = shifts_bin80 * config.binSize_alignment;
    
    % Convert to bin 20 shifts
    shifts_bin20 = round(shifts_pixels / config.binSize_analysis);
    
    % Apply shifts to bin 20 profiles
    combined_profiles_analysis = [WT_profiles_analysis, KO_profiles_analysis];
    x_coords_analysis = WT_x_coords_analysis(:);
    
    valley_y_offsets = mean_valley_intensity - valley_intensities;
    
    aligned_profiles_analysis = zeros(size(combined_profiles_analysis));
    for i = 1:size(combined_profiles_analysis, 2)
        profile = combined_profiles_analysis(:, i);
        shift = shifts_bin20(i);
        y_offset = valley_y_offsets(i);
        
        % Apply X-shift
        if shift > 0
            shifted = [repmat(profile(1), shift, 1); profile(1:end-shift)];
        elseif shift < 0
            shifted = [profile(-shift+1:end); repmat(profile(end), -shift, 1)];
        else
            shifted = profile;
        end
        
        aligned_profiles_analysis(:, i) = shifted + y_offset;
    end
    
    % Also apply shifts to bin 80 profiles for plotting
    aligned_profiles_bin80 = zeros(size(combined_profiles_align));
    for i = 1:size(combined_profiles_align, 2)
        profile = combined_profiles_align(:, i);
        shift = shifts_bin80(i);
        y_offset = valley_y_offsets(i);
        
        if shift > 0
            shifted = [repmat(profile(1), shift, 1); profile(1:end-shift)];
        elseif shift < 0
            shifted = [profile(-shift+1:end); repmat(profile(end), -shift, 1)];
        else
            shifted = profile;
        end
        
        aligned_profiles_bin80(:, i) = shifted + y_offset;
    end
    
    % Create 0-1 normalized versions
    aligned_profiles_01 = zeros(size(aligned_profiles_analysis));
    for i = 1:size(aligned_profiles_analysis, 2)
        prof = aligned_profiles_analysis(:, i);
        aligned_profiles_01(:, i) = (prof - min(prof)) / (max(prof) - min(prof));
    end
    
    aligned_profiles_01_bin80 = zeros(size(aligned_profiles_bin80));
    for i = 1:size(aligned_profiles_bin80, 2)
        prof = aligned_profiles_bin80(:, i);
        aligned_profiles_01_bin80(:, i) = (prof - min(prof)) / (max(prof) - min(prof));
    end
    
    % Convert to microns (centered at valley = 0)
    [x_microns_bin20, ~] = convertToMicrons(x_coords_analysis, mean_valley_pixel, config.pixelSize);
    [x_microns_bin80, ~] = convertToMicrons(x_coords_align, mean_valley_pixel, config.pixelSize);
    
    % Split into WT and KO - bin 20
    WT_aligned = aligned_profiles_analysis(:, 1:num_WT);
    KO_aligned = aligned_profiles_analysis(:, num_WT+1:end);
    WT_aligned_01 = aligned_profiles_01(:, 1:num_WT);
    KO_aligned_01 = aligned_profiles_01(:, num_WT+1:end);
    
    % Split into WT and KO - bin 80
    WT_aligned_bin80 = aligned_profiles_bin80(:, 1:num_WT);
    KO_aligned_bin80 = aligned_profiles_bin80(:, num_WT+1:end);
    WT_aligned_01_bin80 = aligned_profiles_01_bin80(:, 1:num_WT);
    KO_aligned_01_bin80 = aligned_profiles_01_bin80(:, num_WT+1:end);
    
    % Calculate statistics - bin 20
    WT_mean = mean(WT_aligned, 2);
    WT_std = std(WT_aligned, 0, 2);
    WT_sem = WT_std / sqrt(num_WT);
    WT_mean_01 = mean(WT_aligned_01, 2);
    WT_std_01 = std(WT_aligned_01, 0, 2);
    
    KO_mean = mean(KO_aligned, 2);
    KO_std = std(KO_aligned, 0, 2);
    KO_sem = KO_std / sqrt(num_KO);
    KO_mean_01 = mean(KO_aligned_01, 2);
    KO_std_01 = std(KO_aligned_01, 0, 2);
    
    % Calculate statistics - bin 80
    WT_mean_bin80 = mean(WT_aligned_bin80, 2);
    WT_std_bin80 = std(WT_aligned_bin80, 0, 2);
    WT_mean_01_bin80 = mean(WT_aligned_01_bin80, 2);
    WT_std_01_bin80 = std(WT_aligned_01_bin80, 0, 2);
    
    KO_mean_bin80 = mean(KO_aligned_bin80, 2);
    KO_std_bin80 = std(KO_aligned_bin80, 0, 2);
    KO_mean_01_bin80 = mean(KO_aligned_01_bin80, 2);
    KO_std_01_bin80 = std(KO_aligned_01_bin80, 0, 2);
    
    % Identify outliers
    WT_deviations = zeros(1, num_WT);
    for i = 1:num_WT
        WT_deviations(i) = sqrt(mean((WT_aligned(:, i) - WT_mean).^2));
    end
    WT_outliers = find(WT_deviations > mean(WT_deviations) + config.outlierThreshold * std(WT_deviations));
    
    KO_deviations = zeros(1, num_KO);
    for i = 1:num_KO
        KO_deviations(i) = sqrt(mean((KO_aligned(:, i) - KO_mean).^2));
    end
    KO_outliers = find(KO_deviations > mean(KO_deviations) + config.outlierThreshold * std(KO_deviations));
    
    % Store WT data
    alignedData.(normMethod).(channelName).WT = struct(...
        'original_profiles', WT_profiles_analysis, ...
        'aligned_profiles', WT_aligned, ...
        'aligned_profiles_01', WT_aligned_01, ...
        'aligned_profiles_bin80', WT_aligned_bin80, ...
        'aligned_profiles_01_bin80', WT_aligned_01_bin80, ...
        'x_coords', x_coords_analysis, ...
        'x_coords_bin80', x_coords_align, ...
        'x_microns', x_microns_bin20, ...
        'x_microns_bin80', x_microns_bin80, ...
        'mean_valley_pixel', mean_valley_pixel, ...
        'animals', WT_animals, ...
        'shifts', shifts_bin20(1:num_WT), ...
        'shifts_pixels', shifts_pixels(1:num_WT), ...
        'valley_positions', valley_positions(1:num_WT), ...
        'valley_y_offsets', valley_y_offsets(1:num_WT), ...
        'mean_profile', WT_mean, ...
        'std_profile', WT_std, ...
        'sem_profile', WT_sem, ...
        'mean_profile_01', WT_mean_01, ...
        'std_profile_01', WT_std_01, ...
        'mean_profile_bin80', WT_mean_bin80, ...
        'std_profile_bin80', WT_std_bin80, ...
        'mean_profile_01_bin80', WT_mean_01_bin80, ...
        'std_profile_01_bin80', WT_std_01_bin80, ...
        'deviations', WT_deviations, ...
        'outliers', WT_outliers, ...
        'profiles_align', WT_profiles_align, ...
        'x_coords_align', x_coords_align);
    
    % Store KO data
    alignedData.(normMethod).(channelName).KO = struct(...
        'original_profiles', KO_profiles_analysis, ...
        'aligned_profiles', KO_aligned, ...
        'aligned_profiles_01', KO_aligned_01, ...
        'aligned_profiles_bin80', KO_aligned_bin80, ...
        'aligned_profiles_01_bin80', KO_aligned_01_bin80, ...
        'x_coords', x_coords_analysis, ...
        'x_coords_bin80', x_coords_align, ...
        'x_microns', x_microns_bin20, ...
        'x_microns_bin80', x_microns_bin80, ...
        'mean_valley_pixel', mean_valley_pixel, ...
        'animals', KO_animals, ...
        'shifts', shifts_bin20(num_WT+1:end), ...
        'shifts_pixels', shifts_pixels(num_WT+1:end), ...
        'valley_positions', valley_positions(num_WT+1:end), ...
        'valley_y_offsets', valley_y_offsets(num_WT+1:end), ...
        'mean_profile', KO_mean, ...
        'std_profile', KO_std, ...
        'sem_profile', KO_sem, ...
        'mean_profile_01', KO_mean_01, ...
        'std_profile_01', KO_std_01, ...
        'mean_profile_bin80', KO_mean_bin80, ...
        'std_profile_bin80', KO_std_bin80, ...
        'mean_profile_01_bin80', KO_mean_01_bin80, ...
        'std_profile_01_bin80', KO_std_01_bin80, ...
        'deviations', KO_deviations, ...
        'outliers', KO_outliers, ...
        'profiles_align', KO_profiles_align, ...
        'x_coords_align', x_coords_align);
    
    fprintf('    %s: ✓ (WT: %d profiles, KO: %d profiles)\n', ...
            channelName, num_WT, num_KO);
end

fprintf('  Alignment complete\n');
end

function profiles = ensureColumnProfiles(profiles, expected_length)
if size(profiles, 1) ~= expected_length
    profiles = profiles';
end
end

function [valley_positions, valley_intensities] = detectValleys(profiles, x_coords, config)
num_profiles = size(profiles, 2);
valley_positions = zeros(1, num_profiles);
valley_intensities = zeros(1, num_profiles);

search_start_idx = find(x_coords >= config.valleySearchMin, 1);
search_end_idx = find(x_coords <= config.valleySearchMax, 1, 'last');

for i = 1:num_profiles
    profile = profiles(:, i);
    
    if config.smoothWindow > 1
        profile_smooth = smooth(profile, config.smoothWindow);
    else
        profile_smooth = profile;
    end
    
    [~, min_idx] = min(profile_smooth(search_start_idx:search_end_idx));
    valley_idx = search_start_idx + min_idx - 1;
    
    valley_positions(i) = valley_idx;
    valley_intensities(i) = profile(valley_idx);
end
end

function [x_microns, valley_microns] = convertToMicrons(x_pixels, valley_pixel, pixelSize)
valley_microns = 0;
x_microns = (x_pixels - valley_pixel) * pixelSize;
end

function statsData = performMannWhitneyTests(alignedData, config)
% Perform Mann-Whitney U tests comparing per-animal summary statistics

channels = {'Red', 'Green'};
normMethod = 'Linear';
statsData = struct();

for ch = 1:length(channels)
    channelName = channels{ch};
    
    % Get aligned profiles and x coordinates
    x_microns = alignedData.(normMethod).(channelName).WT.x_microns;
    
    WT_profiles = alignedData.(normMethod).(channelName).WT.aligned_profiles;
    KO_profiles = alignedData.(normMethod).(channelName).KO.aligned_profiles;
    
    WT_profiles_01 = alignedData.(normMethod).(channelName).WT.aligned_profiles_01;
    KO_profiles_01 = alignedData.(normMethod).(channelName).KO.aligned_profiles_01;
    
    num_WT = size(WT_profiles, 2);
    num_KO = size(KO_profiles, 2);
    
    % Find indices for negative (x < 0) and positive (x > 0) halves
    negative_idx = x_microns < 0;
    positive_idx = x_microns > 0;
    
    % ========== ORIGINAL DATA ==========
    % Compute per-animal mean for each half
    WT_neg_means = zeros(num_WT, 1);
    WT_pos_means = zeros(num_WT, 1);
    for i = 1:num_WT
        WT_neg_means(i) = mean(WT_profiles(negative_idx, i));
        WT_pos_means(i) = mean(WT_profiles(positive_idx, i));
    end
    
    KO_neg_means = zeros(num_KO, 1);
    KO_pos_means = zeros(num_KO, 1);
    for i = 1:num_KO
        KO_neg_means(i) = mean(KO_profiles(negative_idx, i));
        KO_pos_means(i) = mean(KO_profiles(positive_idx, i));
    end
    
    % Mann-Whitney tests on original data
    [p_neg_orig, ~, stats_neg_orig] = ranksum(WT_neg_means, KO_neg_means);
    [p_pos_orig, ~, stats_pos_orig] = ranksum(WT_pos_means, KO_pos_means);
    
    % Compute effect size (Cohen's d) and power for original data
    [d_neg_orig, power_neg_orig] = computeEffectSizeAndPower(WT_neg_means, KO_neg_means, config.alpha);
    [d_pos_orig, power_pos_orig] = computeEffectSizeAndPower(WT_pos_means, KO_pos_means, config.alpha);
    
    % ========== NORMALIZED DATA ==========
    % Compute per-animal mean for each half
    WT_neg_means_01 = zeros(num_WT, 1);
    WT_pos_means_01 = zeros(num_WT, 1);
    for i = 1:num_WT
        WT_neg_means_01(i) = mean(WT_profiles_01(negative_idx, i));
        WT_pos_means_01(i) = mean(WT_profiles_01(positive_idx, i));
    end
    
    KO_neg_means_01 = zeros(num_KO, 1);
    KO_pos_means_01 = zeros(num_KO, 1);
    for i = 1:num_KO
        KO_neg_means_01(i) = mean(KO_profiles_01(negative_idx, i));
        KO_pos_means_01(i) = mean(KO_profiles_01(positive_idx, i));
    end
    
    % Mann-Whitney tests on normalized data
    [p_neg_norm, ~, stats_neg_norm] = ranksum(WT_neg_means_01, KO_neg_means_01);
    [p_pos_norm, ~, stats_pos_norm] = ranksum(WT_pos_means_01, KO_pos_means_01);
    
    % Compute effect size and power for normalized data
    [d_neg_norm, power_neg_norm] = computeEffectSizeAndPower(WT_neg_means_01, KO_neg_means_01, config.alpha);
    [d_pos_norm, power_pos_norm] = computeEffectSizeAndPower(WT_pos_means_01, KO_pos_means_01, config.alpha);
    
    % Store original data results
    statsData.(normMethod).(channelName).original = struct(...
        'p_negative', p_neg_orig, ...
        'p_positive', p_pos_orig, ...
        'U_negative', stats_neg_orig.ranksum, ...
        'U_positive', stats_pos_orig.ranksum, ...
        'cohens_d_negative', d_neg_orig, ...
        'cohens_d_positive', d_pos_orig, ...
        'power_negative', power_neg_orig, ...
        'power_positive', power_pos_orig, ...
        'WT_negative_mean', mean(WT_neg_means), ...
        'WT_negative_std', std(WT_neg_means), ...
        'KO_negative_mean', mean(KO_neg_means), ...
        'KO_negative_std', std(KO_neg_means), ...
        'WT_positive_mean', mean(WT_pos_means), ...
        'WT_positive_std', std(WT_pos_means), ...
        'KO_positive_mean', mean(KO_pos_means), ...
        'KO_positive_std', std(KO_pos_means), ...
        'WT_negative_values', WT_neg_means, ...
        'KO_negative_values', KO_neg_means, ...
        'WT_positive_values', WT_pos_means, ...
        'KO_positive_values', KO_pos_means, ...
        'num_WT', num_WT, ...
        'num_KO', num_KO);
    
    % Store normalized data results
    statsData.(normMethod).(channelName).normalized = struct(...
        'p_negative', p_neg_norm, ...
        'p_positive', p_pos_norm, ...
        'U_negative', stats_neg_norm.ranksum, ...
        'U_positive', stats_pos_norm.ranksum, ...
        'cohens_d_negative', d_neg_norm, ...
        'cohens_d_positive', d_pos_norm, ...
        'power_negative', power_neg_norm, ...
        'power_positive', power_pos_norm, ...
        'WT_negative_mean', mean(WT_neg_means_01), ...
        'WT_negative_std', std(WT_neg_means_01), ...
        'KO_negative_mean', mean(KO_neg_means_01), ...
        'KO_negative_std', std(KO_neg_means_01), ...
        'WT_positive_mean', mean(WT_pos_means_01), ...
        'WT_positive_std', std(WT_pos_means_01), ...
        'KO_positive_mean', mean(KO_pos_means_01), ...
        'KO_positive_std', std(KO_pos_means_01), ...
        'WT_negative_values', WT_neg_means_01, ...
        'KO_negative_values', KO_neg_means_01, ...
        'WT_positive_values', WT_pos_means_01, ...
        'KO_positive_values', KO_pos_means_01, ...
        'num_WT', num_WT, ...
        'num_KO', num_KO);
    
    % Print results
    fprintf('    %s Original:\n', channelName);
    fprintf('      Negative half: p=%.4f, d=%.2f, power=%.2f\n', p_neg_orig, d_neg_orig, power_neg_orig);
    fprintf('      Positive half: p=%.4f, d=%.2f, power=%.2f\n', p_pos_orig, d_pos_orig, power_pos_orig);
    fprintf('    %s Normalized:\n', channelName);
    fprintf('      Negative half: p=%.4f, d=%.2f, power=%.2f\n', p_neg_norm, d_neg_norm, power_neg_norm);
    fprintf('      Positive half: p=%.4f, d=%.2f, power=%.2f\n', p_pos_norm, d_pos_norm, power_pos_norm);
end

fprintf('  Mann-Whitney tests complete\n');
end

function [cohens_d, power] = computeEffectSizeAndPower(group1, group2, alpha)
% Compute Cohen's d effect size
n1 = length(group1);
n2 = length(group2);
mean1 = mean(group1);
mean2 = mean(group2);
std1 = std(group1);
std2 = std(group2);

% Pooled standard deviation
pooled_std = sqrt(((n1-1)*std1^2 + (n2-1)*std2^2) / (n1 + n2 - 2));

% Cohen's d
if pooled_std > 0
    cohens_d = abs(mean1 - mean2) / pooled_std;
else
    cohens_d = 0;
end

% Compute post-hoc power using non-central t-distribution
df = n1 + n2 - 2;
ncp = cohens_d * sqrt(n1 * n2 / (n1 + n2));  % Non-centrality parameter
t_crit = tinv(1 - alpha/2, df);  % Critical t-value for two-tailed test

% Power = P(reject H0 | H1 is true)
if ncp > 0
    power = 1 - nctcdf(t_crit, df, ncp) + nctcdf(-t_crit, df, ncp);
else
    power = alpha;  % If no effect, power equals alpha
end

% Ensure power is in valid range
power = max(0, min(1, power));
end

%% ═════════════════════════════════════════════════════════════
%% AVERAGE IMAGES AND HEATMAPS
%% ═════════════════════════════════════════════════════════════

function createAverageImages(alignedData, imageInfo, config, projMode)
% Load images and combine across aligned animals per genotype/channel.
% projMode selects the across-animal projection:
%   'mean' (default) - pixel-wise mean across aligned animals (original behavior)
%   'max'            - pixel-wise maximum-intensity projection across animals
%   'std'            - pixel-wise standard deviation across animals
%                      (like ImageJ Z-project -> Standard Deviation)
% Each variant writes to its own output folder and .mat, so all sets coexist.
% NOTE: images are NOT X-shifted (shifting disabled); the aligned-set membership
% still selects which animals contribute.
if nargin < 4 || isempty(projMode); projMode = 'mean'; end

switch lower(projMode)
    case 'mean'
        outDir   = config.averageImagesDir;
        projTag  = '';            % filenames unchanged for backward compatibility
        projWord = 'Mean';
    case 'max'
        outDir   = config.averageImagesMipDir;
        projTag  = '_MIP';
        projWord = 'MIP';
    case 'std'
        outDir   = config.averageImagesStdDir;
        projTag  = '_STD';
        projWord = 'STD';
    otherwise
        error('createAverageImages: unknown projMode "%s"', projMode);
end
if ~exist(char(outDir), 'dir'), mkdir(char(outDir)); end

if ~exist(char(config.normalizedDataFile), 'file')
    warning('normalized_data.mat not found! Skipping average images.');
    return;
end

loaded = load(char(config.normalizedDataFile));
linearData = loaded.data.linearData;

channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};
normMethod = 'Linear';

% First pass: compute all projected images
for ch = 1:length(channels)
    channelName = channels{ch};
    original_images = linearData.(channelName).images;
    [H, W, ~] = size(original_images);
    
    for g = 1:length(genotypes)
        genotype = genotypes{g};
        
        shifts_pixels = alignedData.(normMethod).(channelName).(genotype).shifts_pixels;
        animals = alignedData.(normMethod).(channelName).(genotype).animals;
        
        genotype_idx = find(strcmp({imageInfo.genotype}, genotype));
        
        % Collect shifted images
        shifted_stack = zeros(H, W, length(genotype_idx));
        count = 0;
        
        for i = 1:length(genotype_idx)
            img_idx = genotype_idx(i);
            animal_id = imageInfo(img_idx).animal;
            
            animal_shift_idx = find(animals == animal_id);
            if isempty(animal_shift_idx), continue; end
            
            % NOTE: image shifting is intentionally disabled. The bin-80 valley
            % shift is quantized to whole 80-px bins, which over-shifts when
            % little/no shift is warranted, and the edge-replicate padding it
            % used smears a streak on the far side. Profiles/stats still use the
            % computed alignment (unchanged); only the IMAGES are left unshifted.
            original_img = original_images(:, :, img_idx);
            shifted_img = original_img;   % no shift (was applyXShiftWithPadding)
            
            count = count + 1;
            shifted_stack(:, :, count) = shifted_img;
        end
        
        % Project across aligned animals (mean, max, or std)
        shifted_stack = shifted_stack(:, :, 1:count);
        if strcmpi(projMode, 'max')
            avg_image = max(shifted_stack, [], 3);
        elseif strcmpi(projMode, 'std')
            % Standard deviation across animals (ImageJ "Standard Deviation"
            % Z-project). Uses the sample SD (N-1) via std's default flag 0.
            % With a single contributing animal SD is 0 everywhere.
            if count > 1
                avg_image = std(shifted_stack, 0, 3);
            else
                avg_image = zeros(size(shifted_stack, 1), size(shifted_stack, 2));
            end
        else
            avg_image = mean(shifted_stack, 3);
        end
        
        % Save projected image as 16-bit TIFF (original values)
        avgFilename = sprintf('Average_%s_%s%s.tif', genotype, channelName, projTag);
        avgPath = char(fullfile(outDir, avgFilename));
        imwrite(uint16(avg_image), avgPath);
        
        % Store for later processing
        avgData.(channelName).(genotype) = avg_image;
        
        fprintf('    %s %s: %s of %d images computed\n', genotype, channelName, projWord, count);
    end
end

% Compute global intensity scale across ALL 4 average images
allImages = [avgData.Red.WT(:); avgData.Red.KO(:); avgData.Green.WT(:); avgData.Green.KO(:)];
globalMin = min(allImages);
globalMax = max(allImages);
globalScale.min = globalMin;
globalScale.max = globalMax;

fprintf('    Global intensity scale: [%.1f, %.1f]\n', globalMin, globalMax);

% Second pass: save 8-bit versions with both individual and global scaling
for ch = 1:length(channels)
    channelName = channels{ch};
    
    for g = 1:length(genotypes)
        genotype = genotypes{g};
        avg_image = avgData.(channelName).(genotype);
        
        % Save 8-bit with individual scaling (for reference)
        avg_image_8bit_individual = uint8(255 * mat2gray(avg_image));
        avgFilename8_ind = sprintf('Average_%s_%s%s_8bit.tif', genotype, channelName, projTag);
        avgPath8_ind = char(fullfile(outDir, avgFilename8_ind));
        imwrite(avg_image_8bit_individual, avgPath8_ind);
        
        % Save 8-bit with GLOBAL scaling (for consistent comparison)
        avg_image_normalized = (avg_image - globalMin) / (globalMax - globalMin);
        avg_image_normalized = max(0, min(1, avg_image_normalized));  % Clamp to [0, 1]
        avg_image_8bit_global = uint8(255 * avg_image_normalized);
        avgFilename8_global = sprintf('Average_%s_%s%s_8bit_globalscale.tif', genotype, channelName, projTag);
        avgPath8_global = char(fullfile(outDir, avgFilename8_global));
        imwrite(avg_image_8bit_global, avgPath8_global);
        
        fprintf('    %s %s: 8-bit versions saved (individual and global scale)\n', genotype, channelName);
    end
end

% Save projected data and global scale for heatmaps
save(char(fullfile(outDir, 'average_images_data.mat')), 'avgData', 'globalScale', '-v7.3');

fprintf('  %s images saved with global intensity scale\n', projWord);
end

function createHeatmaps(alignedData, config, projMode)
% Create heatmaps showing intensity differences between WT and KO
% Simplified to key colormaps only:
%   - Greater plots: plasma, white-red/white-blue, white-blue/white-red
%   - Combined diverging: blue-white-red, red-white-blue
% Fixed color scales: 0-150 for greater plots, -150 to +150 for combined
% Global intensity scale applied to all average images for consistent comparison
%
% projMode selects which projected average images to build heatmaps from:
%   'mean' (default) - the mean-across-animals averages (05_Heatmaps)
%   'max'            - the MIP-across-animals images     (05b_Heatmaps_MIP)
%   'std'            - the STD-across-animals images      (05d_Heatmaps_STD)
if nargin < 3 || isempty(projMode); projMode = 'mean'; end

switch lower(projMode)
    case 'mean'
        avgDir  = config.averageImagesDir;
        heatDir = config.heatmapsDir;
    case 'max'
        avgDir  = config.averageImagesMipDir;
        heatDir = config.heatmapsMipDir;
    case 'std'
        avgDir  = config.averageImagesStdDir;
        heatDir = config.heatmapsStdDir;
    otherwise
        error('createHeatmaps: unknown projMode "%s"', projMode);
end
if ~exist(char(heatDir), 'dir'), mkdir(char(heatDir)); end

% Load projected images (mean or MIP across animals)
avgDataFile = fullfile(avgDir, 'average_images_data.mat');
if ~exist(char(avgDataFile), 'file')
    warning('Average images data not found! Skipping heatmaps.');
    return;
end

loaded = load(char(avgDataFile));
avgData = loaded.avgData;
globalScale = loaded.globalScale;  % Load global intensity scale

fprintf('    Using global intensity scale: [%.1f, %.1f]\n', globalScale.min, globalScale.max);

% Label for the quantity being shown, so STD figures aren't mislabeled as mean.
switch lower(projMode)
    case 'max';  projLabel = 'Max Intensity';
    case 'std';  projLabel = 'Std Dev';
    otherwise;   projLabel = 'Mean Intensity';
end

% Get midline position from alignedData
mean_valley_pixel = alignedData.Linear.Red.WT.mean_valley_pixel;

channels = {'Red', 'Green'};

% X-axis tick positions
x_ticks = [-300, -200, -100, 0, 100, 200, 300];

% High resolution for average images
avgImageResolution = 900;

% Difference color-scale limits. Mean/max keep the established fixed 0-150
% (and +/-150) scale. STD values are much smaller than mean intensities, so a
% fixed 150 ceiling would render the STD difference plots nearly black; instead
% derive the STD scale from the actual WT/KO STD difference across both
% channels so the plasma/diverging maps use their full color range.
if strcmpi(projMode, 'std')
    maxAbsDiff = 0;
    for cc = 1:length(channels)
        dW = avgData.(channels{cc}).WT - avgData.(channels{cc}).KO;
        maxAbsDiff = max(maxAbsDiff, max(abs(dW(:))));
    end
    if ~isfinite(maxAbsDiff) || maxAbsDiff <= 0
        maxAbsDiff = 1;
    end
    % Round up to a tidy number for readable colorbar ticks
    diff_scale_max = ceil(maxAbsDiff);
    combined_scale = diff_scale_max;
    fprintf('    STD difference color scale set to +/-%g (data-driven)\n', diff_scale_max);
else
    diff_scale_max = 150;  % For "greater" plots: 0 to 150
    combined_scale = 150;  % For combined plots: -150 to +150
end

% Global intensity limits for grayscale images
intensityLimits = [globalScale.min, globalScale.max];

for ch = 1:length(channels)
    channelName = channels{ch};
    
    % Clean up any existing figures and force memory refresh
    drawnow;
    pause(0.1);
    
    fprintf('    Processing %s channel heatmaps...\n', channelName);
    
    WT_avg = avgData.(channelName).WT;
    KO_avg = avgData.(channelName).KO;
    [H, W] = size(WT_avg);
    
    % Flip images vertically for correct orientation (top of tissue at top of figure)
    WT_avg_flip = flipud(WT_avg);
    KO_avg_flip = flipud(KO_avg);
    
    % Compute difference maps on flipped images
    diff_raw = WT_avg_flip - KO_avg_flip;  % Positive = WT > KO
    diff_WT_greater = max(diff_raw, 0);  % Only where WT > KO
    diff_KO_greater = max(-diff_raw, 0);  % Only where KO > WT
    
    % Create x coordinates centered at midline (0 at aligned valley)
    x_microns = ((0:(W-1)) - mean_valley_pixel) * config.pixelSize;
    
    % Create y coordinates from 0 at bottom to max at top
    y_microns = (0:(H-1)) * config.pixelSize;
    
    % ===== SAVE DIFFERENCE MAPS AS HIGH-RES RGB TIF FILES =====
    % Flip vertically for correct orientation and apply colormaps
    % These can be used to manually replace compressed images in EPS files
    
    % Flip images for correct orientation in saved TIFs
    diff_WT_greater_forTIF = flipud(diff_WT_greater);
    diff_KO_greater_forTIF = flipud(diff_KO_greater);
    diff_raw_forTIF = flipud(diff_raw);
    
    % Get colormaps
    cmap_plasma = plasma(256);
    cmap_white_red = white_to_red_cmap(256);
    cmap_white_blue = white_to_blue_cmap(256);
    cmap_bwr = blue_white_red(256);
    cmap_rwb = red_white_blue(256);
    
    % Helper function to apply colormap to data and create RGB image
    % For "greater" maps: scale 0 to diff_scale_max (150 for mean/max,
    % data-driven for std). For combined maps: -combined_scale to +combined_scale.
    
    % === WT > KO with plasma ===
    rgb_WT_plasma = applyColormapToImage(diff_WT_greater_forTIF, cmap_plasma, 0, diff_scale_max);
    imwrite(rgb_WT_plasma, char(fullfile(heatDir, sprintf('DifferenceMap_%s_WT_greater_plasma.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === KO > WT with plasma ===
    rgb_KO_plasma = applyColormapToImage(diff_KO_greater_forTIF, cmap_plasma, 0, diff_scale_max);
    imwrite(rgb_KO_plasma, char(fullfile(heatDir, sprintf('DifferenceMap_%s_KO_greater_plasma.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === WT > KO with white-red ===
    rgb_WT_whitered = applyColormapToImage(diff_WT_greater_forTIF, cmap_white_red, 0, diff_scale_max);
    imwrite(rgb_WT_whitered, char(fullfile(heatDir, sprintf('DifferenceMap_%s_WT_greater_white_red.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === KO > WT with white-blue ===
    rgb_KO_whiteblue = applyColormapToImage(diff_KO_greater_forTIF, cmap_white_blue, 0, diff_scale_max);
    imwrite(rgb_KO_whiteblue, char(fullfile(heatDir, sprintf('DifferenceMap_%s_KO_greater_white_blue.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === WT > KO with white-blue (reverse scheme) ===
    rgb_WT_whiteblue = applyColormapToImage(diff_WT_greater_forTIF, cmap_white_blue, 0, diff_scale_max);
    imwrite(rgb_WT_whiteblue, char(fullfile(heatDir, sprintf('DifferenceMap_%s_WT_greater_white_blue.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === KO > WT with white-red (reverse scheme) ===
    rgb_KO_whitered = applyColormapToImage(diff_KO_greater_forTIF, cmap_white_red, 0, diff_scale_max);
    imwrite(rgb_KO_whitered, char(fullfile(heatDir, sprintf('DifferenceMap_%s_KO_greater_white_red.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === Combined blue-white-red ===
    rgb_combined_bwr = applyColormapToImage(diff_raw_forTIF, cmap_bwr, -combined_scale, combined_scale);
    imwrite(rgb_combined_bwr, char(fullfile(heatDir, sprintf('DifferenceMap_%s_combined_blue_white_red.tif', channelName))), 'tif', 'Compression', 'none');
    
    % === Combined red-white-blue ===
    rgb_combined_rwb = applyColormapToImage(diff_raw_forTIF, cmap_rwb, -combined_scale, combined_scale);
    imwrite(rgb_combined_rwb, char(fullfile(heatDir, sprintf('DifferenceMap_%s_combined_red_white_blue.tif', channelName))), 'tif', 'Compression', 'none');
    
    % Also save MAT file with raw data for reference
    diffMatFile = fullfile(heatDir, sprintf('DifferenceMap_%s_data.mat', channelName));
    save(char(diffMatFile), 'diff_raw', 'diff_WT_greater', 'diff_KO_greater', 'x_microns', 'y_microns', ...
         'diff_scale_max', 'combined_scale');
    
    fprintf('    %s: Saved %d RGB difference map TIFs and MAT file\n', channelName, 8);
    
    % ===== 1. GREATER PLASMA =====
    cmap_plasma = plasma(256);
    
    fig1 = figure('Position', [100, 100, 1500, 500], 'Visible', config.showFigures);
    set(fig1, 'DefaultAxesFontName', 'Arial');
    set(fig1, 'DefaultTextFontName', 'Arial');
    
    subplot(1,3,1);
    imagesc(x_microns, y_microns, WT_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb1 = colorbar; cb1.Label.String = projLabel; cb1.Label.FontName = 'Arial';
    title(sprintf('WT %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,2);
    imagesc(x_microns, y_microns, KO_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb2 = colorbar; cb2.Label.String = projLabel; cb2.Label.FontName = 'Arial';
    title(sprintf('KO %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,3);
    imagesc(x_microns, y_microns, diff_WT_greater);
    axis image; colormap(gca, cmap_plasma); 
    set(gca, 'YDir', 'normal');
    caxis([0, diff_scale_max]);
    cb3 = colorbar; cb3.Label.String = 'WT − KO'; cb3.Label.FontName = 'Arial';
    title(sprintf('%s: WT > KO - %s (plasma)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig1, heatDir, ...
                         sprintf('Heatmap_%s_WT_greater_plasma', channelName), config, avgImageResolution);
    
    % Clear memory before next figure
    drawnow;
    
    % KO > WT plasma
    fig2 = figure('Position', [100, 100, 1500, 500], 'Visible', config.showFigures);
    set(fig2, 'DefaultAxesFontName', 'Arial');
    set(fig2, 'DefaultTextFontName', 'Arial');
    
    ax1 = subplot(1,3,1);
    imagesc(x_microns, y_microns, WT_avg_flip);
    axis image; colormap(ax1, gray); caxis(intensityLimits);
    set(ax1, 'YDir', 'normal');
    cb1 = colorbar; cb1.Label.String = projLabel; cb1.Label.FontName = 'Arial';
    title(sprintf('WT %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(ax1, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(ax1, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    drawnow;
    
    ax2 = subplot(1,3,2);
    imagesc(x_microns, y_microns, KO_avg_flip);
    axis image; colormap(ax2, gray); caxis(intensityLimits);
    set(ax2, 'YDir', 'normal');
    cb2 = colorbar; cb2.Label.String = projLabel; cb2.Label.FontName = 'Arial';
    title(sprintf('KO %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(ax2, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(ax2, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    drawnow;
    
    ax3 = subplot(1,3,3);
    imagesc(x_microns, y_microns, diff_KO_greater);
    axis image; colormap(ax3, cmap_plasma); 
    set(ax3, 'YDir', 'normal');
    caxis([0, diff_scale_max]);
    cb3 = colorbar; cb3.Label.String = 'KO − WT'; cb3.Label.FontName = 'Arial';
    title(sprintf('%s: KO > WT - %s (plasma)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(ax3, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(ax3, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    drawnow;
    
    % Force complete render before save
    pause(0.3);
    
    saveFigureMultiFormatHighRes(fig2, heatDir, ...
                         sprintf('Heatmap_%s_KO_greater_plasma', channelName), config, avgImageResolution);
    
    % ===== 2. GREATER WHITE-RED / WHITE-BLUE (WT>KO red, KO>WT blue) =====
    cmap_white_red = white_to_red_cmap(256);
    cmap_white_blue = white_to_blue_cmap(256);
    
    fig3 = figure('Position', [100, 100, 1500, 500], 'Visible', config.showFigures);
    set(fig3, 'DefaultAxesFontName', 'Arial');
    set(fig3, 'DefaultTextFontName', 'Arial');
    
    subplot(1,3,1);
    imagesc(x_microns, y_microns, WT_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb1 = colorbar; cb1.Label.String = projLabel; cb1.Label.FontName = 'Arial';
    title(sprintf('WT %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,2);
    imagesc(x_microns, y_microns, KO_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb2 = colorbar; cb2.Label.String = projLabel; cb2.Label.FontName = 'Arial';
    title(sprintf('KO %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,3);
    imagesc(x_microns, y_microns, diff_WT_greater);
    axis image; colormap(gca, cmap_white_red); 
    set(gca, 'YDir', 'normal');
    caxis([0, diff_scale_max]);
    cb3 = colorbar; cb3.Label.String = 'WT − KO'; cb3.Label.FontName = 'Arial';
    title(sprintf('%s: WT > KO - %s (white-red)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'k--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig3, heatDir, ...
                         sprintf('Heatmap_%s_WT_greater_white_red', channelName), config, avgImageResolution);
    
    fig4 = figure('Position', [100, 100, 1500, 500], 'Visible', config.showFigures);
    set(fig4, 'DefaultAxesFontName', 'Arial');
    set(fig4, 'DefaultTextFontName', 'Arial');
    
    subplot(1,3,1);
    imagesc(x_microns, y_microns, WT_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb1 = colorbar; cb1.Label.String = projLabel; cb1.Label.FontName = 'Arial';
    title(sprintf('WT %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,2);
    imagesc(x_microns, y_microns, KO_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb2 = colorbar; cb2.Label.String = projLabel; cb2.Label.FontName = 'Arial';
    title(sprintf('KO %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,3);
    imagesc(x_microns, y_microns, diff_KO_greater);
    axis image; colormap(gca, cmap_white_blue); 
    set(gca, 'YDir', 'normal');
    caxis([0, diff_scale_max]);
    cb3 = colorbar; cb3.Label.String = 'KO − WT'; cb3.Label.FontName = 'Arial';
    title(sprintf('%s: KO > WT - %s (white-blue)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'k--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig4, heatDir, ...
                         sprintf('Heatmap_%s_KO_greater_white_blue', channelName), config, avgImageResolution);
    
    % ===== 3. GREATER WHITE-BLUE / WHITE-RED (reverse: WT>KO blue, KO>WT red) =====
    fig5 = figure('Position', [100, 100, 1500, 500], 'Visible', config.showFigures);
    set(fig5, 'DefaultAxesFontName', 'Arial');
    set(fig5, 'DefaultTextFontName', 'Arial');
    
    subplot(1,3,1);
    imagesc(x_microns, y_microns, WT_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb1 = colorbar; cb1.Label.String = projLabel; cb1.Label.FontName = 'Arial';
    title(sprintf('WT %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,2);
    imagesc(x_microns, y_microns, KO_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb2 = colorbar; cb2.Label.String = projLabel; cb2.Label.FontName = 'Arial';
    title(sprintf('KO %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,3);
    imagesc(x_microns, y_microns, diff_WT_greater);
    axis image; colormap(gca, cmap_white_blue); 
    set(gca, 'YDir', 'normal');
    caxis([0, diff_scale_max]);
    cb3 = colorbar; cb3.Label.String = 'WT − KO'; cb3.Label.FontName = 'Arial';
    title(sprintf('%s: WT > KO - %s (white-blue)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'k--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig5, heatDir, ...
                         sprintf('Heatmap_%s_WT_greater_white_blue', channelName), config, avgImageResolution);
    
    fig6 = figure('Position', [100, 100, 1500, 500], 'Visible', config.showFigures);
    set(fig6, 'DefaultAxesFontName', 'Arial');
    set(fig6, 'DefaultTextFontName', 'Arial');
    
    subplot(1,3,1);
    imagesc(x_microns, y_microns, WT_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb1 = colorbar; cb1.Label.String = projLabel; cb1.Label.FontName = 'Arial';
    title(sprintf('WT %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,2);
    imagesc(x_microns, y_microns, KO_avg_flip);
    axis image; colormap(gca, gray); caxis(intensityLimits);
    set(gca, 'YDir', 'normal');
    cb2 = colorbar; cb2.Label.String = projLabel; cb2.Label.FontName = 'Arial';
    title(sprintf('KO %s - %s', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'w--', 'LineWidth', 1.5); hold off;
    
    subplot(1,3,3);
    imagesc(x_microns, y_microns, diff_KO_greater);
    axis image; colormap(gca, cmap_white_red); 
    set(gca, 'YDir', 'normal');
    caxis([0, diff_scale_max]);
    cb3 = colorbar; cb3.Label.String = 'KO − WT'; cb3.Label.FontName = 'Arial';
    title(sprintf('%s: KO > WT - %s (white-red)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'k--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig6, heatDir, ...
                         sprintf('Heatmap_%s_KO_greater_white_red', channelName), config, avgImageResolution);
    
    % ===== 4. COMBINED DIVERGING BLUE-WHITE-RED =====
    cmap_bwr = blue_white_red(256);
    
    fig7 = figure('Position', [100, 100, 800, 800], 'Visible', config.showFigures);
    set(fig7, 'DefaultAxesFontName', 'Arial');
    set(fig7, 'DefaultTextFontName', 'Arial');
    
    imagesc(x_microns, y_microns, diff_raw);
    axis image;
    set(gca, 'YDir', 'normal');
    colormap(cmap_bwr);
    caxis([-combined_scale, combined_scale]);
    cb = colorbar;
    cb.Label.String = 'WT − KO';
    cb.Label.FontName = 'Arial';
    title(sprintf('%s Difference Map - %s\nBlue (KO>WT) - White (0) - Red (WT>KO)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'k--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig7, heatDir, ...
                         sprintf('Heatmap_%s_combined_blue_white_red', channelName), config, avgImageResolution);
    
    % ===== 5. COMBINED DIVERGING RED-WHITE-BLUE (reverse) =====
    cmap_rwb = red_white_blue(256);
    
    fig8 = figure('Position', [100, 100, 800, 800], 'Visible', config.showFigures);
    set(fig8, 'DefaultAxesFontName', 'Arial');
    set(fig8, 'DefaultTextFontName', 'Arial');
    
    imagesc(x_microns, y_microns, diff_raw);
    axis image;
    set(gca, 'YDir', 'normal');
    colormap(cmap_rwb);
    caxis([-combined_scale, combined_scale]);
    cb = colorbar;
    cb.Label.String = 'WT − KO';
    cb.Label.FontName = 'Arial';
    title(sprintf('%s Difference Map - %s\nRed (KO>WT) - White (0) - Blue (WT>KO)', projLabel, channelName), 'FontSize', 12, 'FontName', 'Arial');
    xlabel('Distance (µm)', 'FontName', 'Arial');
    ylabel('Distance (µm)', 'FontName', 'Arial');
    set(gca, 'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1);
    set(gca, 'XTick', x_ticks);
    hold on; xline(0, 'k--', 'LineWidth', 1.5); hold off;
    
    saveFigureMultiFormatHighRes(fig8, heatDir, ...
                         sprintf('Heatmap_%s_combined_red_white_blue', channelName), config, avgImageResolution);
    
    % Print summary
    fprintf('    %s: Max WT>KO=%.2f, Max KO>WT=%.2f (scale: 0-%d for greater, ±%d for combined)\n', ...
            channelName, max(diff_WT_greater(:)), max(diff_KO_greater(:)), diff_scale_max, combined_scale);
end

fprintf('  Heatmaps created (fixed scales: 0-150, ±150)\n');
end

function saveFigureMultiFormatHighRes(fig, outputDir, baseName, config, resolution)
% Save figure with specified resolution for high-quality heatmaps
% Saves both EPS and PDF for flexibility
outputDir = char(outputDir);
baseName = char(baseName);

% Verify figure is valid
if ~ishandle(fig) || ~isgraphics(fig, 'figure')
    warning('Invalid figure handle for %s, skipping save.', baseName);
    return;
end

% Ensure output directory exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Ensure figure is fully rendered before saving
drawnow;
pause(0.2);

if config.savePNG
    pngFile = fullfile(outputDir, [baseName '.png']);
    print(fig, pngFile, '-dpng', sprintf('-r%d', resolution));

    % Validate the saved PNG isn't silently corrupted (blank/black) - this can
    % happen over Remote Desktop without print() throwing any error at all.
    try
        checkImg = imread(pngFile);
        if std(double(checkImg(:))) < 1
            warning('%s looks blank after saving - retrying once...', pngFile);
            drawnow;
            pause(0.5);
            print(fig, pngFile, '-dpng', sprintf('-r%d', resolution));
            checkImg = imread(pngFile);
            if std(double(checkImg(:))) < 1
                warning('%s STILL blank after retry - figure may need to be regenerated later.', pngFile);
            end
        end
    catch
        % Couldn't read back the file to validate it - not fatal, just skip the check
    end
end

if config.saveFIG
    figFile = fullfile(outputDir, [baseName '.fig']);
    savefig(fig, figFile);
end

% Save as both EPS and PDF
if config.saveEPS
    % Force figure refresh (painters is more stable than opengl over Remote Desktop)
    set(fig, 'Renderer', 'painters');
    drawnow;
    pause(0.1);
    
    % Save PDF (more reliable for embedded images)
    pdfFile = fullfile(outputDir, [baseName '.pdf']);
    try
        exportgraphics(fig, pdfFile, 'ContentType', 'auto', 'Resolution', resolution, ...
                      'BackgroundColor', 'white', 'ColorSpace', 'rgb');
    catch
        try
            exportgraphics(fig, pdfFile, 'ContentType', 'auto', 'Resolution', resolution, ...
                          'BackgroundColor', 'white');
        catch
            print(fig, pdfFile, '-dpdf', sprintf('-r%d', resolution));
        end
    end
    
    % Save EPS (for compatibility)
    epsFile = fullfile(outputDir, [baseName '.eps']);
    try
        exportgraphics(fig, epsFile, 'ContentType', 'auto', 'Resolution', resolution, ...
                      'BackgroundColor', 'white', 'ColorSpace', 'rgb');
    catch
        try
            exportgraphics(fig, epsFile, 'ContentType', 'auto', 'Resolution', resolution, ...
                          'BackgroundColor', 'white');
        catch
            try
                print(fig, epsFile, '-depsc2', sprintf('-r%d', resolution));
            catch ME
                warning('EPS export failed for %s (PNG/FIG still saved OK): %s', epsFile, ME.message);
            end
        end
    end
end

if ~config.showFigures || config.closeFIG
    try
        close(fig);
    catch
        % Figure handle already invalid (can happen over remote desktop) - nothing to close
    end
end
end

function cmap = getColormap(name, n)
% Return colormap by name
switch lower(name)
    case 'viridis'
        cmap = viridis(n);
    case 'plasma'
        cmap = plasma(n);
    case 'jet'
        cmap = jet(n);
    otherwise
        cmap = viridis(n);
end
end

function setMicronAxes(imageSize, pixelSize)
H = imageSize(1);
W = imageSize(2);
xlim([0, W * pixelSize]);
ylim([0, H * pixelSize]);
xlabel('Distance (µm)', 'FontName', 'Arial');
ylabel('Distance (µm)', 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
end

function addScaleBar(pixelSize, scaleBarMicrons, imageSize)
H = imageSize(1);
W = imageSize(2);

xlim_current = xlim;
ylim_current = ylim;
padding = 0.05;
x_range = xlim_current(2) - xlim_current(1);
y_range = ylim_current(2) - ylim_current(1);

bar_x_start = xlim_current(2) - (padding + 0.15) * x_range;
bar_x_end = bar_x_start + scaleBarMicrons;
bar_y = ylim_current(2) - padding * y_range;

hold on;
line([bar_x_start, bar_x_end], [bar_y, bar_y], 'Color', 'black', 'LineWidth', 6);
line([bar_x_start, bar_x_end], [bar_y, bar_y], 'Color', 'white', 'LineWidth', 4);

text_x = (bar_x_start + bar_x_end) / 2;
text_y = bar_y - 0.03 * y_range;
text(text_x, text_y, sprintf('%d µm', scaleBarMicrons), ...
    'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
hold off;
end

function generateDocumentationFiles(alignedData, statsData, config)
% Generate documentation text files: figure legends, methods summaries

docDir = fullfile(config.alignmentDir, 'Documentation');
if ~exist(docDir, 'dir'), mkdir(docDir); end

normMethod = 'Linear';

% Get sample sizes
num_WT = size(alignedData.(normMethod).Red.WT.aligned_profiles, 2);
num_KO = size(alignedData.(normMethod).Red.KO.aligned_profiles, 2);

% Get alignment info for documentation
mean_valley_WT = alignedData.(normMethod).Red.WT.mean_valley_pixel;
mean_valley_KO = alignedData.(normMethod).Red.KO.mean_valley_pixel;

% ===== 1. FIGURE LEGENDS =====
legendFile = fullfile(docDir, 'Figure_Legends.txt');
fid = fopen(char(legendFile), 'w');

fprintf(fid, 'FIGURE LEGENDS\n');
fprintf(fid, '==============\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, '--- STEP 1: LINEAR STRETCH MONTAGES ---\n\n');

fprintf(fid, 'Figure: Montage_WT_Red_LinearStretch / Montage_WT_Green_LinearStretch\n');
fprintf(fid, 'Montage of wild-type (WT) images after linear histogram stretching (1%%-99%% saturation). ');
fprintf(fid, 'Each panel shows an individual animal labeled by number. ');
fprintf(fid, 'Linear stretching was applied using group-specific intensity limits to normalize within-group variation while preserving relative intensity differences. ');
fprintf(fid, 'N = %d WT animals.\n\n', num_WT);

fprintf(fid, 'Figure: Montage_KO_Red_LinearStretch / Montage_KO_Green_LinearStretch\n');
fprintf(fid, 'Montage of knockout (KO) images after linear histogram stretching (1%%-99%% saturation). ');
fprintf(fid, 'Each panel shows an individual animal labeled by number. ');
fprintf(fid, 'Linear stretching was applied using group-specific intensity limits to normalize within-group variation while preserving relative intensity differences. ');
fprintf(fid, 'N = %d KO animals.\n\n', num_KO);

fprintf(fid, '--- STEP 2: HISTOGRAM MONTAGES ---\n\n');

fprintf(fid, 'Figure: Montage_WT_Red / Montage_WT_Green / Montage_KO_Red / Montage_KO_Green\n');
fprintf(fid, 'One-dimensional intensity profiles computed by binning images along the X-axis (bin size = %d pixels, %.2f µm). ', ...
        config.binSize_analysis, config.binSize_analysis * config.pixelSize);
fprintf(fid, 'Each colored trace represents an individual animal. ');
fprintf(fid, 'Black solid line indicates the group mean; black dashed lines indicate ± 1 standard deviation. ');
fprintf(fid, 'Profiles were computed prior to alignment.\n\n');

fprintf(fid, '--- STEP 3: ALIGNMENT AND COMPARISON ---\n\n');

fprintf(fid, 'Figure: BeforeAfter_[Channel]_[Genotype]\n');
fprintf(fid, 'Comparison of intensity profiles before and after valley-based alignment. ');
fprintf(fid, 'Left panel: Original profiles in pixel coordinates. ');
fprintf(fid, 'Right panel: Aligned profiles in micron coordinates, centered at the midline (x = 0). ');
fprintf(fid, 'Gray traces show individual animals; black lines show mean ± SD. ');
fprintf(fid, 'Alignment was performed using bin %d data for robust valley detection, then applied to bin %d data for visualization.\n\n', ...
        config.binSize_alignment, config.binSize_analysis);

fprintf(fid, 'Figure: WTvsKO_[Channel]_original_bin20_[colorscheme] / WTvsKO_[Channel]_original_bin80_[colorscheme]\n');
fprintf(fid, 'Comparison of aligned intensity profiles between WT and KO groups. ');
fprintf(fid, 'Available color schemes: blue_red (WT=blue, KO=red), red_blue (WT=red, KO=blue). ');
fprintf(fid, 'Solid line: mean; shaded region: ± SD. ');
fprintf(fid, 'Vertical dashed line indicates the midline (x = 0). ');
fprintf(fid, 'Title shows Mann-Whitney U test p-values and post-hoc power for negative (x < 0) and positive (x > 0) halves. ');
fprintf(fid, 'Asterisks (*) indicate significant differences where both p < %.2f and power > %.1f. ', config.alpha, config.powerThreshold);
fprintf(fid, 'WT N = %d, KO N = %d.\n\n', num_WT, num_KO);

fprintf(fid, 'Figure: WTvsKO_[Channel]_normalized_bin20_[colorscheme] / WTvsKO_[Channel]_normalized_bin80_[colorscheme]\n');
fprintf(fid, 'Same as above, but with profiles normalized to 0-1 range within each animal to emphasize shape differences independent of absolute intensity.\n\n');

fprintf(fid, 'Figure: Stats_[Channel]_original / Stats_[Channel]_normalized\n');
fprintf(fid, 'Statistical summary for WT vs KO comparison. ');
fprintf(fid, 'Top left: Group means ± SD for negative and positive halves. ');
fprintf(fid, 'Top right: Mann-Whitney U test p-values with significance threshold (α = %.2f). ', config.alpha);
fprintf(fid, 'Bottom left: Cohen''s d effect sizes with conventional thresholds (small = 0.2, medium = 0.5, large = 0.8). ');
fprintf(fid, 'Bottom right: Post-hoc statistical power with threshold (%.1f). ', config.powerThreshold);
fprintf(fid, 'Asterisks indicate results meeting both significance and power criteria.\n\n');

fprintf(fid, '--- QUALITY CONTROL PLOTS ---\n\n');

fprintf(fid, 'Figure: ValleyDetect_[Channel]_[Genotype]\n');
fprintf(fid, 'Valley detection diagnostic showing intensity profiles with detected valley positions (circles). ');
fprintf(fid, 'Blue dashed line indicates the mean valley position used for alignment. ');
fprintf(fid, 'Green dotted lines indicate the search range for valley detection (%d-%d pixels). ', ...
        config.valleySearchMin, config.valleySearchMax);
fprintf(fid, 'Valley detection was performed on bin %d data for robustness.\n\n', config.binSize_alignment);

fprintf(fid, 'Figure: AlignmentQC_[Channel]_[Genotype]\n');
fprintf(fid, 'Alignment quality control showing shift magnitudes and profile deviations. ');
fprintf(fid, 'Left: Pixel shifts applied to each animal for alignment. ');
fprintf(fid, 'Right: Root-mean-square deviation of each profile from the group mean. ');
fprintf(fid, 'Red bars indicate potential outliers exceeding %.1f standard deviations from the mean deviation.\n\n', config.outlierThreshold);

fprintf(fid, 'Figure: Top5_[Channel]_[Genotype]\n');
fprintf(fid, 'Top 5 profiles with largest deviations from the group mean, shown for quality control. ');
fprintf(fid, 'Black line: group mean. Colored lines: individual outlier profiles.\n\n');

fprintf(fid, '--- HEATMAPS ---\n\n');

fprintf(fid, 'All heatmaps display average images using a GLOBAL intensity scale computed across all four images ');
fprintf(fid, '(Red WT, Red KO, Green WT, Green KO) to enable direct visual comparison between channels and genotypes. ');
fprintf(fid, 'Heatmaps are saved as PNG, FIG, EPS, and PDF formats.\n\n');

fprintf(fid, 'Figure: Heatmap_[Channel]_WT_greater_plasma / Heatmap_[Channel]_KO_greater_plasma\n');
fprintf(fid, 'Spatial visualization of regions where WT or KO mean intensity exceeds the other. ');
fprintf(fid, 'Left: WT mean intensity image (grayscale, global scale). ');
fprintf(fid, 'Center: KO mean intensity image (grayscale, global scale). ');
fprintf(fid, 'Right: Difference map using plasma colormap. ');
fprintf(fid, 'Difference color scale fixed at 0-150 for consistency across all channels. ');
fprintf(fid, 'X-axis centered at midline (0), ticks at -300 to +300 µm. ');
fprintf(fid, 'Y-axis runs from 0 at bottom to maximum at top. ');
fprintf(fid, 'White dashed vertical line indicates the aligned midline (x = 0). ');
fprintf(fid, 'Resolution: 900 DPI.\n\n');

fprintf(fid, 'Figure: Heatmap_[Channel]_WT_greater_white_red / Heatmap_[Channel]_KO_greater_white_blue\n');
fprintf(fid, 'Difference maps using white-to-color scheme (Scheme 1). ');
fprintf(fid, 'WT > KO shown with white-to-red colormap; KO > WT shown with white-to-blue colormap. ');
fprintf(fid, 'White indicates zero difference, color intensity indicates magnitude of difference. ');
fprintf(fid, 'Difference color scale fixed at 0-150 for consistency. ');
fprintf(fid, 'Black dashed vertical line indicates midline. Resolution: 900 DPI.\n\n');

fprintf(fid, 'Figure: Heatmap_[Channel]_WT_greater_white_blue / Heatmap_[Channel]_KO_greater_white_red\n');
fprintf(fid, 'Difference maps using reversed white-to-color scheme (Scheme 2). ');
fprintf(fid, 'WT > KO shown with white-to-blue colormap; KO > WT shown with white-to-red colormap. ');
fprintf(fid, 'White indicates zero difference, color intensity indicates magnitude of difference. ');
fprintf(fid, 'Difference color scale fixed at 0-150 for consistency. ');
fprintf(fid, 'Black dashed vertical line indicates midline. Resolution: 900 DPI.\n\n');

fprintf(fid, 'Figure: Heatmap_[Channel]_combined_blue_white_red\n');
fprintf(fid, 'Combined difference map using blue-white-red diverging colormap. ');
fprintf(fid, 'Blue indicates KO > WT (negative values), white indicates equal intensity (0), red indicates WT > KO (positive values). ');
fprintf(fid, 'Difference color scale fixed at -150 to +150 for consistency and symmetry around zero. ');
fprintf(fid, 'Black dashed vertical line indicates midline. Resolution: 900 DPI.\n\n');

fprintf(fid, 'Figure: Heatmap_[Channel]_combined_red_white_blue\n');
fprintf(fid, 'Combined difference map using red-white-blue diverging colormap (reverse). ');
fprintf(fid, 'Red indicates KO > WT (negative values), white indicates equal intensity (0), blue indicates WT > KO (positive values). ');
fprintf(fid, 'Difference color scale fixed at -150 to +150 for consistency and symmetry around zero. ');
fprintf(fid, 'Black dashed vertical line indicates midline. Resolution: 900 DPI.\n\n');

fprintf(fid, '--- AVERAGE IMAGES ---\n\n');

fprintf(fid, 'All average images use a GLOBAL intensity scale computed across all four images ');
fprintf(fid, '(Red WT, Red KO, Green WT, Green KO) to enable direct visual comparison.\n\n');

fprintf(fid, 'File: Average_[Genotype]_[Channel].tif\n');
fprintf(fid, 'Pixel-wise mean intensity image computed from all aligned (shifted) images within each group. ');
fprintf(fid, 'Saved as 16-bit TIFF with original intensity values. ');
fprintf(fid, 'WT: N = %d, KO: N = %d.\n\n', num_WT, num_KO);

fprintf(fid, 'File: Average_[Genotype]_[Channel]_8bit.tif\n');
fprintf(fid, '8-bit version scaled individually (each image uses its own min-max range).\n\n');

fprintf(fid, 'File: Average_[Genotype]_[Channel]_8bit_globalscale.tif\n');
fprintf(fid, '8-bit version scaled using GLOBAL min-max range across all four average images. ');
fprintf(fid, 'These images are directly comparable and match the heatmap display scaling.\n\n');

fprintf(fid, '--- DIFFERENCE MAP IMAGES (for manual figure assembly) ---\n\n');

fprintf(fid, 'All difference map TIF files are saved as uncompressed RGB images with colormaps applied, ');
fprintf(fid, 'suitable for importing into vector graphics software (e.g., Illustrator, Inkscape). ');
fprintf(fid, 'Images are correctly oriented (tissue top at image top).\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_WT_greater_plasma.tif\n');
fprintf(fid, 'RGB TIF showing regions where WT > KO using plasma colormap. Scale: 0-150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_KO_greater_plasma.tif\n');
fprintf(fid, 'RGB TIF showing regions where KO > WT using plasma colormap. Scale: 0-150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_WT_greater_white_red.tif\n');
fprintf(fid, 'RGB TIF showing regions where WT > KO using white-to-red colormap. Scale: 0-150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_KO_greater_white_blue.tif\n');
fprintf(fid, 'RGB TIF showing regions where KO > WT using white-to-blue colormap. Scale: 0-150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_WT_greater_white_blue.tif\n');
fprintf(fid, 'RGB TIF showing regions where WT > KO using white-to-blue colormap (reverse scheme). Scale: 0-150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_KO_greater_white_red.tif\n');
fprintf(fid, 'RGB TIF showing regions where KO > WT using white-to-red colormap (reverse scheme). Scale: 0-150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_combined_blue_white_red.tif\n');
fprintf(fid, 'RGB TIF showing full WT-KO difference using blue-white-red diverging colormap. ');
fprintf(fid, 'Blue = KO > WT, White = equal, Red = WT > KO. Scale: -150 to +150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_combined_red_white_blue.tif\n');
fprintf(fid, 'RGB TIF showing full WT-KO difference using red-white-blue diverging colormap. ');
fprintf(fid, 'Red = KO > WT, White = equal, Blue = WT > KO. Scale: -150 to +150.\n\n');

fprintf(fid, 'File: DifferenceMap_[Channel]_data.mat\n');
fprintf(fid, 'MATLAB data file containing raw difference values and scale parameters: ');
fprintf(fid, 'diff_raw (WT-KO), diff_WT_greater, diff_KO_greater, x_microns, y_microns, diff_scale_max, combined_scale.\n\n');

fclose(fid);
fprintf('    Saved: Figure_Legends.txt\n');

% ===== 2. IMAGE ANALYSIS METHODS =====
imageMethodsFile = fullfile(docDir, 'Methods_ImageAnalysis.txt');
fid = fopen(char(imageMethodsFile), 'w');

fprintf(fid, 'IMAGE ANALYSIS METHODS\n');
fprintf(fid, '======================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'For use in Methods section:\n\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'COMPREHENSIVE DESCRIPTION\n\n');

fprintf(fid, '1. IMAGE ACQUISITION AND PREPROCESSING\n');
fprintf(fid, '--------------------------------------\n\n');

fprintf(fid, 'Fluorescence images were processed using a custom MATLAB pipeline (R2020a or later). ');
fprintf(fid, 'Two-channel images containing red and green fluorescence signals were processed. ');
fprintf(fid, 'Each channel was extracted and analyzed separately to enable channel-specific quantification of signal intensity and spatial distribution.\n\n');

fprintf(fid, '2. INTENSITY NORMALIZATION (LINEAR HISTOGRAM STRETCHING)\n');
fprintf(fid, '--------------------------------------------------------\n\n');

fprintf(fid, 'To account for variability in staining intensity and imaging conditions across samples, ');
fprintf(fid, 'images were normalized using linear histogram stretching. ');
fprintf(fid, 'This approach rescales pixel intensity values to span the full dynamic range of the image format.\n\n');

fprintf(fid, 'Procedure:\n');
fprintf(fid, '  a) For each genotype group (WT and KO) and each channel (Red and Green) separately:\n');
fprintf(fid, '     - All images within the group were pooled\n');
fprintf(fid, '     - The 1st percentile (P1) and 99th percentile (P99) intensity values were computed across all pooled pixels\n');
fprintf(fid, '     - These percentile limits were used rather than absolute min/max to exclude extreme outliers (hot pixels, dead pixels, artifacts)\n\n');

fprintf(fid, '  b) Linear transformation applied to each image:\n');
fprintf(fid, '     - Output = (Input - P1) / (P99 - P1)\n');
fprintf(fid, '     - Values below P1 were clipped to 0\n');
fprintf(fid, '     - Values above P99 were clipped to 1 (or maximum bit depth)\n\n');

fprintf(fid, 'Rationale:\n');
fprintf(fid, '  - Group-specific normalization preserves relative intensity differences between genotypes\n');
fprintf(fid, '  - Within-group normalization reduces variability from imaging conditions (exposure, illumination)\n');
fprintf(fid, '  - Linear stretching (vs. nonlinear methods like histogram equalization) preserves the relative relationships between intensity values\n');
fprintf(fid, '  - 1%%-99%% saturation limits provide robustness to outlier pixels\n\n');

fprintf(fid, '3. ONE-DIMENSIONAL INTENSITY PROFILE COMPUTATION\n');
fprintf(fid, '------------------------------------------------\n\n');

fprintf(fid, 'To quantify the spatial distribution of fluorescence intensity across the tissue, ');
fprintf(fid, 'two-dimensional images were reduced to one-dimensional intensity profiles along the medial-lateral axis.\n\n');

fprintf(fid, 'Procedure:\n');
fprintf(fid, '  a) Images were divided into vertical bins along the X-axis (medial-lateral direction)\n');
fprintf(fid, '  b) Within each bin, all pixel intensities were averaged along the Y-axis (dorsal-ventral direction)\n');
fprintf(fid, '  c) This produced a single intensity value per bin, yielding a 1D profile\n\n');

fprintf(fid, 'Two bin sizes were used:\n');
fprintf(fid, '  - Alignment bin size: %d pixels (%.2f µm) - larger bins for robust valley detection\n', ...
        config.binSize_alignment, config.binSize_alignment * config.pixelSize);
fprintf(fid, '  - Analysis/plotting bin size: %d pixels (%.2f µm) - smaller bins for higher spatial resolution in final plots\n\n', ...
        config.binSize_analysis, config.binSize_analysis * config.pixelSize);

fprintf(fid, 'The dual-bin strategy ensures:\n');
fprintf(fid, '  - Robust alignment using smoothed profiles (less noise sensitivity)\n');
fprintf(fid, '  - High-resolution visualization of spatial patterns in final outputs\n\n');

fprintf(fid, '4. MIDLINE ALIGNMENT (VALLEY DETECTION)\n');
fprintf(fid, '---------------------------------------\n\n');

fprintf(fid, 'To enable valid comparisons across animals, all images were spatially aligned to a common anatomical reference point: ');
fprintf(fid, 'the midline, identified as the local minimum (valley) in the intensity profile.\n\n');

fprintf(fid, 'Procedure:\n');
fprintf(fid, '  a) Valley detection was performed on smoothed profiles (bin %d data)\n', config.binSize_alignment);
fprintf(fid, '  b) The intensity minimum was identified within a constrained search window (%d-%d pixels from left edge)\n', ...
        config.valleySearchMin, config.valleySearchMax);
fprintf(fid, '     - This constraint prevents false detection of minima at image edges or in unexpected regions\n');
fprintf(fid, '  c) The valley position was recorded for each animal (in pixel coordinates)\n');
fprintf(fid, '  d) The mean valley position was computed across all animals within each genotype group\n');
fprintf(fid, '     - WT mean valley position: %.1f pixels\n', mean_valley_WT);
fprintf(fid, '     - KO mean valley position: %.1f pixels\n', mean_valley_KO);
fprintf(fid, '  e) Individual images were shifted horizontally to align their valleys to the group mean position\n');
fprintf(fid, '     - Shift magnitude = (individual valley position) - (mean valley position)\n');
fprintf(fid, '     - Positive shifts move the image left; negative shifts move it right\n\n');

fprintf(fid, 'Post-alignment coordinate system:\n');
fprintf(fid, '  - X = 0 corresponds to the aligned midline (valley)\n');
fprintf(fid, '  - X < 0: positions to the left of the midline (negative half)\n');
fprintf(fid, '  - X > 0: positions to the right of the midline (positive half)\n');
fprintf(fid, '  - Left/right are IMAGE coordinates. Ipsilateral/contralateral are defined relative\n');
fprintf(fid, '    to the injected eye, and because each channel reports a different eye the mapping\n');
fprintf(fid, '    from image half to laterality is channel-specific:\n');
fprintf(fid, '      Red channel:   ipsilateral = x < 0,  contralateral = x > 0\n');
fprintf(fid, '      Green channel: ipsilateral = x > 0,  contralateral = x < 0\n');
fprintf(fid, '  - Coordinates converted from pixels to microns using calibration factor: %.4f µm/pixel\n\n', config.pixelSize);

fprintf(fid, '5. PROFILE NORMALIZATION (0-1 SCALING)\n');
fprintf(fid, '--------------------------------------\n\n');

fprintf(fid, 'In addition to the original intensity profiles, a secondary normalization was applied to enable ');
fprintf(fid, 'comparison of spatial distribution patterns independent of absolute intensity levels.\n\n');

fprintf(fid, 'Procedure (per animal, per channel):\n');
fprintf(fid, '  a) For each aligned intensity profile, identify:\n');
fprintf(fid, '     - Profile minimum: min_val = min(profile)\n');
fprintf(fid, '     - Profile maximum: max_val = max(profile)\n');
fprintf(fid, '  b) Apply min-max normalization:\n');
fprintf(fid, '     - Normalized = (Original - min_val) / (max_val - min_val)\n');
fprintf(fid, '  c) Result: all profiles span the range [0, 1]\n\n');

fprintf(fid, 'Interpretation:\n');
fprintf(fid, '  - 0-1 normalized profiles reveal SHAPE differences (spatial distribution patterns)\n');
fprintf(fid, '  - Original profiles reveal MAGNITUDE differences (absolute intensity levels)\n');
fprintf(fid, '  - A genotype difference in original but not normalized data suggests intensity change without distribution change\n');
fprintf(fid, '  - A genotype difference in normalized but not original data suggests distribution change without overall intensity change\n');
fprintf(fid, '  - Differences in both suggest changes in both intensity and distribution\n\n');

fprintf(fid, '6. AVERAGE IMAGE COMPUTATION\n');
fprintf(fid, '----------------------------\n\n');

fprintf(fid, 'To visualize the mean spatial pattern of fluorescence for each group, pixel-wise average images were computed.\n\n');

fprintf(fid, 'Procedure:\n');
fprintf(fid, '  a) All aligned (shifted) images within a genotype group were stacked\n');
fprintf(fid, '  b) For each pixel position (x, y), the mean intensity across all images was computed\n');
fprintf(fid, '  c) Result: a single average image representing the group mean spatial pattern\n\n');

fprintf(fid, 'Output files:\n');
fprintf(fid, '  - 16-bit TIFF: preserves full dynamic range for quantitative analysis\n');
fprintf(fid, '  - 8-bit TIFF: scaled for visualization purposes\n\n');

fprintf(fid, '7. DIFFERENCE MAP COMPUTATION\n');
fprintf(fid, '-----------------------------\n\n');

fprintf(fid, 'To visualize spatial regions where genotypes differ in mean intensity, difference maps were computed.\n\n');

fprintf(fid, 'Procedure:\n');
fprintf(fid, '  a) Difference = WT_average - KO_average (pixel-wise subtraction)\n');
fprintf(fid, '  b) Positive values indicate WT > KO (higher intensity in WT)\n');
fprintf(fid, '  c) Negative values indicate KO > WT (higher intensity in KO)\n');
fprintf(fid, '  d) Zero values indicate equal mean intensity between groups\n\n');

fprintf(fid, 'Visualization:\n');
fprintf(fid, '  - "Greater" maps: show only regions where one genotype exceeds the other (unidirectional)\n');
fprintf(fid, '  - "Combined" maps: show both directions using diverging colormaps centered at zero\n');
fprintf(fid, '  - Fixed color scales (0-150 for greater; ±150 for combined) enable comparison across channels\n\n');

fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'CONCISE METHODS PARAGRAPH (for publication):\n\n');

fprintf(fid, 'Fluorescence images were processed using a custom MATLAB pipeline. ');
fprintf(fid, 'Red and green channel images were extracted and normalized separately for each genotype group ');
fprintf(fid, 'using linear histogram stretching with 1%% and 99%% intensity saturation limits, ');
fprintf(fid, 'which normalizes within-group intensity variation while preserving relative differences between groups. ');
fprintf(fid, 'One-dimensional intensity profiles were computed by averaging pixel intensities along the Y-axis ');
fprintf(fid, 'within spatial bins of %d pixels (%.2f µm) for alignment and %d pixels (%.2f µm) for visualization. ', ...
        config.binSize_alignment, config.binSize_alignment * config.pixelSize, ...
        config.binSize_analysis, config.binSize_analysis * config.pixelSize);
fprintf(fid, 'Images were aligned to the midline using a valley-detection algorithm that identified the intensity minimum ');
fprintf(fid, 'within a search window of %d-%d pixels. ', config.valleySearchMin, config.valleySearchMax);
fprintf(fid, 'The mean valley position across all animals was computed, and individual images were shifted horizontally ');
fprintf(fid, 'to align their valleys to this common reference point. ');
fprintf(fid, 'To assess spatial distribution patterns independent of absolute intensity, profiles were additionally normalized to the 0-1 range within each animal. ');
fprintf(fid, 'Average intensity images were generated by computing the pixel-wise mean across all aligned images within each group, ');
fprintf(fid, 'and difference maps were calculated by subtracting the KO group average from the WT group average. ');
fprintf(fid, 'Spatial coordinates were converted from pixels to microns using a calibration factor of %.4f µm/pixel, ', config.pixelSize);
fprintf(fid, 'with the X-axis centered at the aligned midline (x = 0).\n\n');

fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'PARAMETERS USED:\n\n');
fprintf(fid, '- Pixel size: %.4f µm/pixel\n', config.pixelSize);
fprintf(fid, '- Alignment bin size: %d pixels (%.2f µm)\n', config.binSize_alignment, config.binSize_alignment * config.pixelSize);
fprintf(fid, '- Analysis bin size: %d pixels (%.2f µm)\n', config.binSize_analysis, config.binSize_analysis * config.pixelSize);
fprintf(fid, '- Valley search range: %d-%d pixels\n', config.valleySearchMin, config.valleySearchMax);
fprintf(fid, '- Linear stretch saturation: 1%%-99%%\n');
fprintf(fid, '- Outlier threshold: %.1f standard deviations\n', config.outlierThreshold);
fprintf(fid, '- Heatmap difference scale: 0-150 (greater maps), ±150 (combined maps)\n');
fprintf(fid, '- Heatmap resolution: 900 DPI\n');
fprintf(fid, '- Sample sizes: WT N = %d, KO N = %d\n', num_WT, num_KO);

fclose(fid);
fprintf('    Saved: Methods_ImageAnalysis.txt\n');

% ===== 3. STATISTICAL ANALYSIS METHODS =====
statsMethodsFile = fullfile(docDir, 'Methods_StatisticalAnalysis.txt');
fid = fopen(char(statsMethodsFile), 'w');

fprintf(fid, 'STATISTICAL ANALYSIS METHODS\n');
fprintf(fid, '============================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'For use in Methods section:\n\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'COMPREHENSIVE DESCRIPTION\n\n');

fprintf(fid, '1. RATIONALE FOR STATISTICAL APPROACH\n');
fprintf(fid, '-------------------------------------\n\n');

fprintf(fid, 'The statistical design addresses several key considerations:\n\n');

fprintf(fid, '  a) Biological replication: Each animal represents an independent biological replicate.\n');
fprintf(fid, '     - Multiple measurements within an animal (e.g., multiple bins in a profile) are technical replicates\n');
fprintf(fid, '     - Statistical tests must be performed at the level of biological replicates (animals) to avoid pseudoreplication\n\n');

fprintf(fid, '  b) Non-parametric approach: The Mann-Whitney U test was chosen because:\n');
fprintf(fid, '     - Sample sizes are relatively small (WT N = %d, KO N = %d)\n', num_WT, num_KO);
fprintf(fid, '     - Normality of the underlying distribution cannot be assumed\n');
fprintf(fid, '     - The test is robust to outliers and skewed distributions\n\n');

fprintf(fid, '  c) Spatial subdivision: Profiles were divided at the aligned midline into two halves:\n');
fprintf(fid, '     - Negative half (x < 0): left of the midline\n');
fprintf(fid, '     - Positive half (x > 0): right of the midline\n');
fprintf(fid, '     - Laterality is channel-specific, since each channel reports a different eye:\n');
fprintf(fid, '         Red channel:   ipsilateral = x < 0,  contralateral = x > 0\n');
fprintf(fid, '         Green channel: ipsilateral = x > 0,  contralateral = x < 0\n');
fprintf(fid, '     - This enables detection of region-specific differences\n\n');

fprintf(fid, '  d) Dual analysis (original and normalized):\n');
fprintf(fid, '     - Original intensity data: detects differences in absolute fluorescence levels\n');
fprintf(fid, '     - 0-1 normalized data: detects differences in spatial distribution patterns\n');
fprintf(fid, '     - Together, these distinguish between intensity vs. distribution phenotypes\n\n');

fprintf(fid, '2. SUMMARY STATISTIC COMPUTATION\n');
fprintf(fid, '--------------------------------\n\n');

fprintf(fid, 'For each animal, a single summary statistic was computed for each region to ensure statistical independence.\n\n');

fprintf(fid, 'Procedure:\n');
fprintf(fid, '  a) For each aligned profile, identify bins with x < 0 (negative half) and x > 0 (positive half)\n');
fprintf(fid, '  b) Compute the arithmetic mean of intensity values across all bins in each half\n');
fprintf(fid, '  c) Result: one value per animal per region (2 values total per animal)\n\n');

fprintf(fid, 'This approach:\n');
fprintf(fid, '  - Avoids pseudoreplication by collapsing multiple bins into one summary per animal\n');
fprintf(fid, '  - Provides a robust central tendency measure\n');
fprintf(fid, '  - Enables direct comparison between genotypes at the animal level\n\n');

fprintf(fid, '3. MANN-WHITNEY U TEST\n');
fprintf(fid, '----------------------\n\n');

fprintf(fid, 'The Mann-Whitney U test (also called Wilcoxon rank-sum test) was used to compare WT and KO groups.\n\n');

fprintf(fid, 'Test characteristics:\n');
fprintf(fid, '  - Non-parametric: does not assume normal distribution\n');
fprintf(fid, '  - Compares the ranks of values between groups rather than the values themselves\n');
fprintf(fid, '  - Null hypothesis: the distribution of values is the same in both groups\n');
fprintf(fid, '  - Alternative hypothesis: one group tends to have larger values than the other\n\n');

fprintf(fid, 'Implementation:\n');
fprintf(fid, '  - MATLAB function: ranksum()\n');
fprintf(fid, '  - Two-tailed test (either group could have higher values)\n');
fprintf(fid, '  - Significance threshold (α): %.2f\n\n', config.alpha);

fprintf(fid, 'Test outputs:\n');
fprintf(fid, '  - U statistic: sum of ranks in one group minus the minimum possible sum\n');
fprintf(fid, '  - p-value: probability of observing a result as extreme as the data under the null hypothesis\n\n');

fprintf(fid, 'Interpretation:\n');
fprintf(fid, '  - p < %.2f: reject null hypothesis; groups differ significantly\n', config.alpha);
fprintf(fid, '  - p ≥ %.2f: fail to reject null hypothesis; no significant difference detected\n\n', config.alpha);

fprintf(fid, '4. EFFECT SIZE (COHEN''S D)\n');
fprintf(fid, '--------------------------\n\n');

fprintf(fid, 'Effect size quantifies the magnitude of the difference between groups, independent of sample size.\n\n');

fprintf(fid, 'Cohen''s d formula:\n');
fprintf(fid, '  d = |mean(WT) - mean(KO)| / pooled_SD\n\n');

fprintf(fid, 'Where pooled standard deviation:\n');
fprintf(fid, '  pooled_SD = sqrt(((n_WT - 1) * SD_WT^2 + (n_KO - 1) * SD_KO^2) / (n_WT + n_KO - 2))\n\n');

fprintf(fid, 'Conventional interpretation thresholds (Cohen, 1988):\n');
fprintf(fid, '  - d = 0.2: small effect\n');
fprintf(fid, '  - d = 0.5: medium effect\n');
fprintf(fid, '  - d = 0.8: large effect\n\n');

fprintf(fid, 'Note: The absolute value is used for d, making it always positive. ');
fprintf(fid, 'The direction of the effect (which group is higher) is determined separately.\n\n');

fprintf(fid, '5. POST-HOC POWER ANALYSIS\n');
fprintf(fid, '--------------------------\n\n');

fprintf(fid, 'Statistical power was computed post-hoc to assess the reliability of the findings.\n\n');

fprintf(fid, 'Power definition:\n');
fprintf(fid, '  - Power = 1 - β, where β is the probability of Type II error (false negative)\n');
fprintf(fid, '  - Power = probability of detecting a true effect if one exists\n\n');

fprintf(fid, 'Computation method:\n');
fprintf(fid, '  - Based on non-central t-distribution\n');
fprintf(fid, '  - Uses observed effect size (Cohen''s d) and sample sizes\n');
fprintf(fid, '  - Non-centrality parameter: δ = d * sqrt(n_WT * n_KO / (n_WT + n_KO))\n');
fprintf(fid, '  - Power computed as area under non-central t-distribution beyond critical t-value\n\n');

fprintf(fid, 'Interpretation:\n');
fprintf(fid, '  - Power > %.1f: adequate power; reliable result (either significant or non-significant)\n', config.powerThreshold);
fprintf(fid, '  - Power < %.1f: low power; non-significant results should be interpreted with caution\n', config.powerThreshold);
fprintf(fid, '  - Low power with significant result: still valid, but effect may be overestimated\n\n');

fprintf(fid, 'Note: Post-hoc power analysis is most useful for interpreting non-significant results. ');
fprintf(fid, 'A non-significant result with low power may indicate insufficient sample size rather than absence of effect.\n\n');

fprintf(fid, '6. COMBINED SIGNIFICANCE AND POWER CRITERION\n');
fprintf(fid, '--------------------------------------------\n\n');

fprintf(fid, 'Results are flagged as "robust" (marked with *) only when BOTH criteria are met:\n');
fprintf(fid, '  - p-value < %.2f (statistically significant)\n', config.alpha);
fprintf(fid, '  - power > %.1f (adequately powered)\n\n', config.powerThreshold);

fprintf(fid, 'This dual criterion:\n');
fprintf(fid, '  - Reduces false positives from underpowered significant results\n');
fprintf(fid, '  - Highlights findings with both statistical and practical significance\n');
fprintf(fid, '  - Provides more reliable basis for biological interpretation\n\n');

fprintf(fid, '7. INTERPRETATION OF ORIGINAL VS. NORMALIZED COMPARISONS\n');
fprintf(fid, '--------------------------------------------------------\n\n');

fprintf(fid, 'Two parallel analyses were performed: on original intensity data and on 0-1 normalized data.\n\n');

fprintf(fid, 'Original intensity analysis:\n');
fprintf(fid, '  - Tests for differences in absolute fluorescence intensity\n');
fprintf(fid, '  - Sensitive to changes in: protein expression levels, staining efficiency, imaging parameters\n');
fprintf(fid, '  - Significant difference = genotypes differ in how much signal is present\n\n');

fprintf(fid, 'Normalized (0-1) intensity analysis:\n');
fprintf(fid, '  - Tests for differences in spatial distribution pattern (shape)\n');
fprintf(fid, '  - Each animal''s profile is scaled to span [0, 1], removing absolute intensity information\n');
fprintf(fid, '  - Significant difference = genotypes differ in WHERE signal is distributed, independent of total amount\n\n');

fprintf(fid, 'Combined interpretation:\n\n');

fprintf(fid, '  | Original | Normalized | Interpretation                                      |\n');
fprintf(fid, '  |----------|------------|----------------------------------------------------|\n');
fprintf(fid, '  | Sig      | Sig        | Both intensity and distribution differ              |\n');
fprintf(fid, '  | Sig      | Non-sig    | Intensity differs, but distribution pattern similar |\n');
fprintf(fid, '  | Non-sig  | Sig        | Distribution differs, but total intensity similar   |\n');
fprintf(fid, '  | Non-sig  | Non-sig    | No detectable difference in either metric           |\n\n');

fprintf(fid, 'This dual approach provides more nuanced biological insight than either analysis alone.\n\n');

fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'CONCISE METHODS PARAGRAPH (for publication):\n\n');

fprintf(fid, 'To compare intensity distributions between WT (N = %d) and KO (N = %d) groups, ', num_WT, num_KO);
fprintf(fid, 'aligned intensity profiles were divided at the midline into a negative (x < 0) and a positive (x > 0) half. ');
fprintf(fid, 'Because the two eyes were injected with spectrally distinct tracers, the half corresponding to the SCN ipsilateral ');
fprintf(fid, 'to the injected eye differs between channels: the ipsilateral half is x < 0 for the red channel and x > 0 for the green channel. ');
fprintf(fid, 'For each animal, a single summary statistic was computed as the mean intensity across all bins within each half, ');
fprintf(fid, 'ensuring that each biological replicate contributed exactly one value to the statistical comparison. ');
fprintf(fid, 'Group differences were assessed using the Mann-Whitney U test (Wilcoxon rank-sum test), ');
fprintf(fid, 'a non-parametric test appropriate for comparing independent samples without assuming normality. ');
fprintf(fid, 'Effect sizes were quantified using Cohen''s d, calculated as the absolute difference between group means divided by the pooled standard deviation, ');
fprintf(fid, 'with conventional thresholds for interpretation: small (d = 0.2), medium (d = 0.5), and large (d = 0.8). ');
fprintf(fid, 'Post-hoc statistical power was estimated using the non-central t-distribution based on the observed effect size and sample sizes. ');
fprintf(fid, 'Results were considered statistically robust when both p < %.2f and power > %.1f. ', config.alpha, config.powerThreshold);
fprintf(fid, 'To distinguish between changes in absolute intensity versus spatial distribution, analyses were performed on both original intensity data and ');
fprintf(fid, '0-1 normalized data (each profile independently scaled to the range [0, 1]). ');
fprintf(fid, 'All statistical analyses were performed in MATLAB.\n\n');

fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'PARAMETERS USED:\n\n');
fprintf(fid, '- Statistical test: Mann-Whitney U (Wilcoxon rank-sum test)\n');
fprintf(fid, '- Significance threshold (α): %.2f\n', config.alpha);
fprintf(fid, '- Power threshold: %.1f\n', config.powerThreshold);
fprintf(fid, '- Effect size metric: Cohen''s d (absolute value)\n');
fprintf(fid, '- Effect size thresholds: small (0.2), medium (0.5), large (0.8)\n');
fprintf(fid, '- Sample sizes: WT N = %d, KO N = %d\n', num_WT, num_KO);
fprintf(fid, '- Spatial regions: Negative half (x < 0), Positive half (x > 0)\n');
fprintf(fid, '- Data types analyzed: Original intensities, 0-1 normalized intensities\n');
fprintf(fid, '- Bin sizes analyzed: %d pixels (%.2f µm), %d pixels (%.2f µm)\n\n', ...
        config.binSize_analysis, config.binSize_analysis * config.pixelSize, ...
        config.binSize_alignment, config.binSize_alignment * config.pixelSize);

fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'STATISTICAL RESULTS SUMMARY:\n\n');

channels = {'Red', 'Green'};
dataTypes = {'original', 'normalized'};

for ch = 1:length(channels)
    fprintf(fid, '=== %s Channel ===\n\n', upper(channels{ch}));
    
    for dt = 1:length(dataTypes)
        stats = statsData.(normMethod).(channels{ch}).(dataTypes{dt});
        
        fprintf(fid, '%s Data:\n', dataTypes{dt});
        fprintf(fid, '  Negative half (x < 0) [%s]:\n', lateralityFor(channels{ch}, 'negative', config));
        fprintf(fid, '    WT mean ± SD: %.4f ± %.4f\n', stats.WT_negative_mean, stats.WT_negative_std);
        fprintf(fid, '    KO mean ± SD: %.4f ± %.4f\n', stats.KO_negative_mean, stats.KO_negative_std);
        fprintf(fid, '    Mann-Whitney U: %.1f\n', stats.U_negative);
        fprintf(fid, '    p-value: %.4f %s\n', stats.p_negative, getSigString(stats.p_negative, config.alpha));
        fprintf(fid, '    Cohen''s d: %.3f (%s effect)\n', stats.cohens_d_negative, getEffectSizeLabel(stats.cohens_d_negative));
        fprintf(fid, '    Power: %.3f %s\n', stats.power_negative, getPowerString(stats.power_negative, config.powerThreshold));
        fprintf(fid, '    Robust result: %s\n\n', getRobustString(stats.p_negative, stats.power_negative, config));
        
        fprintf(fid, '  Positive half (x > 0) [%s]:\n', lateralityFor(channels{ch}, 'positive', config));
        fprintf(fid, '    WT mean ± SD: %.4f ± %.4f\n', stats.WT_positive_mean, stats.WT_positive_std);
        fprintf(fid, '    KO mean ± SD: %.4f ± %.4f\n', stats.KO_positive_mean, stats.KO_positive_std);
        fprintf(fid, '    Mann-Whitney U: %.1f\n', stats.U_positive);
        fprintf(fid, '    p-value: %.4f %s\n', stats.p_positive, getSigString(stats.p_positive, config.alpha));
        fprintf(fid, '    Cohen''s d: %.3f (%s effect)\n', stats.cohens_d_positive, getEffectSizeLabel(stats.cohens_d_positive));
        fprintf(fid, '    Power: %.3f %s\n', stats.power_positive, getPowerString(stats.power_positive, config.powerThreshold));
        fprintf(fid, '    Robust result: %s\n\n', getRobustString(stats.p_positive, stats.power_positive, config));
    end
end

fclose(fid);
fprintf('    Saved: Methods_StatisticalAnalysis.txt\n');

fprintf('  Documentation files saved to: %s\n', docDir);
end

% Helper functions for documentation
function keepIdx = applyExclusions(idxList, channelName, imageInfo, config)
% Given a list of image indices (into imageInfo) for one genotype, return the
% subset to KEEP for the given channel after removing animals listed in
% config.exclude.(channelName). If no exclusion list exists for the channel,
% all indices are kept. Order is preserved.
    keepIdx = idxList;
    if ~isfield(config, 'exclude') || ~isfield(config.exclude, channelName)
        return;
    end
    excludeAnimals = config.exclude.(channelName);
    if isempty(excludeAnimals)
        return;
    end
    animalsForIdx = [imageInfo(idxList).animal];
    drop = ismember(animalsForIdx, excludeAnimals);
    keepIdx = idxList(~drop);
end

function label = lateralityFor(channelName, half, config)
% Return 'ipsilateral' or 'contralateral' for a given channel and image half.
% half is 'negative' (x < 0) or 'positive' (x > 0). The mapping is channel-specific
% because each channel reports a different eye (see config.ipsiHalf).
    if strcmp(config.ipsiHalf.(channelName), half)
        label = 'ipsilateral';
    else
        label = 'contralateral';
    end
end

function str = getSigString(p, alpha)
    if p < alpha
        str = '(significant)';
    else
        str = '(not significant)';
    end
end

function str = getEffectSizeLabel(d)
    if d < 0.2
        str = 'negligible';
    elseif d < 0.5
        str = 'small';
    elseif d < 0.8
        str = 'medium';
    else
        str = 'large';
    end
end

function str = getPowerString(power, threshold)
    if power >= threshold
        str = '(adequate)';
    else
        str = '(low)';
    end
end

function str = getRobustString(p, power, config)
    if p < config.alpha && power >= config.powerThreshold
        str = 'YES *';
    else
        str = 'NO';
    end
end

function cmap = bluewhitered(n)
% Create a blue-white-red colormap
if nargin < 1
    n = 64;
end

% Create colormap: blue -> white -> red
half_n = ceil(n/2);
blue_to_white = [linspace(0, 1, half_n)', linspace(0, 1, half_n)', ones(half_n, 1)];
white_to_red = [ones(n - half_n, 1), linspace(1, 0, n - half_n)', linspace(1, 0, n - half_n)'];

cmap = [blue_to_white; white_to_red];
end

function cmap = viridis(n)
% Create viridis colormap (perceptually uniform, colorblind-friendly)
if nargin < 1
    n = 256;
end

% Viridis colormap data points (from matplotlib)
viridis_data = [
    0.267004, 0.004874, 0.329415;
    0.282327, 0.140926, 0.457517;
    0.253935, 0.265254, 0.529983;
    0.206756, 0.371758, 0.553117;
    0.163625, 0.471133, 0.558148;
    0.127568, 0.566949, 0.550556;
    0.134692, 0.658636, 0.517649;
    0.266941, 0.748751, 0.440573;
    0.477504, 0.821444, 0.318195;
    0.741388, 0.873449, 0.149561;
    0.993248, 0.906157, 0.143936
];

% Interpolate to get n colors
x_orig = linspace(0, 1, size(viridis_data, 1));
x_new = linspace(0, 1, n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(x_orig, viridis_data(:, i), x_new, 'pchip');
end

% Ensure values are in [0, 1]
cmap = max(0, min(1, cmap));
end

function rgb_img = applyColormapToImage(data, cmap, min_val, max_val)
% Apply a colormap to grayscale data and return RGB image
% data: 2D matrix of values
% cmap: Nx3 colormap matrix
% min_val, max_val: range for scaling data to colormap indices

% Normalize data to 0-1 range based on min/max
normalized = (data - min_val) / (max_val - min_val);
normalized = max(0, min(1, normalized));  % Clamp to [0, 1]

% Convert to colormap indices (1 to size(cmap,1))
n_colors = size(cmap, 1);
indices = round(normalized * (n_colors - 1)) + 1;
indices = max(1, min(n_colors, indices));  % Ensure valid indices

% Create RGB image
[H, W] = size(data);
rgb_img = zeros(H, W, 3, 'uint8');

for c = 1:3
    channel = zeros(H, W);
    for i = 1:n_colors
        channel(indices == i) = cmap(i, c);
    end
    rgb_img(:,:,c) = uint8(channel * 255);
end
end

function cmap = plasma(n)
% Create plasma colormap (perceptually uniform, high contrast)
if nargin < 1
    n = 256;
end

% Plasma colormap data points (from matplotlib)
plasma_data = [
    0.050383, 0.029803, 0.527975;
    0.254627, 0.013882, 0.615419;
    0.417642, 0.000564, 0.658390;
    0.562738, 0.051545, 0.641509;
    0.692840, 0.165141, 0.564522;
    0.798216, 0.280197, 0.469538;
    0.881443, 0.392529, 0.383229;
    0.949217, 0.517763, 0.295662;
    0.988260, 0.652325, 0.211364;
    0.988648, 0.809579, 0.145357;
    0.940015, 0.975158, 0.131326
];

% Interpolate to get n colors
x_orig = linspace(0, 1, size(plasma_data, 1));
x_new = linspace(0, 1, n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(x_orig, plasma_data(:, i), x_new, 'pchip');
end

% Ensure values are in [0, 1]
cmap = max(0, min(1, cmap));
end

function cmap = blue_white_red(n)
% Create blue-white-red diverging colormap
% Blue = negative (KO > WT), White = 0, Red = positive (WT > KO)
if nargin < 1
    n = 256;
end

half_n = floor(n/2);

% Blue to white (for negative values)
blue_to_white = zeros(half_n, 3);
blue_to_white(:,1) = linspace(0, 1, half_n);      % R: 0 -> 1
blue_to_white(:,2) = linspace(0, 1, half_n);      % G: 0 -> 1
blue_to_white(:,3) = ones(half_n, 1);              % B: 1 -> 1

% White to red (for positive values)
white_to_red = zeros(n - half_n, 3);
white_to_red(:,1) = ones(n - half_n, 1);           % R: 1 -> 1
white_to_red(:,2) = linspace(1, 0, n - half_n);   % G: 1 -> 0
white_to_red(:,3) = linspace(1, 0, n - half_n);   % B: 1 -> 0

cmap = [blue_to_white; white_to_red];
end

function cmap = red_white_blue(n)
% Create red-white-blue diverging colormap (reverse of blue_white_red)
% Red = negative (KO > WT), White = 0, Blue = positive (WT > KO)
if nargin < 1
    n = 256;
end

half_n = floor(n/2);

% Red to white (for negative values)
red_to_white = zeros(half_n, 3);
red_to_white(:,1) = ones(half_n, 1);               % R: 1 -> 1
red_to_white(:,2) = linspace(0, 1, half_n);       % G: 0 -> 1
red_to_white(:,3) = linspace(0, 1, half_n);       % B: 0 -> 1

% White to blue (for positive values)
white_to_blue = zeros(n - half_n, 3);
white_to_blue(:,1) = linspace(1, 0, n - half_n);  % R: 1 -> 0
white_to_blue(:,2) = linspace(1, 0, n - half_n);  % G: 1 -> 0
white_to_blue(:,3) = ones(n - half_n, 1);          % B: 1 -> 1

cmap = [red_to_white; white_to_blue];
end

function cmap = white_to_red_cmap(n)
% Create white-to-red sequential colormap for WT > KO
if nargin < 1
    n = 256;
end

cmap = zeros(n, 3);
cmap(:,1) = ones(n, 1);                    % R: 1 -> 1
cmap(:,2) = linspace(1, 0, n)';           % G: 1 -> 0
cmap(:,3) = linspace(1, 0, n)';           % B: 1 -> 0
end

function cmap = white_to_blue_cmap(n)
% Create white-to-blue sequential colormap for KO > WT
if nargin < 1
    n = 256;
end

cmap = zeros(n, 3);
cmap(:,1) = linspace(1, 0, n)';           % R: 1 -> 0
cmap(:,2) = linspace(1, 0, n)';           % G: 1 -> 0
cmap(:,3) = ones(n, 1);                    % B: 1 -> 1
end

function cmap = magenta_white_green(n)
% Create magenta-white-green diverging colormap
% Magenta = negative (KO > WT), White = 0, Green = positive (WT > KO)
if nargin < 1
    n = 256;
end

half_n = floor(n/2);

% Magenta to white (for negative values)
magenta_to_white = zeros(half_n, 3);
magenta_to_white(:,1) = linspace(1, 1, half_n);   % R: 1 -> 1
magenta_to_white(:,2) = linspace(0, 1, half_n);   % G: 0 -> 1
magenta_to_white(:,3) = linspace(1, 1, half_n);   % B: 1 -> 1

% White to green (for positive values)
white_to_green = zeros(n - half_n, 3);
white_to_green(:,1) = linspace(1, 0, n - half_n); % R: 1 -> 0
white_to_green(:,2) = linspace(1, 0.8, n - half_n); % G: 1 -> 0.8 (keep some green)
white_to_green(:,3) = linspace(1, 0, n - half_n); % B: 1 -> 0

cmap = [magenta_to_white; white_to_green];
end

function cmap = green_white_magenta(n)
% Create green-white-magenta diverging colormap
% Green = negative (KO > WT), White = 0, Magenta = positive (WT > KO)
if nargin < 1
    n = 256;
end

half_n = floor(n/2);

% Green to white (for negative values)
green_to_white = zeros(half_n, 3);
green_to_white(:,1) = linspace(0, 1, half_n);     % R: 0 -> 1
green_to_white(:,2) = linspace(0.8, 1, half_n);   % G: 0.8 -> 1
green_to_white(:,3) = linspace(0, 1, half_n);     % B: 0 -> 1

% White to magenta (for positive values)
white_to_magenta = zeros(n - half_n, 3);
white_to_magenta(:,1) = linspace(1, 1, n - half_n);   % R: 1 -> 1
white_to_magenta(:,2) = linspace(1, 0, n - half_n);   % G: 1 -> 0
white_to_magenta(:,3) = linspace(1, 1, n - half_n);   % B: 1 -> 1

cmap = [green_to_white; white_to_magenta];
end

function cmap = plasma_diverging_mirror(n)
% Create diverging colormap using mirrored plasma
% Option A: reversed plasma (yellow→purple) for negative, white at center,
% regular plasma (purple→yellow) for positive
if nargin < 1
    n = 256;
end

% Get base plasma colormap
plasma_base = plasma(128);

half_n = floor(n/2);

% Reversed plasma to white (for negative values)
neg_half = zeros(half_n, 3);
plasma_reversed = flipud(plasma_base);
for i = 1:3
    % Interpolate reversed plasma and blend to white
    plasma_interp = interp1(linspace(0, 1, 128), plasma_reversed(:,i), linspace(0, 1, half_n));
    % Blend towards white at the end
    blend = linspace(0, 1, half_n)';
    neg_half(:,i) = plasma_interp' .* (1 - blend) + blend;
end

% White to plasma (for positive values)
pos_half = zeros(n - half_n, 3);
for i = 1:3
    plasma_interp = interp1(linspace(0, 1, 128), plasma_base(:,i), linspace(0, 1, n - half_n));
    % Blend from white
    blend = linspace(0, 1, n - half_n)';
    pos_half(:,i) = (1 - blend) + plasma_interp' .* blend;
end

cmap = [neg_half; pos_half];
cmap = max(0, min(1, cmap));
end

function cmap = plasma_diverging_complement(n)
% Create diverging colormap with cool tones for negative, plasma for positive
% Option B: Blue/cyan scale for negative (KO>WT), white at center,
% plasma for positive (WT>KO)
if nargin < 1
    n = 256;
end

% Get base plasma colormap
plasma_base = plasma(128);

half_n = floor(n/2);

% Cool blue/cyan to white (for negative values)
neg_half = zeros(half_n, 3);
% Start from deep blue/cyan, go to white
neg_half(:,1) = linspace(0.1, 1, half_n);    % R: dark -> white
neg_half(:,2) = linspace(0.3, 1, half_n);    % G: teal -> white
neg_half(:,3) = linspace(0.7, 1, half_n);    % B: blue -> white

% White to plasma (for positive values)
pos_half = zeros(n - half_n, 3);
for i = 1:3
    plasma_interp = interp1(linspace(0, 1, 128), plasma_base(:,i), linspace(0, 1, n - half_n));
    % Blend from white
    blend = linspace(0, 1, n - half_n)';
    pos_half(:,i) = (1 - blend) + plasma_interp' .* blend;
end

cmap = [neg_half; pos_half];
cmap = max(0, min(1, cmap));
end

%% ═════════════════════════════════════════════════════════════
%% PLOTTING FUNCTIONS
%% ═════════════════════════════════════════════════════════════

function createBeforeAfterPlots(alignedData, config)
channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};
normMethod = 'Linear';

for ch = 1:length(channels)
    for g = 1:length(genotypes)
        original = alignedData.(normMethod).(channels{ch}).(genotypes{g}).original_profiles;
        aligned = alignedData.(normMethod).(channels{ch}).(genotypes{g}).aligned_profiles;
        x_coords = alignedData.(normMethod).(channels{ch}).(genotypes{g}).x_coords;
        x_microns = alignedData.(normMethod).(channels{ch}).(genotypes{g}).x_microns;
        
        orig_mean = mean(original, 2);
        orig_std = std(original, 0, 2);
        aligned_mean = mean(aligned, 2);
        aligned_std = std(aligned, 0, 2);
        
        n_samples = size(aligned, 2);
        
        fig = figure('Position', [100, 100, 1400, 500], 'Visible', config.showFigures);
        set(fig, 'DefaultAxesFontName', 'Arial');
        set(fig, 'DefaultTextFontName', 'Arial');
        
        subplot(1, 2, 1);
        plot(x_coords, original, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        hold on;
        plot(x_coords, orig_mean, 'k-', 'LineWidth', 2);
        plot(x_coords, orig_mean + orig_std, 'k--', 'LineWidth', 1);
        plot(x_coords, orig_mean - orig_std, 'k--', 'LineWidth', 1);
        xlabel('Position (pixels)', 'FontName', 'Arial');
        ylabel('Mean Intensity', 'FontName', 'Arial');
        title(sprintf('Before - %s %s', genotypes{g}, channels{ch}), 'FontName', 'Arial');
        legend({'Traces', sprintf('Mean ± SD, N=%d', n_samples)}, 'Location', 'best', 'FontName', 'Arial');
        set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');
        
        subplot(1, 2, 2);
        plot(x_microns, aligned, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        hold on;
        plot(x_microns, aligned_mean, 'k-', 'LineWidth', 2);
        plot(x_microns, aligned_mean + aligned_std, 'k--', 'LineWidth', 1);
        plot(x_microns, aligned_mean - aligned_std, 'k--', 'LineWidth', 1);
        xline(0, 'b--', 'Midline', 'LineWidth', 2);
        xlabel('Position (µm)', 'FontName', 'Arial');
        ylabel('Mean Intensity', 'FontName', 'Arial');
        title(sprintf('After - %s %s', genotypes{g}, channels{ch}), 'FontName', 'Arial');
        legend({'Traces', sprintf('Mean ± SD, N=%d', n_samples)}, 'Location', 'best', 'FontName', 'Arial');
        set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');
        
        saveFigureMultiFormat(fig, fullfile(config.alignmentDir, 'Before_After_Comparison'), ...
                             sprintf('BeforeAfter_%s_%s', channels{ch}, genotypes{g}), config);
    end
end

fprintf('  Before/after plots created\n');
end

function createWTvsKOComparisonPlots(alignedData, statsData, config)
channels = {'Red', 'Green'};
normMethod = 'Linear';
dataTypes = {'original', 'normalized'};
binSizes = {'bin20', 'bin80'};

% Define color schemes: {WT_mean, WT_SD, KO_mean, KO_SD, scheme_name}
colorSchemes = {
    {[0 0 1], [0.6 0.6 1], [1 0 0], [1 0.6 0.6], 'blue_red'};          % WT=blue, KO=red
    {[1 0 0], [1 0.6 0.6], [0 0 1], [0.6 0.6 1], 'red_blue'};          % WT=red, KO=blue
};

for ch = 1:length(channels)
    for dt = 1:length(dataTypes)
        for bs = 1:length(binSizes)
            
            % Get appropriate data based on bin size
            if strcmp(binSizes{bs}, 'bin20')
                if strcmp(dataTypes{dt}, 'original')
                    WT_mean = alignedData.(normMethod).(channels{ch}).WT.mean_profile;
                    WT_std = alignedData.(normMethod).(channels{ch}).WT.std_profile;
                    KO_mean = alignedData.(normMethod).(channels{ch}).KO.mean_profile;
                    KO_std = alignedData.(normMethod).(channels{ch}).KO.std_profile;
                else
                    WT_mean = alignedData.(normMethod).(channels{ch}).WT.mean_profile_01;
                    WT_std = alignedData.(normMethod).(channels{ch}).WT.std_profile_01;
                    KO_mean = alignedData.(normMethod).(channels{ch}).KO.mean_profile_01;
                    KO_std = alignedData.(normMethod).(channels{ch}).KO.std_profile_01;
                end
                x_microns = alignedData.(normMethod).(channels{ch}).WT.x_microns;
                binLabel = 'Bin 20';
            else
                if strcmp(dataTypes{dt}, 'original')
                    WT_mean = alignedData.(normMethod).(channels{ch}).WT.mean_profile_bin80;
                    WT_std = alignedData.(normMethod).(channels{ch}).WT.std_profile_bin80;
                    KO_mean = alignedData.(normMethod).(channels{ch}).KO.mean_profile_bin80;
                    KO_std = alignedData.(normMethod).(channels{ch}).KO.std_profile_bin80;
                else
                    WT_mean = alignedData.(normMethod).(channels{ch}).WT.mean_profile_01_bin80;
                    WT_std = alignedData.(normMethod).(channels{ch}).WT.std_profile_01_bin80;
                    KO_mean = alignedData.(normMethod).(channels{ch}).KO.mean_profile_01_bin80;
                    KO_std = alignedData.(normMethod).(channels{ch}).KO.std_profile_01_bin80;
                end
                x_microns = alignedData.(normMethod).(channels{ch}).WT.x_microns_bin80;
                binLabel = 'Bin 80';
            end
            
            % Get stats (same for both bin sizes since computed on bin20)
            stats_results = statsData.(normMethod).(channels{ch}).(dataTypes{dt});
            num_WT = stats_results.num_WT;
            num_KO = stats_results.num_KO;
            
            % Calculate overall mean and SD for legend
            WT_overall_mean = mean(WT_mean);
            WT_overall_std = mean(WT_std);
            KO_overall_mean = mean(KO_mean);
            KO_overall_std = mean(KO_std);
            
            % Determine significance markers
            neg_marker = getSignificanceMarker(stats_results.p_negative, stats_results.power_negative, config);
            pos_marker = getSignificanceMarker(stats_results.p_positive, stats_results.power_positive, config);
            
            % Generate plot for each color scheme
            for cs = 1:length(colorSchemes)
                WT_color = colorSchemes{cs}{1};
                WT_SD_color = colorSchemes{cs}{2};
                KO_color = colorSchemes{cs}{3};
                KO_SD_color = colorSchemes{cs}{4};
                scheme_name = colorSchemes{cs}{5};
                
                fig = figure('Position', [100, 100, 1200, 700], 'Visible', config.showFigures);
                set(fig, 'DefaultAxesFontName', 'Arial');
                set(fig, 'DefaultTextFontName', 'Arial');
                hold on;
                
                % Shaded SD regions
                fill([x_microns; flipud(x_microns)], [WT_mean + WT_std; flipud(WT_mean - WT_std)], ...
                     WT_SD_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                fill([x_microns; flipud(x_microns)], [KO_mean + KO_std; flipud(KO_mean - KO_std)], ...
                     KO_SD_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
                % Plot means with enhanced legend info
                plot(x_microns, WT_mean, 'Color', WT_color, 'LineWidth', 2, ...
                     'DisplayName', sprintf('WT (Mean=%.2f, SD=%.2f, N=%d)', WT_overall_mean, WT_overall_std, num_WT));
                plot(x_microns, KO_mean, 'Color', KO_color, 'LineWidth', 2, ...
                     'DisplayName', sprintf('KO (Mean=%.2f, SD=%.2f, N=%d)', KO_overall_mean, KO_overall_std, num_KO));
                
                % Midline
                xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                
                xlabel('Position (µm)', 'FontName', 'Arial');
                if strcmp(dataTypes{dt}, 'original')
                    ylabel('Mean Intensity', 'FontName', 'Arial');
                else
                    ylabel('Normalized [0-1]', 'FontName', 'Arial');
                end
                
                % Title with stats
                title(sprintf('%s - %s (%s)\nNegative: p=%.4f, power=%.2f%s | Positive: p=%.4f, power=%.2f%s', ...
                      channels{ch}, dataTypes{dt}, binLabel, ...
                      stats_results.p_negative, stats_results.power_negative, neg_marker, ...
                      stats_results.p_positive, stats_results.power_positive, pos_marker), ...
                      'FontSize', 11, 'FontName', 'Arial');
                legend('Location', 'best', 'FontName', 'Arial');
                set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');

                % Symmetric, padded x-limits centered on the midline (x=0). The
                % SCN midline is right-of-center in the crop, so the raw data
                % extent is asymmetric (e.g. ~-400 to +300) and MATLAB's
                % autoscale would run the curve right up to the edge. Center on
                % 0 and pad ~5% so the whole profile is visible with margin.
                xspan = max(abs(x_microns));
                if isempty(xspan) || ~isfinite(xspan) || xspan <= 0
                    xspan = 1;
                end
                xlim([-1, 1] * xspan * 1.05);

                saveFigureMultiFormat(fig, fullfile(config.alignmentDir, 'Aligned_Profiles'), ...
                                     sprintf('WTvsKO_%s_%s_%s_%s', channels{ch}, dataTypes{dt}, binSizes{bs}, scheme_name), config);
            end
        end
        
        % Create stats summary plot (one per channel/dataType)
        createStatsSummaryPlot(statsData.(normMethod).(channels{ch}).(dataTypes{dt}), ...
                              channels{ch}, dataTypes{dt}, config);
    end
end

fprintf('  WT vs KO plots created (2 color schemes each)\n');
end

function marker = getSignificanceMarker(p, power, config)
% Returns asterisk only if p < alpha AND power > threshold
if p < config.alpha && power > config.powerThreshold
    marker = ' *';
else
    marker = '';
end
end

function createStatsSummaryPlot(stats_results, channelName, dataType, config)
fig = figure('Position', [100, 100, 1000, 600], 'Visible', config.showFigures);
set(fig, 'DefaultAxesFontName', 'Arial');
set(fig, 'DefaultTextFontName', 'Arial');

% Subplot 1: Bar plot of means with error bars
subplot(2, 2, 1);
means = [stats_results.WT_negative_mean, stats_results.KO_negative_mean; ...
         stats_results.WT_positive_mean, stats_results.KO_positive_mean];
stds = [stats_results.WT_negative_std, stats_results.KO_negative_std; ...
        stats_results.WT_positive_std, stats_results.KO_positive_std];
    
b = bar(means, 'grouped');
b(1).FaceColor = [0.2 0.2 0.2];  % WT dark gray
b(2).FaceColor = [0.8 0.2 0.2];  % KO red
hold on;

% Add error bars
ngroups = size(means, 1);
nbars = size(means, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, means(:,i), stds(:,i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
end

xticks(1:2);
xticklabels({'Negative (x<0)', 'Positive (x>0)'});
ylabel('Mean Intensity', 'FontName', 'Arial');
title('Group Means ± SD', 'FontName', 'Arial');
legend({'WT', 'KO'}, 'Location', 'best', 'FontName', 'Arial');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');

% Subplot 2: p-values
subplot(2, 2, 2);
p_vals = [stats_results.p_negative, stats_results.p_positive];
b2 = bar(1:2, p_vals, 0.6);
b2.FaceColor = 'flat';
b2.CData(1,:) = [0.2 0.4 0.8];
b2.CData(2,:) = [0.8 0.4 0.2];
hold on;
yline(config.alpha, 'r--', sprintf('α=%.2f', config.alpha), 'LineWidth', 2);
xticks(1:2);
xticklabels({'Negative', 'Positive'});
ylabel('p-value', 'FontName', 'Arial');
title('Mann-Whitney p-values', 'FontName', 'Arial');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');

% Subplot 3: Effect size (Cohen's d)
subplot(2, 2, 3);
d_vals = [stats_results.cohens_d_negative, stats_results.cohens_d_positive];
b3 = bar(1:2, d_vals, 0.6);
b3.FaceColor = 'flat';
b3.CData(1,:) = [0.2 0.4 0.8];
b3.CData(2,:) = [0.8 0.4 0.2];
hold on;
yline(0.8, 'g--', 'Large (0.8)', 'LineWidth', 1.5);
yline(0.5, 'y--', 'Medium (0.5)', 'LineWidth', 1.5);
yline(0.2, 'c--', 'Small (0.2)', 'LineWidth', 1.5);
xticks(1:2);
xticklabels({'Negative', 'Positive'});
ylabel("Cohen's d", 'FontName', 'Arial');
title('Effect Size', 'FontName', 'Arial');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');

% Subplot 4: Power with significance markers
subplot(2, 2, 4);
power_vals = [stats_results.power_negative, stats_results.power_positive];
b4 = bar(1:2, power_vals, 0.6);
b4.FaceColor = 'flat';
b4.CData(1,:) = [0.2 0.4 0.8];
b4.CData(2,:) = [0.8 0.4 0.2];
hold on;
yline(config.powerThreshold, 'r--', sprintf('Threshold (%.1f)', config.powerThreshold), 'LineWidth', 2);
ylim([0 1.1]);

% Add asterisks if significant AND powered
neg_marker = getSignificanceMarker(stats_results.p_negative, stats_results.power_negative, config);
pos_marker = getSignificanceMarker(stats_results.p_positive, stats_results.power_positive, config);
if ~isempty(neg_marker)
    text(1, power_vals(1) + 0.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Arial');
end
if ~isempty(pos_marker)
    text(2, power_vals(2) + 0.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Arial');
end

xticks(1:2);
xticklabels({'Negative', 'Positive'});
ylabel('Power', 'FontName', 'Arial');
title('Post-hoc Power (* = p<0.05 & power>0.8)', 'FontName', 'Arial');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');

sgtitle(sprintf('Statistical Summary - %s (%s)\nWT N=%d, KO N=%d', ...
        channelName, dataType, stats_results.num_WT, stats_results.num_KO), 'FontSize', 12, 'FontName', 'Arial');

saveFigureMultiFormat(fig, fullfile(config.alignmentDir, 'Statistical_Analysis'), ...
                     sprintf('Stats_%s_%s', channelName, dataType), config);
end

function createValleyDiagnosticPlots(alignedData, allProfileData, config)
channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};
normMethod = 'Linear';

for ch = 1:length(channels)
    for g = 1:length(genotypes)
        profiles = alignedData.(normMethod).(channels{ch}).(genotypes{g}).profiles_align;
        x_coords = alignedData.(normMethod).(channels{ch}).(genotypes{g}).x_coords_align;
        animals = alignedData.(normMethod).(channels{ch}).(genotypes{g}).animals;
        valley_positions = alignedData.(normMethod).(channels{ch}).(genotypes{g}).valley_positions;
        
        x_coords = x_coords(:);
        profiles = ensureColumnProfiles(profiles, length(x_coords));
        
        unique_animals = unique(animals);
        colors = lines(length(unique_animals));
        
        mean_valley_pixel = x_coords(round(mean(valley_positions)));
        
        fig = figure('Position', [100, 100, 1200, 600], 'Visible', config.showFigures);
        set(fig, 'DefaultAxesFontName', 'Arial');
        set(fig, 'DefaultTextFontName', 'Arial');
        hold on;
        
        for a = 1:length(unique_animals)
            animal_idx = find(animals == unique_animals(a));
            for i = 1:length(animal_idx)
                idx = animal_idx(i);
                plot(x_coords, profiles(:, idx), 'Color', colors(a,:), 'LineWidth', 1);
                plot(x_coords(round(valley_positions(idx))), profiles(round(valley_positions(idx)), idx), ...
                     'o', 'Color', colors(a,:), 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
        
        xline(mean_valley_pixel, 'b--', 'LineWidth', 2);
        xline([config.valleySearchMin config.valleySearchMax], 'g:', 'LineWidth', 1.5);
        
        xlabel('Position (pixels)', 'FontName', 'Arial');
        ylabel('Mean Intensity', 'FontName', 'Arial');
        title(sprintf('Valley Detection - %s %s (bin %d)', genotypes{g}, channels{ch}, config.binSize_alignment), 'FontName', 'Arial');
        set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');
        
        saveFigureMultiFormat(fig, fullfile(config.alignmentDir, 'QC_Plots'), ...
                             sprintf('ValleyDetect_%s_%s', channels{ch}, genotypes{g}), config);
    end
end

fprintf('  Valley diagnostic plots created\n');
end

function createAlignmentPlots(alignedData, config)
channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};
normMethod = 'Linear';

for ch = 1:length(channels)
    for g = 1:length(genotypes)
        shifts_pixels = alignedData.(normMethod).(channels{ch}).(genotypes{g}).shifts_pixels;
        deviations = alignedData.(normMethod).(channels{ch}).(genotypes{g}).deviations;
        outliers = alignedData.(normMethod).(channels{ch}).(genotypes{g}).outliers;
        animals = alignedData.(normMethod).(channels{ch}).(genotypes{g}).animals;
        
        fig = figure('Position', [100, 100, 1400, 500], 'Visible', config.showFigures);
        set(fig, 'DefaultAxesFontName', 'Arial');
        set(fig, 'DefaultTextFontName', 'Arial');
        
        subplot(1, 2, 1);
        bar(shifts_pixels);
        hold on;
        if ~isempty(outliers)
            bar(outliers, shifts_pixels(outliers), 'r');
        end
        xticks(1:length(animals));
        xticklabels(arrayfun(@(x) sprintf('%.0f', x), animals, 'UniformOutput', false));
        xtickangle(45);
        xlabel('Animal', 'FontName', 'Arial');
        ylabel('Shift (pixels)', 'FontName', 'Arial');
        title(sprintf('Shifts - %s %s', genotypes{g}, channels{ch}), 'FontName', 'Arial');
        set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');
        
        subplot(1, 2, 2);
        bar(deviations);
        hold on;
        if ~isempty(outliers)
            bar(outliers, deviations(outliers), 'r');
        end
        yline(mean(deviations) + config.outlierThreshold * std(deviations), 'r--', 'LineWidth', 2);
        xticks(1:length(animals));
        xticklabels(arrayfun(@(x) sprintf('%.0f', x), animals, 'UniformOutput', false));
        xtickangle(45);
        xlabel('Animal', 'FontName', 'Arial');
        ylabel('Deviation', 'FontName', 'Arial');
        title('Deviations', 'FontName', 'Arial');
        set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');
        
        saveFigureMultiFormat(fig, fullfile(config.alignmentDir, 'QC_Plots'), ...
                             sprintf('AlignmentQC_%s_%s', channels{ch}, genotypes{g}), config);
        
        createTop5OutlierPlots(alignedData, normMethod, channels{ch}, genotypes{g}, config);
    end
end

fprintf('  QC plots created\n');
end

function createTop5OutlierPlots(alignedData, normMethod, channelName, genotype, config)
aligned = alignedData.(normMethod).(channelName).(genotype).aligned_profiles;
mean_prof = alignedData.(normMethod).(channelName).(genotype).mean_profile;
deviations = alignedData.(normMethod).(channelName).(genotype).deviations;
x_microns = alignedData.(normMethod).(channelName).(genotype).x_microns;
animals = alignedData.(normMethod).(channelName).(genotype).animals;

[~, sorted_idx] = sort(deviations, 'descend');
num_to_plot = min(5, length(deviations));
top_idx = sorted_idx(1:num_to_plot);

fig = figure('Position', [100, 100, 1000, 600], 'Visible', config.showFigures);
set(fig, 'DefaultAxesFontName', 'Arial');
set(fig, 'DefaultTextFontName', 'Arial');

plot(x_microns, mean_prof, 'k-', 'LineWidth', 3, 'DisplayName', 'Mean');
hold on;

colors = lines(num_to_plot);
for i = 1:num_to_plot
    idx = top_idx(i);
    plot(x_microns, aligned(:, idx), 'LineWidth', 1.5, 'Color', colors(i,:), ...
         'DisplayName', sprintf('Animal %.0f', animals(idx)));
end

xline(0, 'b--', 'LineWidth', 1.5);

xlabel('Position (µm)', 'FontName', 'Arial');
ylabel('Mean Intensity', 'FontName', 'Arial');
title(sprintf('Top %d Outliers - %s %s', num_to_plot, genotype, channelName), 'FontName', 'Arial');
legend('show', 'Location', 'best', 'FontName', 'Arial');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontName', 'Arial');

saveFigureMultiFormat(fig, fullfile(config.alignmentDir, 'Outliers'), ...
                     sprintf('Top5_%s_%s', channelName, genotype), config);
end

function saveShiftedImages(alignedData, imageInfo, config)
if ~exist(char(config.normalizedDataFile), 'file')
    warning('normalized_data.mat not found! Skipping image shifting.');
    return;
end

loaded = load(char(config.normalizedDataFile));
linearData = loaded.data.linearData;

channels = {'Red', 'Green'};
genotypes = {'WT', 'KO'};
normMethod = 'Linear';

totalImages = 0;

for ch = 1:length(channels)
    original_images = linearData.(channels{ch}).images;
    
    for g = 1:length(genotypes)
        shifts_pixels = alignedData.(normMethod).(channels{ch}).(genotypes{g}).shifts_pixels;
        animals = alignedData.(normMethod).(channels{ch}).(genotypes{g}).animals;
        
        genotype_idx = find(strcmp({imageInfo.genotype}, genotypes{g}));
        
        for i = 1:length(genotype_idx)
            img_idx = genotype_idx(i);
            animal_id = imageInfo(img_idx).animal;
            
            animal_shift_idx = find(animals == animal_id);
            if isempty(animal_shift_idx), continue; end
            
            % Image shifting disabled (see createAverageImages note). Save the
            % aligned-set images without any X-shift; filename says 'unshifted'
            % so it doesn't misrepresent the contents.
            original_img = original_images(:, :, img_idx);
            shifted_img = original_img;   % no shift (was applyXShiftWithPadding)
            
            filename = sprintf('%s_%s_Animal%d_unshifted.%s', genotypes{g}, channels{ch}, animal_id, config.imageFormat);
            outPath = char(fullfile(config.shiftedImagesDir, channels{ch}, genotypes{g}, filename));
            
            imwrite(uint16(shifted_img), outPath);
            totalImages = totalImages + 1;
        end
    end
end

fprintf('  Saved %d shifted images\n', totalImages);
end

function shifted_img = applyXShiftWithPadding(img, shift_pixels)
[height, width] = size(img);

if shift_pixels == 0
    shifted_img = img;
    return;
end

shifted_img = zeros(height, width, 'like', img);

if shift_pixels > 0
    if shift_pixels >= width
        shifted_img = repmat(img(:, 1), 1, width);
    else
        shifted_img(:, shift_pixels+1:end) = img(:, 1:width-shift_pixels);
        shifted_img(:, 1:shift_pixels) = repmat(img(:, 1), 1, shift_pixels);
    end
else
    shift_pixels = abs(shift_pixels);
    if shift_pixels >= width
        shifted_img = repmat(img(:, end), 1, width);
    else
        shifted_img(:, 1:width-shift_pixels) = img(:, shift_pixels+1:end);
        shifted_img(:, width-shift_pixels+1:end) = repmat(img(:, end), 1, shift_pixels);
    end
end
end

function saveFigureMultiFormat(fig, outputDir, baseName, config)
% Save figure in multiple formats
% EPS files are saved in RGB colorspace
outputDir = char(outputDir);
baseName = char(baseName);

% Verify figure is valid
if ~ishandle(fig) || ~isgraphics(fig, 'figure')
    warning('Invalid figure handle for %s, skipping save.', baseName);
    return;
end

% Ensure output directory exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Ensure figure is fully rendered before saving
drawnow;

if config.savePNG
    pngFile = fullfile(outputDir, [baseName '.png']);
    print(fig, pngFile, '-dpng', sprintf('-r%d', config.figureResolution));
end

if config.saveFIG
    figFile = fullfile(outputDir, [baseName '.fig']);
    savefig(fig, figFile);
end

if config.saveEPS
    epsFile = fullfile(outputDir, [baseName '.eps']);
    % Try exportgraphics first (R2020a+) with RGB colorspace
    try
        exportgraphics(fig, epsFile, 'ContentType', 'vector', 'BackgroundColor', 'white', 'ColorSpace', 'rgb');
    catch
        try
            exportgraphics(fig, epsFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
        catch
            % Fallback: use print with painters renderer (RGB by default, no -cmyk flag)
            try
                print(fig, epsFile, '-depsc2', '-painters', sprintf('-r%d', config.figureResolution));
            catch ME
                warning('EPS export failed for %s (PNG/FIG still saved OK): %s', epsFile, ME.message);
            end
        end
    end
end

if ~config.showFigures || config.closeFIG
    try
        close(fig);
    catch
        % Figure handle already invalid (can happen over remote desktop) - nothing to close
    end
end
end

function saveAlignedDataWithStats(alignedData, statsData, imageInfo, config)
data.alignedData = alignedData;
data.statsData = statsData;
data.imageInfo = imageInfo;
data.config = config;

data.metadata = struct(...
    'alignmentMethod', sprintf('Valley-based (bin %d)', config.binSize_alignment), ...
    'plottingBinSize', config.binSize_analysis, ...
    'pixelSize', config.pixelSize, ...
    'binSize_alignment', config.binSize_alignment, ...
    'binSize_analysis', config.binSize_analysis, ...
    'statisticalTest', 'Mann-Whitney U test on per-animal means');

save(char(config.alignedDataFile), 'data', '-v7.3');
fprintf('  Saved: %s\n', config.alignedDataFile);

createSummaryFiles(data, config);
end

function createSummaryFiles(data, config)
channels = {'Red', 'Green'};
normMethod = 'Linear';
dataTypes = {'original', 'normalized'};

for ch = 1:length(channels)
    for dt = 1:length(dataTypes)
        stats_results = data.statsData.(normMethod).(channels{ch}).(dataTypes{dt});
        
        if strcmp(dataTypes{dt}, 'original')
            WT_mean = data.alignedData.(normMethod).(channels{ch}).WT.mean_profile;
            WT_std = data.alignedData.(normMethod).(channels{ch}).WT.std_profile;
            KO_mean = data.alignedData.(normMethod).(channels{ch}).KO.mean_profile;
            KO_std = data.alignedData.(normMethod).(channels{ch}).KO.std_profile;
        else
            WT_mean = data.alignedData.(normMethod).(channels{ch}).WT.mean_profile_01;
            WT_std = data.alignedData.(normMethod).(channels{ch}).WT.std_profile_01;
            KO_mean = data.alignedData.(normMethod).(channels{ch}).KO.mean_profile_01;
            KO_std = data.alignedData.(normMethod).(channels{ch}).KO.std_profile_01;
        end
        
        x_microns = data.alignedData.(normMethod).(channels{ch}).WT.x_microns;
        
        % Save profile data CSV
        T = table(x_microns, WT_mean, WT_std, KO_mean, KO_std, ...
                  'VariableNames', {'Position_um', 'WT_Mean', 'WT_SD', 'KO_Mean', 'KO_SD'});
        
        csvFile = char(fullfile(config.alignmentDir, 'Statistical_Analysis', ...
                           sprintf('Profiles_%s_%s.csv', channels{ch}, dataTypes{dt})));
        writetable(T, csvFile);
        
        % Save statistical test results
        stats_file = char(fullfile(config.alignmentDir, 'Statistical_Analysis', ...
                                   sprintf('MannWhitney_%s_%s.txt', channels{ch}, dataTypes{dt})));
        fid = fopen(stats_file, 'w');
        fprintf(fid, 'Mann-Whitney U Test Results - %s (%s)\n', channels{ch}, dataTypes{dt});
        fprintf(fid, '================================================\n\n');
        fprintf(fid, 'Test: Mann-Whitney U (Wilcoxon rank-sum) test\n');
        fprintf(fid, 'Comparison: Per-animal mean intensities for each half\n\n');
        fprintf(fid, 'Sample sizes (biological replicates): WT N=%d, KO N=%d\n\n', ...
                stats_results.num_WT, stats_results.num_KO);
        
        fprintf(fid, 'NEGATIVE HALF (x < 0, %s):\n', lateralityFor(channels{ch}, 'negative', config));
        fprintf(fid, '  WT: Mean=%.4f, SD=%.4f\n', stats_results.WT_negative_mean, stats_results.WT_negative_std);
        fprintf(fid, '  KO: Mean=%.4f, SD=%.4f\n', stats_results.KO_negative_mean, stats_results.KO_negative_std);
        fprintf(fid, '  U statistic: %.1f\n', stats_results.U_negative);
        fprintf(fid, '  p-value: %.6f\n', stats_results.p_negative);
        fprintf(fid, "  Cohen's d: %.4f\n", stats_results.cohens_d_negative);
        fprintf(fid, '  Post-hoc power: %.4f\n', stats_results.power_negative);
        neg_sig = stats_results.p_negative < config.alpha && stats_results.power_negative > config.powerThreshold;
        fprintf(fid, '  Significant (p<%.2f & power>%.1f): %s\n\n', config.alpha, config.powerThreshold, mat2str(neg_sig));
        
        fprintf(fid, 'POSITIVE HALF (x > 0, %s):\n', lateralityFor(channels{ch}, 'positive', config));
        fprintf(fid, '  WT: Mean=%.4f, SD=%.4f\n', stats_results.WT_positive_mean, stats_results.WT_positive_std);
        fprintf(fid, '  KO: Mean=%.4f, SD=%.4f\n', stats_results.KO_positive_mean, stats_results.KO_positive_std);
        fprintf(fid, '  U statistic: %.1f\n', stats_results.U_positive);
        fprintf(fid, '  p-value: %.6f\n', stats_results.p_positive);
        fprintf(fid, "  Cohen's d: %.4f\n", stats_results.cohens_d_positive);
        fprintf(fid, '  Post-hoc power: %.4f\n', stats_results.power_positive);
        pos_sig = stats_results.p_positive < config.alpha && stats_results.power_positive > config.powerThreshold;
        fprintf(fid, '  Significant (p<%.2f & power>%.1f): %s\n', config.alpha, config.powerThreshold, mat2str(pos_sig));
        
        fclose(fid);
    end
end

fprintf('  CSV and statistics summaries saved\n');
end
