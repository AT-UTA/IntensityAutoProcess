function [averagedc1_intensities, averagedc2_intensities, averagedc3_intensities] = process_cell_imagev2(mainfolder)
   
%Preallocating variables to store the intensity data
    c1_intensities = [];
    c2_intensities = [];
    c3_intensities = [];
    collapsedc1 = [];
    collapsedc2 = [];
    collapsedc3 = [];
    averagedc1_intensities = [];
    averagedc2_intensities = [];
    averagedc3_intensities = [];
    subfolders = dir(mainfolder);
    subfolders = subfolders(~ismember({subfolders.name},{'.','..'}));
    groupnum = length(subfolders);

% Ask the user for the names of the channels
    channel1_name = input('Please enter the name for channel 1: ', 's');
    channel2_name = input('Please enter the name for channel 2: ', 's');
    channel3_name = input('Please enter the name for channel 3: ', 's');

for k = 1:groupnum
   subfolders_path = fullfile(mainfolder, subfolders(k).name);
   %getting images from group
   % Get a list of all the .nd2 files in the current subfolder
   nd2_files = dir(fullfile(subfolders_path, '*.nd2'));
   % Load the images using bioformats
   num_images = length(nd2_files);  
 
   for i = 1:num_images
        c1_intensity = [];
        c2_intensity = [];
        c3_intensity = [];
        nd2_files_path = fullfile(subfolders_path, nd2_files(i).name);
        reader = bfGetReader(nd2_files_path);
        channel_num = reader.getSizeC();
        num_z_planes = reader.getSizeZ();
        for j = 1:num_z_planes
            c1 = bfGetPlane(reader, (j - 1) * channel_num + 1);   
            c2 = bfGetPlane(reader, (j - 1) * channel_num + 2);
            c3 = bfGetPlane(reader, (j - 1) * channel_num + 3);

            % Apply Gaussian filter
            c1 = imgaussfilt(c1, optimize_sigma(c1));
            c2 = imgaussfilt(c2, optimize_sigma(c2));
            c3 = imgaussfilt(c3, optimize_sigma(c3));

            image(:,:,1) = c1;
            image(:,:,2) = c2;
            image(:,:,3) = c3;
            % Normalize the pixel intensities for each channel with respect to the global intensity of that channel
            c1_intensity(:,:,j) = mat2gray(image(:,:,1));
            c2_intensity(:,:,j) = mat2gray(image(:,:,2));
            c3_intensity(:,:,j) = mat2gray(image(:,:,3));
        end
        %Sum of Z-plane intensities for each channel and image into a few
        %master variables
        collapsedc1(:,:,i) = sum(c1_intensity,3);
        collapsedc2(:,:,i) = sum(c2_intensity,3);
        collapsedc3(:,:,i) = sum(c3_intensity,3);
       
   end
   averagedc1_intensities(:,:,k) = (sum(collapsedc1,3) / num_images);
   averagedc2_intensities(:,:,k) = (sum(collapsedc2,3) / num_images);
   averagedc3_intensities(:,:,k) = (sum(collapsedc3,3) / num_images);
   % Append the intensity data for each channel to the preallocated variables
   c1_intensities = [c1_intensities, reshape(averagedc1_intensities(:,:,k), [], 1)];
   c2_intensities = [c2_intensities, reshape(averagedc2_intensities(:,:,k), [], 1)];
   c3_intensities = [c3_intensities, reshape(averagedc3_intensities(:,:,k), [], 1)];
   end

% Perform ANOVA test for each channel
[p_c1, tbl_c1, stats_c1] = anova1([c1_intensities], [], 'off');
[p_c2, tbl_c2, stats_c2] = anova1([c2_intensities], [], 'off');
[p_c3, tbl_c3, stats_c3] = anova1([c3_intensities], [], 'off');

% Determine which image had the highest intensity for each channel
[max_c1, idx_c1] = max(sum( averagedc1_intensities, [1, 2]));
[max_c2, idx_c2] = max(sum( averagedc2_intensities, [1, 2]));
[max_c3, idx_c3] = max(sum( averagedc3_intensities, [1, 2]));

disp([channel1_name ': ' subfolders(idx_c1).name ' had the highest intensity (p-value = ' num2str(p_c1) ')']);
disp([channel2_name ': ' subfolders(idx_c2).name ' had the highest intensity (p-value = ' num2str(p_c2) ')']);
disp([channel3_name ': ' subfolders(idx_c3).name ' had the highest intensity (p-value = ' num2str(p_c3) ')']);

% Perform the Tukey-Kramer post-hoc test for each channel
% Channel 1
multcompare_c1 = multcompare(stats_c1, 'Alpha', 0.05, 'CType', 'tukey-kramer');

% Channel 2
multcompare_c2 = multcompare(stats_c2, 'Alpha', 0.05, 'CType', 'tukey-kramer');

% Channel 3
multcompare_c3 = multcompare(stats_c3, 'Alpha', 0.05, 'CType', 'tukey-kramer');

% Display the results of the post-hoc tests for each channel
% Channel 1
fprintf('Post-hoc tests for %s:\n', channel1_name);
fprintf('%-8s%-8s%-12s%-12s%-12s%s\n', 'Group 1', 'Group 2', 'Diff', 'Lower', 'Upper', 'p-Value');
for i = 1:size(multcompare_c1, 1)
    fprintf('%-8s%-8s%-12.4f%-12.4f%-12.4f%.4f\n', subfolders(multcompare_c1(i,1)).name, subfolders(multcompare_c1(i,2)).name, multcompare_c1(i,3:end));
end

% Channel 2
fprintf('\nPost-hoc tests for %s:\n', channel2_name);
fprintf('%-8s%-8s%-12s%-12s%-12s%s\n', 'Group 1', 'Group 2', 'Diff', 'Lower', 'Upper', 'p-Value');
for i = 1:size(multcompare_c2, 1)
    fprintf('%-8s%-8s%-12.4f%-12.4f%-12.4f%.4f\n', subfolders(multcompare_c2(i,1)).name, subfolders(multcompare_c2(i,2)).name, multcompare_c2(i,3:end));
end

% Channel 3
fprintf('\nPost-hoc tests for %s:\n', channel3_name);
fprintf('%-8s%-8s%-12s%-12s%-12s%s\n', 'Group 1', 'Group 2', 'Diff', 'Lower', 'Upper', 'p-Value');
for i = 1:size(multcompare_c3, 1)
    fprintf('%-8s%-8s%-12.4f%-12.4f%-12.4f%.4f\n', subfolders(multcompare_c3(i,1)).name, subfolders(multcompare_c3(i,2)).name, multcompare_c3(i,3:end));
end

% Calculate mean and standard deviation for each group and channel
averagedc1_intensities_mean = squeeze(mean(averagedc1_intensities, [1, 2]));
averagedc1_intensities_std = squeeze(std(averagedc1_intensities, [], [1, 2]));

averagedc2_intensities_mean = squeeze(mean(averagedc2_intensities, [1, 2]));
averagedc2_intensities_std = squeeze(std(averagedc2_intensities, [], [1, 2]));

averagedc3_intensities_mean = squeeze(mean(averagedc3_intensities, [1, 2]));
averagedc3_intensities_std = squeeze(std(averagedc3_intensities, [], [1, 2]));

% Display mean and standard deviation for each group and channel
disp('Mean and standard deviation for each channel and folder:');
for m = 1:numel(subfolders)
    disp([subfolders(m).name ':']);
    disp([channel1_name ' - Mean: ' num2str(averagedc1_intensities_mean(m)) '  Std: ' num2str(averagedc1_intensities_std(m))]);
    disp([channel2_name ' - Mean: ' num2str(averagedc2_intensities_mean(m)) '  Std: ' num2str(averagedc2_intensities_std(m))]);
    disp([channel3_name ' - Mean: ' num2str(averagedc3_intensities_mean(m)) '  Std: ' num2str(averagedc3_intensities_std(m))]);
end
% Create a bar graph with mean intensities for each group and channel
figure;
hold on;
bar_width = 0.2;
group_labels = {subfolders.name};
group_x_positions = 1:groupnum;
c1_x_positions = group_x_positions - bar_width;
c2_x_positions = group_x_positions;
c3_x_positions = group_x_positions + bar_width;

% Create a bar for each channel
c1_bar = bar(c1_x_positions, averagedc1_intensities_mean, bar_width, 'FaceColor', 'b');
c2_bar = bar(c2_x_positions, averagedc2_intensities_mean, bar_width, 'FaceColor', 'g');
c3_bar = bar(c3_x_positions, averagedc3_intensities_mean, bar_width, 'FaceColor', 'r');

% Add error bars to the graph
errorbar(c1_x_positions, averagedc1_intensities_mean, averagedc1_intensities_std, 'k', 'LineStyle', 'none');
errorbar(c2_x_positions, averagedc2_intensities_mean, averagedc2_intensities_std, 'k', 'LineStyle', 'none');
errorbar(c3_x_positions, averagedc3_intensities_mean, averagedc3_intensities_std, 'k', 'LineStyle', 'none');

% Add x-axis and y-axis labels, and a title to the graph
xlabel('Groups');
ylabel('Mean intensity');
title('Mean intensity for each group and channel');

% Set the x-tick positions and labels
xticks(group_x_positions);
xticklabels(group_labels);

% Add a legend to the graph
legend([c1_bar, c2_bar, c3_bar], {channel1_name, channel2_name, channel3_name}, 'Location', 'northwest');

% Box plot for channel 1
figure
boxplot(c1_intensities)
title(['Box Plot for ' channel1_name])
xlabel('Folder Name')
xticklabels({subfolders.name})
ylabel('Intensity')

% Box plot for channel 2
figure
boxplot(c2_intensities)
title(['Box Plot for ' channel2_name])
xlabel('Folder Name')
xticklabels({subfolders.name})
ylabel('Intensity')

% Box plot for channel 3
figure
boxplot(c3_intensities)
title(['Box Plot for ' channel3_name])
xlabel('Folder Name')
xticklabels({subfolders.name})
ylabel('Intensity')
figure
    %Histogram of intensities for each channel
    for l = 1:groupnum
        subplot(groupnum, 1, l)
        histogram(c1_intensities(:,l));
        title(['Histogram for ' subfolders(l).name ' in ' channel1_name])
        xlabel('Intensity')
        ylabel('Frequency')
    end
    figure
    for l = 1:groupnum
        subplot(groupnum, 1, l)
        histogram(c2_intensities(:,l));
        title(['Histogram for ' subfolders(l).name ' in ' channel2_name])
        xlabel('Intensity')
        ylabel('Frequency')
    end
    figure
    for l = 1:groupnum
        subplot(groupnum, 1, l)
        histogram(c3_intensities(:,l));
        title(['Histogram for ' subfolders(l).name ' in ' channel3_name])
        xlabel('Intensity')
        ylabel('Frequency')
    end
end