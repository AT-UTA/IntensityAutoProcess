function sigma = optimize_sigma(img)
% Define a range of sigma values
sigma_range = 0.1:0.1:3;
% Initialize variables to keep track of the best sigma value and its corresponding score
best_sigma = sigma_range(1);
best_score = -Inf;
% Loop over all sigma values in the range
for sigma = sigma_range
% Apply the Gaussian filter with the current sigma value
filtered_img = imgaussfilt(img, sigma);
% Calculate the entropy of the filtered image
entropy_score = entropy(filtered_img);
% If the entropy score is higher than the best score found so far, update the best score and the best sigma value
if entropy_score > best_score
best_score = entropy_score;
best_sigma = sigma;
end
end
% Return the best sigma value
sigma = best_sigma;
end