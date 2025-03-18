function [F1] = F1_score(est_x, true_supp, thr)
    % Calculate the F1 score between estimated support and true support.
    %
    % Parameters:
    % est_x - Estimated signals (matrix)
    % true_supp - True support indices (vector)
    % thr - Threshold for determining significant entries in est_x
    %
    % Returns:
    % F1 - Calculated F1 score

    % Find the indices of estimated support based on the threshold
    estimated_support = find(vecnorm(est_x, 2, 2) > thr);
   % [~,I]=sort(vecnorm(est_x, 2, 2),'descend');
%estimated_support = I(1:length(true_supp));
% Calculate true positives (TP)
    true_positives = length(intersect(estimated_support, true_supp));

    % Handle the case when there are no true positives
    if true_positives == 0
        F1 = 0;
        return; % Exit the function early
    end

    % Calculate false negatives (FN) and false positives (FP)
    false_negatives = length(setdiff(true_supp, estimated_support));
    false_positives = length(setdiff(estimated_support, true_supp));

    % Calculate precision and recall
    precision = true_positives / (true_positives + false_positives);
    recall = true_positives / (true_positives + false_negatives);

    % Calculate F1 score, handling division by zero
    if (precision + recall) == 0
        F1 = 0; % If both precision and recall are zero, set F1 to zero
    else
        F1 = 2 * (precision * recall) / (precision + recall);
    end
end