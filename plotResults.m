%% Helper function for plotting
function plotResults(xaxis, sparsity_type, NMSE_DOL, NMSE_MSBL, NMSE_aditia, NMSE_PC, NMSE_Oracle_MMSE, ...
    SR_D, SR_SBL, SR_ad, SR_PC)

    % Set figure properties for publication quality
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    
    % Create color schemes and markers for consistency
    colors = lines(5);
    markers = {'o', 's', 'd', '^', 'x'};
    line_styles = {'--', ':', '--', '--', '-'};
    
    % Set figure size for publication
    figure_width = 7;  % inches
    figure_height = 5; % inches
    
    % Plot NMSE vs SNR
    fig1 = figure(1);
    set(fig1, 'Units', 'inches', 'Position', [1, 1, figure_width, figure_height]);
    clf;
    hold on;
    
    % Create plots with consistent styling
    h1 = plot(xaxis, NMSE_DOL, [markers{1} line_styles{1}], 'Color', colors(1,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h2 = plot(xaxis, NMSE_MSBL, [markers{2} line_styles{2}], 'Color', colors(2,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h3 = plot(xaxis, NMSE_aditia, [markers{3} line_styles{3}], 'Color', colors(3,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h4 = plot(xaxis, NMSE_PC, [markers{4} line_styles{4}], 'Color', colors(4,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h5 = plot(xaxis, NMSE_Oracle_MMSE, [markers{5} line_styles{5}], 'Color', colors(5,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    
    % Customize appearance
    grid on;
    box on;
    title_text = [sparsity_type ' Sparse Signal Recovery'];
    
    % Create legend
    legend([h1, h2, h3, h4, h5], 'Proposed', 'M-SBL', 'MSBL-DoL', 'PC-MSBL', 'Oracle MMSE', ...
        'FontSize', 10, 'Location', 'best', 'Interpreter', 'latex');
    
    % Labels
    xlabel('SNR [dB]', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('NMSE [dB]', 'FontSize', 12, 'Interpreter', 'latex');
    title(title_text, 'FontSize', 14, 'Interpreter', 'latex');
    
    % Adjust axes
    ax = gca;
    ax.FontSize = 10;
    ax.TickLabelInterpreter = 'latex';
    ylim([min([NMSE_DOL, NMSE_MSBL, NMSE_aditia, NMSE_PC, NMSE_Oracle_MMSE]) - 5, ...
          max([NMSE_DOL, NMSE_MSBL, NMSE_aditia, NMSE_PC]) + 5]);
    
    % Save figure
    saveas(fig1, sprintf('NMSE_%s_sparsity.fig', lower(sparsity_type)));
    print(fig1, sprintf('NMSE_%s_sparsity.pdf', lower(sparsity_type)), '-dpdf', '-r300');
    
    % Plot F1 Score vs SNR
    fig2 = figure(2);
    set(fig2, 'Units', 'inches', 'Position', [9, 1, figure_width, figure_height]);
    clf;
    hold on;
    
    % Create plots with consistent styling
    h1 = plot(xaxis, SR_D, [markers{1} line_styles{1}], 'Color', colors(1,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h2 = plot(xaxis, SR_SBL, [markers{2} line_styles{2}], 'Color', colors(2,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h3 = plot(xaxis, SR_ad, [markers{3} line_styles{3}], 'Color', colors(3,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    h4 = plot(xaxis, SR_PC, [markers{4} line_styles{4}], 'Color', colors(4,:), 'MarkerSize', 8, 'LineWidth', 2.5);
    
    % Customize appearance
    grid on;
    box on;
    
    % Create legend
    legend([h1, h2, h3, h4], 'Proposed', 'M-SBL', 'MSBL-DoL', 'PC-MSBL', ...
        'FontSize', 10, 'Location', 'best', 'Interpreter', 'latex');
    
    % Labels
    xlabel('SNR [dB]', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('F1 Score', 'FontSize', 12, 'Interpreter', 'latex');
    title(title_text, 'FontSize', 14, 'Interpreter', 'latex');
    
    % Adjust axes
    ax = gca;
    ax.FontSize = 10;
    ax.TickLabelInterpreter = 'latex';
    ylim([0, 1.05]);
    
    % Save figure
 %   saveas(fig2, sprintf('F1Score_%s_sparsity.fig', lower(sparsity_type)));
  %  print(fig2, sprintf('F1Score_%s_sparsity.pdf', lower(sparsity_type)), '-dpdf', '-r300');
end