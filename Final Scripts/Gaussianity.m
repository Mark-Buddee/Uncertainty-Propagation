% Example data
val = randn(1000,6);
valRel = randn(1000,6);
valRoe = randn(1000,6);

% Plot Cartesian
figure;
for i = 1:6
    subplot(2, 3, i);
    qqplot(val(:, i));
    plot_latex([], 'Theoretical Quantiles', 'Sample Quantiles', '', ...
               sprintf('Q-Q Plot: Dim %d', i), '');
end
% sgtitle('Cartesian', 'FontWeight', 'bold');

% Plot Hillframe
figure;
for i = 1:6
    subplot(2, 3, i);
    qqplot(valRel(:, i));
    plot_latex([], 'Theoretical Quantiles', 'Sample Quantiles', '', ...
               sprintf('Q-Q Plot: Dim %d', i), '');
end
% sgtitle('Hillframe', 'FontWeight', 'bold');
% 
% Plot ROE
figure;
for i = 1:6
    subplot(2, 3, i);
    qqplot(valRoe(:, i));
    plot_latex([], 'Theoretical Quantiles', 'Sample Quantiles', '', ...
               sprintf('Q-Q Plot: Dim %d', i), '');
end
% sgtitle('ROE', 'FontWeight', 'bold');
