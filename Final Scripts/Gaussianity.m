% for i = 1: nPoints  
%     val(i, :)
%     valRel(i, :)
%     valRoe(i, :)
% end

data = randn(1000,1);  % Normally distributed example

figure;
for i = 1:6
    subplot(2, 3, i);              % 2 rows × 3 columns of subplots
    qqplot(val(:, i));              % Q-Q plot for each axis
    title(sprintf('Q-Q Plot: Dim %d', i));
    xlabel('Theoretical Quantiles');
    ylabel('Sample Quantiles');
    grid on;
end

sgtitle('Cartesian');

figure;
for i = 1:6
    subplot(2, 3, i);              % 2 rows × 3 columns of subplots
    qqplot(valRel(:, i));              % Q-Q plot for each axis
    title(sprintf('Q-Q Plot: Dim %d', i));
    xlabel('Theoretical Quantiles');
    ylabel('Sample Quantiles');
    grid on;
end

sgtitle('Hillframe');

figure;
for i = 1:6
    subplot(2, 3, i);              % 2 rows × 3 columns of subplots
    qqplot(valRoe(:, i));              % Q-Q plot for each axis
    title(sprintf('Q-Q Plot: Dim %d', i));
    xlabel('Theoretical Quantiles');
    ylabel('Sample Quantiles');
    grid on;
end

sgtitle('ROE');