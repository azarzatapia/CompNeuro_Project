
% Plot firing rate dynamics
figure;
plot(T, rE, 'b', 'LineWidth', 1.5); hold on;
plot(T, rI, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Firing Rate (a.u.)', 'FontSize', 12);
legend('r_E (Excitatory)', 'r_I (Inhibitory)', 'Location', 'Best', 'FontSize', 10);
title('Experiment 4: Firing Rate Model (Sigmoid Wilson-Cowan)', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12, 'Box', 'off');
xlim([T(1), T(end)]);
