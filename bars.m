% x = ["$(\mathbf{x}_\mathrm{Q}, \mathbf{x})$" "Sparsity-based" "PHADQ-B" 'PHADQ-U' 'PHADQ-B (oracle)'];
% y = [ODGq CP_odg B_odg U_odg oracle_odg];
% 
% 
% % Create the bar plot
% figure;
% b = bar(y);
% 
% xtips1 = b(1).XEndPoints;
% ytips1 = b(1).YEndPoints;
% labels1 = string(b(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% b.FaceColor = 'flat';
% b.CData(2,:) = [.5 0 .5];
% b.CData(3,:) = [0 0 .5];
% b.CData(1,:) = [.5 0 0];
% 
% ax = gca;  % Get current axis
% ax.XTickLabel = x;  % Set the X-tick labels
% ax.XTickLabelRotation = 45;  % Optional: Rotate X-tick labels for better visibility
% 
% % Set the interpreter for the X-axis labels to LaTeX
% ax.XAxis.TickLabelInterpreter = 'latex';
% title('ODG (6\,bps)', 'Interpreter','latex');

figure;
plot(x(48000:48500));
hold on;
plot(xcp(48000:48500));
plot(x_hat(48000:48500));
legend('original', 'Sparsity-based', 'PHADQ-B')