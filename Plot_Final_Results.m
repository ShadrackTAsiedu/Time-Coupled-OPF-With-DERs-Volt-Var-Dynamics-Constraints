% Data from your table
DER_limit = [290 300 335 350 400];              % (A)
DER_inj   = [290 300 335 335 335];              % (A)
Vrec_pu   = [0.9436 0.9482 0.9550 0.9550 0.9550];  % (pu)

% Plot
figure('Color','w'); hold on; grid on
set(gcf,'Position',[100 100 820 420])

% Left axis: DER injection
yyaxis left
plot(DER_limit, DER_inj, '-o', 'LineWidth', 2, 'MarkerSize', 6)
ylabel('DER Injection (A)', 'FontSize', 14)
ylim([0, max(DER_inj)*1.1])
ax = gca; 
ax.YColor = 'k';              % left y-axis color (ticks & label) to black

% Right axis: Voltage recovery
yyaxis right
plot(DER_limit, Vrec_pu, '-s', 'LineWidth', 2, 'MarkerSize', 6)
ylabel('Recovered Voltage (pu)', 'FontSize', 14)
ylim([min(Vrec_pu)*0.98, 1.00])
yline(0.95,'k--','0.95 pu','LineWidth',1,'LabelHorizontalAlignment','left');
ax = gca;
ax.YColor = 'k';              % right y-axis color (ticks & label) to black

% Common settings
xlabel('DER Limit (A)', 'FontSize', 12)
xticks(DER_limit)
legend({'DER Injection','Recovered Voltage'}, 'Location','southeast')

ax.FontSize=14;
box on