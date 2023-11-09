% Adds labels to a plot from TDTR_Bidirectional_MAIN_FIT.m
% change n in figure(n) to the correct value

figure(10)
xlabel('\it{t}_{d} (s)')
ylabel('-{\it V}_{in}/{\it V}_{out} (a.u.)')
%title(['CuI 80deg 1 min, {\it \lambda}_{CuI} = ', num2str(Lambda(3), '%.3g,'), '{\it C}_{Al} = ', num2str(C(1), '%.3g'), '{\it C}_{CuI} = ', num2str(C(3), '%.3g')])
%title(['Si, {\it h}_{CuI} =', num2str(h(3), ' %.3g')])
%title(['Al / Int / CuI / Int / Si, {\it \lambda}_{CuI} = ', num2str(Lambda(3), '%.3g,'), '{\it C}_{CuI} = ', num2str(C(3), '%.3g')])
title(['Au / Ti / CuI / Int / Si, {\it h}_{Au} = ', num2str(h(1), '%.3g,'), '{\it h}_{Ti} = ', num2str(h(2), '%.3g'), '{\it h}_{CuI} = ', num2str(h(3), '%.3g')])
fontsize(11,'points')
%legend('data','fit')

% Outputs a pdf of the figure
fig=gcf;
exportgraphics(fig, 'output.pdf', 'ContentType', 'vector');