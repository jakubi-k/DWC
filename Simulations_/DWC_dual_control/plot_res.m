figure;

sgtitle("Dual side control with DWC")

subplot(2,1,1);
plot(dwc_dual{6}.Values);
title("Mutual inductance (M) over time")
ylabel("Inductance (H)")
grid on;
xlim([0 0.25])

subplot(2,1,2);
plot(dwc_dual{7}.Values);
title("Output current over time")
ylabel("Current (A)")
grid on;
xlim([0 0.25])
