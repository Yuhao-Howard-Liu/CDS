plot(0:0.0005:0.0025,ccpdft0.*10000,'+-')
grid on 
xlabel("\theta_{\lambda_2}(0)")
ylabel("Upfront/bps")

plot(0.0025:0.0015:0.01,ccpdft1.*10000,'+-')
grid on 
xlabel("\theta_{\lambda_2}(1)")
ylabel("Upfront/bps")


x_points = [0.65, 0.65, 1.5, 1.5];  
y_points = [-0.02, 0.14, 0.14, -0.02];
color = [0, 0, 1];
hold on;
a = fill(x_points, y_points, color,'EdgeColor', 'none');
a.FaceAlpha = 0.1;
plot(0:0.05:2,upft,'+-')
grid on
xlabel("Time")
ylabel("Upfront")





[X,Y]=meshgrid([0.85 0.87 0.89 0.91 0.93 0.95],[0.75 0.77 0.79 0.81 0.83 0.85]);
surf(X,Y,rfdft)
colorbar
xlabel("E^Q[\psi_{1,0}]")
ylabel("E^Q[\psi_{1,1}]")
zlabel("Upfront")

[X,Y]=meshgrid([0 0.0004 0.0008 0.0012 0.0016 0.002],[0.002 0.003 0.004 0.005 0.006 0.007]);
surf(X,Y,ccpdft)
colorbar
xlabel("\theta_{\lambda_2}(0)")
ylabel("\theta_{\lambda_2}(1)")
zlabel("Upfront")