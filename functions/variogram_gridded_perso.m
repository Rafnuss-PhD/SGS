function [gamma_x, gamma_y] = variogram_gridded_perso(X)


% X=(X-mean(X(:)))./std(X(:));

[ny,nx] = size(X);
gamma_x=nan(nx,1);
gamma_x_std=nan(nx,1);
gamma_y=nan(ny,1);
gamma_y_std=nan(ny,1);
% vertical
for j=1:ny
   gamma_y(j) = mean2((X(1:end-j+1,:)-X(j:end,:)).^2)/2;
   gamma_y_std(j) = std2((X(1:end-j+1,:)-X(j:end,:)).^2)/2;
end

for i=1:nx
   gamma_x(i) = mean2((X(:,1:end-i+1)-X(:,i:end)).^2)/2;
   gamma_x_std(i) = std2((X(:,1:end-i+1)-X(:,i:end)).^2)/2;
end

% figure; subplot(1,2,1); hold on;
% plot(gamma_y); plot(gamma_y+gamma_y_std,'--k'); plot(gamma_y-gamma_y_std,'--k')
%  subplot(1,2,2); hold on;
% plot(gamma_x); plot(gamma_x+gamma_x_std,'--k'); plot(gamma_x-gamma_x_std,'--k')






end