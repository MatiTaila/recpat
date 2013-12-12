function [y,V] = recpat_lda(data, labels, N)

sigma = 0;
for i=1:3
	sigma = sigma + ( mean(data((i-1)*N+1:i*N)) - mean(data) ) * ( mean(data((i-1)*N+1:i*N)) - mean(data) )' ;
end

sigma = sigma/3;

[U,S,V] = svds(data',size(data,1));
Y = inv(V)*(sigma);

data = Y*data;
descriptores = size(data,2);

for p=1:descriptores-2
	figure,
	plot3(ldaData(1:N,p),ldaData(1:N,p+1),ldaData(1:N,p+2),'*', 'color', colors{1});
	hold on
	grid on
	plot3(ldaData(N+1:2*N,p),ldaData(N+1:2*N,p+1),ldaData(N+1:2*N,p+2),'*', 'color', colors{5}),title('LDA'),xlabel('Primera dirección'),ylabel('Segunda dirección'),zlabel('Tercera dirección');
	plot3(ldaData(2*N+1:end,p),ldaData(2*N+1:end,p+1),ldaData(2*N+1:end,p+2),'*', 'color', colors{3})
	hold off
	pause
	close all
end