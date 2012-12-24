%%
CF = [1 1.25 1.5 2 2.5 3 4];
m = [3.25 2.95 2.75 2.45 2.2 2.05 1.8];
b = m;
x = 0.5:0.01:2;
colors = 'rbgmckrbgmck';
figure(1), hold on;
for i=1:length(CF)
    if i>1
        prevFreqs = CF(i-1)*x;
        prevOffsets = -m(i-1)*x+b(i-1);
        previousIntercept = interp1(prevFreqs,offsets(:,i-1),CF(i));
    else
        previousIntercept = 0;
    end
    
    offsets(:,i) = -m(i)*x+b(i)+previousIntercept;
    plot(CF(i)*x,offsets(:,i),colors(i));
end
hold off;

hold on; 
for i=1:length(CF)
    plot(repmat(CF(i),2,1),[-7 2]',[colors(i) ':']); 
end
hold off;

