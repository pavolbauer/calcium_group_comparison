clear all
Cent1 = ones(1,9);
Cent2 = 1:1:9;
sig = 0.5;
nboot = 500;

tic
for a = 1:length(Cent1),
    xpdf(a,:) = sort([Cent1(a)+randn(1,20) Cent2(a)+randn(1,20)+rand(1,20).*10]); %allocate 200 points in each
    [dip(a), p(a)] = hartigansdipsigniftest(xpdf(a,:), nboot);
    
    subplot(3,3,a)
    hist(xpdf(a,:),-2:0.25:14)    
    title(['dip=',num2str(dip(a),3), ', p=',num2str(p(a),3)])
    xlim([-2 12])
end
% FixAxis([-2 12]);
toc