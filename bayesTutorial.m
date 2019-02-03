function bayesTutorial
close all
X = -10:0.1:10;
prior = normpdf(X,0,3);
like = normpdf(X,3,3);
post = (prior.*like);
post = post./trapz(X, post);
post2 = normpdf(X,1.5,sqrt(3^2/2));
figure; plot(X,prior); hold on; plot(X,like); plot(X,post);
plot(X,post2);



end