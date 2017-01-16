trainX = load('input_1');
trainY = load('output_1_new');
% normalise the data to set it to mean = 0 and variance = 1
meany = mean(trainY);
vary = sqrt(var(trainY));
trainY = trainY - meany;
trainY = trainY / vary;

% set up the initial value for the hyper-parameters
% loghyp stores the logarithm value of the hyper-parameters
% loghyp.cov contains lengthscale_1 ... lengthscale_d and noise variance,
% totoally d + 1 value and d is the input dimension
% loghyp.lik is the signal variance
% loghyp.mean is set to be empty
loghyp.cov = zeros(1,size(trainX, 2)+1);
for i = 1:size(trainX,2)
    avg1 = log(2) - 2*log(sqrt(mean(abs(diff(trainX(:,i))))));
    avg2 = log(2) - 2*log(0.1*sqrt(mean(abs(diff(trainX(:,i))))));
    loghyp.cov(i) = max(avg1, avg2) + rand * 2;
end
loghyp.cov(end) = log(1);
loghyp.lik = log(2) - 2*log(sqrt(mean(abs(diff(trainY)))));
loghyp.mean = [];

% start training the hyper-paramters
% third parameter is -1000 where 1000 indicates the maxinum iteration time
loghyp = minimize(loghyp, @gp, -1000, @infExact, @meanZero, @covSEard, @likGauss, trainX, trainY);
% doing prediction on the test data, return the predict mean and variance
%[m s2] = gp(loghyp, @infExact, @meanZero, @covSEard, @likGauss, trainX, trainY, testX);
% since training data is normalised, in order to get the real mean and
% real variance, need to scale them back
%m = m * vary + meany;
%s2 = s2 * vary * vary;
%save('hyper.mat', 'loghyp', 'm', 's2');
