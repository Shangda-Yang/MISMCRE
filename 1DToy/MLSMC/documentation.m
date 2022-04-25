% ----------------------------------------------------------- %
% Results are obtained with the following settings:
%   * PDE: - \Delta^2 u(X) ) = 100X,   on (0,1)
%                       u(X) = 0,      on \partial [0,1]
%   * U[-1,1] prior
%   * Noise N(0,2^2) of statistical model
%   * Quantity of interest: x^2 x ~ \pi(dx) (posterior)
%   * 2 observations on [1/3; 2/3]
%   * Emperical values are computed with N = 10000 and 
%     M = 100 (realisations) and regression on L = 2:6
%   * MLSMC and SMC are implemented with M = 100
%   * SMC sampler tempering step \lambda = 0.5
%   * MCMC transition kernel RWM \beta = 0.5
% ----------------------------------------------------------- %