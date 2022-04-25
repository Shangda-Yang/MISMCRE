% ----------------------------------------------------------------------- %
% Results are obtained with the following settings:
%   * PDE: - \Delta \cdot a(X) \Delta u(X) ) = 100, on (0,1)^2
%                                       u(X) = 0,   on \partial [0,1]^2
%   * a(X)(Z) = 3 + x_1cos(3z_1)sin(3z_2) + x_2cos(z_1)sin(z_2)
%   * U[-1,1]^2 prior
%   * Noise N(0,0.5^2) of statistical model
%   * Quantity of interest: x_1^2 + x_2^2, x_1, x_2 ~ \pi(dx) (posterior)
%   * 4 observations on [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75]
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% run data_generator --> generator.m first for data and reference solution
% first time --> run emperical test before run complesity test
% For MISMC:
%   MISMC_2D_nlin_fixedtemp --> emperical_test.m for emperical mean and var
%                           --> complexity_test.m for complexity test
% For SLSMC:
%   MISMC_2D_nlin_fixedtemp --> smc_test.m for complexity test
% For MLSMC:
%   MLSMC_2D_nlin --> mlsmc_emperical.m for emperical mean and var
%                 --> complexity_test.m for complexity test
% ----------------------------------------------------------------------- %