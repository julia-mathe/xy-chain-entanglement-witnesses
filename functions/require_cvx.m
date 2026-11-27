function require_cvx()
% REQUIRE_CVX  Stop with a message if CVX is not installed.

if ~exist('cvx_begin', 'file')
    error('require_cvx:missing', ...
          ['This calculation needs CVX, which was not found.\n' ...
           'Install it from http://cvxr.com/cvx (unzip it, open its folder in MATLAB,\n' ...
           'type cvx_setup), then run this script again.']);
end
end
