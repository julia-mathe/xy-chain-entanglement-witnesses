function map = viridis_colormap(n)
% VIRIDIS_COLORMAP  The viridis colour map.

if nargin < 1
    n = 256;
end
reference = [0.267004 0.004874 0.329415
             0.278826 0.175490 0.483397
             0.229739 0.322361 0.545706
             0.172719 0.448791 0.557885
             0.127568 0.566949 0.550556
             0.157851 0.683765 0.501686
             0.369214 0.788888 0.382914
             0.678489 0.863742 0.189503
             0.993248 0.906157 0.143936];
positions = linspace(0, 1, size(reference, 1));
map = interp1(positions, reference, linspace(0, 1, n));
end
