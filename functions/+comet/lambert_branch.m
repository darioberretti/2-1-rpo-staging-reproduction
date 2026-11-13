function [velocity,arrival_velocity,a,e,exitflag] = lambert_branch(r1,r2,tof_days,mu,revolutions,branch)
%LAMBERT_BRANCH Map one historical Lambert orientation/branch to the bundled solver.
%
% PROTOTYPE:
%   [velocity, arrival_velocity, a, e, exitflag] = comet.lambert_branch( ...
%       r1, r2, tof_days, mu, revolutions, branch)
%
% DESCRIPTION:
%   Use the z sign of cross(r1,r2) to sign the flight time, with
%   zero cross-product sign mapped to +1. Branch 0 uses negative
%   revolution count and branch 1 positive. Derive orbital elements
%   from the returned departure state only when the solver succeeds.
%
% INPUT:
%   r1, r2       [3 elements] Heliocentric inertial positions [km].
%   tof_days     Positive flight time [days].
%   mu           Solar gravitational parameter [km^3/s^2].
%   revolutions  Nonnegative integer revolution count.
%   branch       0 for the first branch, 1 for the second (multi-rev).
%
% OUTPUT:
%   velocity, arrival_velocity  [3 x 1] Inertial velocities [km/s].
%   a, e         Semimajor axis [km] and eccentricity [-], NaN on failure.
%   exitflag     FEX-Lambert status; 1 means success. Other values
%                  are returned to build_pair to stop branch enumeration.

persistent dependency_ready
if isempty(dependency_ready)
    if isempty(which('lambert'))
        root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        addpath(fullfile(root,'third_party','FEX-Lambert'));
    end
    dependency_ready = true;
end

% Translate the historical orientation and branch to signed FEX inputs.
cross_product = cross(r1,r2);
orientation = sign(cross_product(3));
if orientation==0, orientation=1; end
signed_revolutions = revolutions;
if branch==0, signed_revolutions=-revolutions; end
[velocity,arrival_velocity,~,exitflag] = lambert(r1(:).',r2(:).', ...
    orientation*tof_days,signed_revolutions,mu);
velocity = velocity(:); arrival_velocity = arrival_velocity(:);
a = NaN; e = NaN;
if exitflag~=1, return; end

% Recover conic elements from the successful inertial departure state.
radius = norm(r1);
energy = dot(velocity,velocity)/2-mu/radius;
a = -mu/(2*energy);
eccentricity = cross(velocity,cross(r1(:),velocity))/mu-r1(:)/radius;
e = norm(eccentricity);
end
