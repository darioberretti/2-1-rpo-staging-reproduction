function [rECI, rrMoonECI] = RTH2ECI_moon(rMoonECI, rRTH, ttRTH, type)
%RTH2ECI_MOON Apply the historical Moon-RTH position or velocity mapping to ECI.
%
% PROTOTYPE:
%   rECI = RTH2ECI_moon(rMoonECI, rRTH, [], type)
%   [rECI, rrMoonECI] = RTH2ECI_moon(rMoonECI, rRTH, ttRTH, type)
%
% DESCRIPTION:
%   Assume a circular equatorial Moon orbit, with radial axis
%   from Earth toward Moon, H along +Z_ECI, and T=H cross R. Position
%   mode rotates and translates by the Moon position. Velocity mode
%   preserves the literal legacy offset [0,sqrt(rEM/muE),0] before
%   rotation, with rEM=384000 and muE=398600; it is not recomputed from
%   a supplied ephemeris. All four arguments are required.
%
% INPUT:
%   rMoonECI  [3 elements] Initial Earth-to-Moon position [km].
%   rRTH      [N x 2 or 3] Moon-RTH position [km] for type=true, or
%               legacy velocity components [km/s] for type=false.
%               Missing H components are filled with zero.
%   ttRTH     [] for a fixed Moon position, or N relative epochs [s]
%               for the internal circular rotation n=sqrt(muE/rEM^3).
%   type      Logical/numeric mode: true position, false legacy velocity.
%
% OUTPUT:
%   rECI      [N x 3] ECI positions [km] or legacy velocities [km/s].
%   rrMoonECI [N x 3] Rotated Moon positions [km], assigned only with
%               nonempty ttRTH; do not request this output for ttRTH=[].

    rMoonECI = rMoonECI(:).';
    if numel(rMoonECI) ~= 3
        error('rMoonECI must be a 1x3 or 3x1 vector.');
    end

    % Handle 2D vs 3D RTH input
    if size(rRTH,2) == 2
        rRTH = [rRTH, zeros(size(rRTH,1),1)];  % append H=0 if missing
    elseif size(rRTH,2) ~= 3
        error('rRTH must be an Nx2 or Nx3 array.');
    end

    rEM = 384000;
    muE = 398600;
    n = sqrt(muE/rEM^3); % Rotation velocity of the Moon around Earth [rad/s]

    % Changing Moon position wrt time
    if ~isempty(ttRTH)
        rrMoonECI = zeros(numel(ttRTH), 3);
        rECI = zeros(numel(ttRTH), 3);
        for k = 1:numel(ttRTH)
            angle = ttRTH(k) * n;
            % Define the rotation matrix
            R = [cos(angle) -sin(angle) 0;
                 sin(angle)  cos(angle) 0;
                 0           0          1];
            % Assuming rMoonECI is defined and k is the index
            rrMoonECI(k, :) = (R * rMoonECI')';  % Rotate Moon position

            % Unit vectors of the RTH frame (planar, equatorial orbit)
            uR = rrMoonECI(k,:) / norm(rrMoonECI(k,:));  % radial: from Earth to Moon center
            uH = [0, 0, 1];                  % out-of-plane (fixed by assumption)
            uT = cross(uH, uR);              % prograde tangential
        
            % Build rotation matrix from RTH to ECI (columns are basis vectors)
            C = [uR(:), uT(:), uH(:)];       % 3x3
        
            if type
                % Rotate and translate all points
                rECI(k,:) = rRTH(k,:) * C.' + rrMoonECI(k,:);    % Nx3
            else
                % Rotate all points
                rECI(k,:) = (rRTH(k,:) + [0 sqrt(rEM/muE) 0]) * C.';
            end
        end
    else
        % Unit vectors of the RTH frame (planar, equatorial orbit)
        uR = rMoonECI / norm(rMoonECI);  % radial: from Earth to Moon center
        uH = [0, 0, 1];                  % out-of-plane (fixed by assumption)
        uT = cross(uH, uR);              % prograde tangential
    
        % Build rotation matrix from RTH to ECI (columns are basis vectors)
        C = [uR(:), uT(:), uH(:)];       % 3x3
    
        if type
            % Rotate and translate all points
            rECI = rRTH * C.' + rMoonECI;    % Nx3
        else
            % Rotate all points
            rECI = (rRTH + [0 sqrt(rEM/muE) 0])* C.';
        end
    end
end
