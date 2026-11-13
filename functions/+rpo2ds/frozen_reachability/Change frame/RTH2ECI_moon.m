function [rECI, rrMoonECI] = RTH2ECI_moon(rMoonECI, rRTH, ttRTH, type)
%RTH2ECI_MOON  Transform Moon-centered RTH positions into ECI coordinates
%
% USAGE
%   rECI = RTH2ECI_moon(rMoonECI, rRTH)
%
% INPUTS
%   rMoonECI : 1x3 or 3x1 vector. Moon position in ECI [same units as rRTH].
%              Assumes the Moon moves on a planar circular, equatorial orbit,
%              so the orbital angular momentum is along +Z_ECI.
%
%   rRTH     : Nx2 or Nx3 array. Positions expressed in the Moon RTH frame:
%              - Column 1: R  (radial, from Moon toward Earth opposite? see note)
%              - Column 2: T  (tangential, prograde direction of Moon motion)
%              - Column 3: H  (out-of-plane, +Z of the orbital plane) [optional]
%              If Nx2 is provided, H is assumed zero.
%
% OUTPUT
%   rECI     : Nx3 array. The same points expressed in ECI.
%
% FRAME DEFINITIONS (planar circular orbit assumption)
%   - uR = rMoonECI / |rMoonECI| (points from Earth toward Moon's center)
%   - uH = [0 0 1]  (constant, because the Moon's orbit is taken as equatorial)
%   - uT = uH x uR  (prograde direction of the Moon's motion)
%
% NOTES
%   • This function just rotates and translates: r_ECI = r_Moon + R*[R T H]^T.
%   • If your ECI z-axis isn't aligned with the Moon's orbital angular momentum,
%     this planar/equatorial assumption is not valid; you would then need the
%     instantaneous Moon velocity (or orbital elements) to build uT/uH.
%
% EXAMPLE
%   rMoonECI = [384400, 0, 0];        % km
%   rRTH     = [1000, 0; 0, 1000];    % two points: +R 1000 km, +T 1000 km
%   rECI     = RTH2ECI_moon(rMoonECI, rRTH);
%

    % Ensure row vector for the Moon position
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
