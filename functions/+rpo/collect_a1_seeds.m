function seeds = collect_a1_seeds(x0_values, vy0_values, times, velocities, merge_tolerance)
%COLLECT_A1_SEEDS Find and de-duplicate velocity sign-change seeds.
%
% PROTOTYPE:
%   seeds = rpo.collect_a1_seeds(x0_values, vy0_values, times, velocities, ...
%       merge_tolerance)
%
% DESCRIPTION:
%   Locate zeros or sign changes in crossing vx along the sampled vy0 grid.
%   Interpolate seed velocity and time, then merge nearby velocities at each
%   x0 across crossing numbers. Inputs use the Earth-Moon rotating frame.
%
% INPUT:
%   x0_values        [Nx x 1] Initial x coordinates [DU].
%   vy0_values       [Nv x 1] Initial y velocities [VU].
%   times            [Nv x Nx x K] Crossing times [TU].
%   velocities       [Nv x Nx x K] Crossing vx values [VU].
%   merge_tolerance  Scalar seed-velocity merge tolerance [VU].
%
% OUTPUT:
%   seeds            Struct array with x0 [DU], vy0 [VU], t_half [TU], and
%                      crossing_number fields. K is the crossing count.

seeds = struct('x0', {}, 'vy0', {}, 't_half', {}, 'crossing_number', {});
zero_tolerance = 1e-6;

% Scan adjacent vy0 samples, keeping a shared accepted-velocity list for each x0.
for x_index = 1:numel(x0_values)
    accepted_vy = [];
    for crossing = 1:size(velocities, 3)
        vx = velocities(:, x_index, crossing);
        t = times(:, x_index, crossing);
        for vy_index = 1:numel(vy0_values)-1
            if ~all(isfinite([vx(vy_index:vy_index+1); t(vy_index:vy_index+1)]))
                continue
            end
            v1 = vx(vy_index);
            v2 = vx(vy_index + 1);
            s1 = local_sign(v1, zero_tolerance);
            s2 = local_sign(v2, zero_tolerance);

            % Keep a sampled zero or interpolate a bracketed sign change.
            if s1 == 0
                seed_vy = vy0_values(vy_index);
                seed_t = t(vy_index);
            elseif s1 * s2 < 0
                vy1 = vy0_values(vy_index);
                vy2 = vy0_values(vy_index + 1);
                seed_vy = vy1 - v1 * (vy2 - vy1) / (v2 - v1);
                seed_t = t(vy_index) + (seed_vy - vy1) ...
                    * (t(vy_index + 1) - t(vy_index)) / (vy2 - vy1);
            else
                continue
            end

            % Merge nearby seed velocities across all crossing numbers at this x0.
            if any(abs(seed_vy - accepted_vy) <= merge_tolerance)
                continue
            end
            seeds(end + 1) = struct( ...
                'x0', x0_values(x_index), ...
                'vy0', seed_vy, ...
                't_half', seed_t, ...
                'crossing_number', crossing); %#ok<AGROW>
            accepted_vy(end + 1) = seed_vy; %#ok<AGROW>
        end
    end
end
end

function output = local_sign(value, tolerance)
%LOCAL_SIGN Treat values within the absolute tolerance as zero.

if abs(value) <= tolerance
    output = 0;
else
    output = sign(value);
end
end
