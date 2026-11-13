clc; close all; clear;
addpath(genpath('functions'));

%% Purpose and instructions
% Build the SEL2-to-Moon encounter table used by Run 04's database proxy.
% Run from the repository root or this script's folder. Base MATLAB and the
% tracked Run 04 orientation input are sufficient; parallel workers are optional.
%
% The halo is reconstructed internally. Each seed is propagated for up to
% 120 days in Sun-Earth CR3BP, then 500 days in Earth-Moon BCR4BP. Only the
% encounter table, configuration and compact diagnostics are saved; no plots
% or dense trajectory histories are retained.
%
% Output: data/generated/sel2_to_moon/. To use it, select
% experiment_database_source = "generated" in Run 04 and also generate its
% Moon-to-RPO table with build_moon_to_rpo_database.m.
% Compatible progress in cache/sel2_to_moon_checkpoint.mat resumes automatically.
% If progress is incompatible, explicitly discard that checkpoint to restart.
% Existing generated data requires confirmation, or
% overwrite_generated_database = true in batch execution. References are untouched.
% See docs/DATA.md for input locations and time conventions.

%% Execution settings
use_parallel = true;                    % Use workers when available; otherwise serial
batch_size = 1000;                      % Seeds per saved checkpoint
overwrite_generated_database = false;  % Authorize replacing an existing generated file

%% Scientific settings
settings.n_seed = 50000;                 % Uniform phases; period endpoint excluded
settings.theta_halo0_nd = 0;             % Initial halo phase [TU_SE]
settings.epsilon_vel_nd = 1e-6;           % Plus-branch velocity perturbation [VU_SE]
settings.theta_sun0_rad = 0;             % Sun phase offset [rad]
settings.t_first_days = 120;             % SE-CR3BP propagation horizon [days]
settings.t_second_days = 500;            % Subsequent BCR4BP horizon [days]
settings.n_samples_first = 600;
settings.n_samples_second = 5000;
settings.mu_em_nd = 0.0121506683;
settings.ws_sun_em_nd = -0.925195985;
settings.du_em_km = 384400;
settings.omega_em_rad_s = 2.665e-6;
settings.tu_em_days = (1/settings.omega_em_rad_s)/86400;
settings.alpha_plus_theta_const = pi;
settings.r_moon_soi_nd = (settings.mu_em_nd/(1-settings.mu_em_nd))^(2/5);
settings.earth_radius_km = 6378.1363;
settings.moon_radius_km = 1737.4;
settings.stop_earth_altitude_km = 500;
settings.stop_moon_altitude_km = 100;
settings.r_outer_guard_nd = 8;            % Radius in each model's own DU
settings.r_earth_guard_nd = (settings.earth_radius_km ...
    + settings.stop_earth_altitude_km)/settings.du_em_km;
settings.r_moon_guard_nd = (settings.moon_radius_km ...
    + settings.stop_moon_altitude_km)/settings.du_em_km;
settings.rel_tol = 1e-12;
settings.abs_tol = 1e-12;

%% Resolve inputs and output locations
repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repository_root,'functions')));
addpath(fullfile(repository_root,'config'));
cfg = reproduction_config();
settings.halo = cfg.sel2_to_rpo.halo;
settings.orientation_input = 'data/inputs/run_04_reference_inputs.mat';
loaded = load(fullfile(repository_root,settings.orientation_input),'inputs');
settings.reference_direction = loaded.inputs.manifold_reference_direction(:);
settings.orientation_sha256 = rpo2ds.sha256_file( ...
    fullfile(repository_root,settings.orientation_input));
execution = struct('use_parallel',use_parallel,'batch_size',batch_size, ...
    'overwrite',overwrite_generated_database);
database_file = fullfile(repository_root,'data/generated/sel2_to_moon', ...
    'SEL2MoonDatabase_BCR4BP_50000_1e-6_thetaSun0.mat');
checkpoint_file = fullfile(repository_root,'cache/sel2_to_moon_checkpoint.mat');

%% Build database
result_explorer = sel2moon.build_database(settings,execution,database_file,checkpoint_file);
