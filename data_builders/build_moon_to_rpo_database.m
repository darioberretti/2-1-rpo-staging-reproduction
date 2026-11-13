clc; close all; clear;
addpath(genpath('functions'));

%% Purpose and instructions
% Build the Moon-to-RPO encounter table used by Run 04's database proxy.
% Run from the repository root or this script's folder. The frozen database seed
% and base MATLAB are sufficient; parallel workers are optional.
%
% The recipe below propagates backward from 50 RPO phases with 0-5 m/s
% impulses for up to 200 days. Stored event times remain negative;
% Run 04 converts them to positive leg durations.
% The completed table and metadata go to data/generated/moon_to_rpo/.
% Select experiment_database_source = "generated" in Run 04; its SEL2-Moon
% database must also be present in data/generated/sel2_to_moon/.
% No figures or full trajectory histories are saved.
%
% Compatible progress in cache/moon_to_rpo_checkpoint.mat resumes automatically.
% If progress is incompatible, explicitly discard that checkpoint to restart.
% Existing generated data requires confirmation, or
% overwrite_generated_database = true in batch execution.
% Reference files are never overwritten. See docs/DATA.md for data locations
% and encounter definitions.

%% Execution settings
use_parallel = true;
batch_size = 1000;
overwrite_generated_database = false;

%% Scientific settings
% Preserve this database's normalization and inclusive sampling endpoints.
settings.mu = 0.0121506683;
settings.distance_unit_km = 384405;
settings.time_unit_days = 4.34811305;
settings.sun_rate_rad_tu = -0.925195985;
orbit_period_tu = (-2*pi/settings.sun_rate_rad_tu)/2;
maximum_impulse_m_s = 5;
maximum_impulse_vu = (maximum_impulse_m_s/1000/settings.distance_unit_km) ...
    * 3600 * 24 * settings.time_unit_days;
settings.impulses_vu = linspace(0, maximum_impulse_vu, 10);
settings.azimuths_rad = mod(linspace(0, 2*pi, 10), 2*pi);
settings.elevations_rad = linspace(-pi/2, pi/2, 10);
settings.phases_tu = linspace(0, orbit_period_tu, 50);
settings.direction = -1;
settings.horizon_days = 200;
settings.earth_guard_km = 6878; % Earth radius plus a 500-km margin
settings.outer_guard_km = 2000000;
settings.moon_radius_km = 1737.4;
settings.moon_soi_km = 384405 * (7.34767309e22/5.97219e24)^(2/5);
settings.relative_tolerance = 1e-10;
settings.absolute_tolerance = 1e-10;
settings.retry_factor = 50;

%% Resolve inputs and output locations
repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repository_root, 'functions')));
settings.seed_file = 'data/inputs/rpo_moon_database_seed.mat';
loaded = load(fullfile(repository_root, settings.seed_file), 'rpo_moon_database_seed');
settings.seed_state_2d = loaded.rpo_moon_database_seed.state_2d(:);
settings.seed_sha256 = rpo2ds.sha256_file(fullfile(repository_root, settings.seed_file));
execution = struct('use_parallel', use_parallel, 'batch_size', batch_size, ...
    'overwrite', overwrite_generated_database, ...
    'producer_script', 'data_builders/build_moon_to_rpo_database.m');
database_file = fullfile(repository_root, 'data/generated/moon_to_rpo', ...
    'RPOMoonDatabase10x3_50phases_backward_5dv.mat');
checkpoint_file = fullfile(repository_root, 'cache/moon_to_rpo_checkpoint.mat');

%% Build database
[RPOMoonDatabase, metadata] = rpomoon.build_database( ...
    settings, execution, database_file, checkpoint_file);
