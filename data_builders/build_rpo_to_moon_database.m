clc; close all; clear;
addpath(genpath('functions'));

%% Purpose and instructions
% Build the phase-aware forward RPO-to-Moon encounter table for Runs 02 and 06.
% Run from the repository root or this script's folder. The frozen database
% seed and base MATLAB are sufficient; no downloaded data or earlier run is needed.
%
% The frozen recipe below uses 0-50 m/s impulses and an 85-day horizon.
% Parallel workers are optional, with serial fallback; progress is printed.
% The completed table and metadata go to data/generated/rpo_to_moon/.
% Select experiment_database_source = "generated" in the consuming runner.
% This builder does not create the legacy dense table or save figures.
%
% Compatible progress in cache/rpo_to_moon_checkpoint.mat resumes automatically.
% If progress is incompatible, explicitly discard that checkpoint to restart.
% Existing generated data requires confirmation, or
% overwrite_generated_database = true in batch execution.
% Reference files are never overwritten. See docs/DATA.md for units
% and encounter definitions, and README.md for shared workflow instructions.

%% Execution settings
use_parallel = true;
batch_size = 1000;
overwrite_generated_database = false;

%% Scientific settings
% Preserve the database's original normalization, including its time unit.
settings.mu = 0.0121506683;
settings.distance_unit_km = 384405;
settings.time_unit_days = 4.34811305;
settings.sun_rate_rad_tu = -0.925195985;
orbit_period_tu = (-2*pi/settings.sun_rate_rad_tu)/2;
maximum_impulse_m_s = 50;
maximum_impulse_vu = (maximum_impulse_m_s/1000/settings.distance_unit_km) ...
    * 3600 * 24 * settings.time_unit_days;
settings.impulses_vu = linspace(0, maximum_impulse_vu, 10);
settings.azimuths_rad = mod(linspace(0, 2*pi, 10), 2*pi);
settings.elevations_rad = linspace(-pi/2, pi/2, 10);
settings.phases_tu = linspace(0, orbit_period_tu, 50);
settings.direction = 1;
settings.horizon_days = 85;
settings.earth_guard_km = 6378;
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
    'producer_script', 'data_builders/build_rpo_to_moon_database.m');
database_file = fullfile(repository_root, 'data/generated/rpo_to_moon', ...
    'RPOMoonDatabase10x3_50phases.mat');
checkpoint_file = fullfile(repository_root, 'cache/rpo_to_moon_checkpoint.mat');

%% Build database
[RPOMoonDatabase, metadata] = rpomoon.build_database( ...
    settings, execution, database_file, checkpoint_file);
