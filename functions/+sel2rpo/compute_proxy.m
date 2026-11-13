function result = compute_proxy(cfg,repository_root,database_source)
%COMPUTE_PROXY Find the exact frontier of all duration-capped encounter pairs.
%
% PROTOTYPE:
%   result = sel2rpo.compute_proxy(cfg, repository_root, database_source)
%
% DESCRIPTION:
%   Pair SEL2-to-Moon and backward RPO-to-Moon database rows in
%   chunks. Minimize summed encounter duration and absolute difference
%   of lunar v-infinity speed magnitudes. This scalar mismatch is a
%   proxy, not a vector patch or optimized transfer cost. The frontier
%   uses all admissible pairs; random sampling affects only the cloud.
%
% INPUT:
%   cfg  Configuration with meta and sel2_to_rpo fields from
%        config/reproduction_config.m; physical units follow field suffixes.
%   cfg.sel2_to_rpo.proxy defines database filenames, duration cap [days],
%       chunk/sample sizes and RNG seed.
%   repository_root  Repository path containing data/{reference,generated}.
%   database_source  zenodo selects reference; generated selects generated.
%
% OUTPUT:
%   result  Result with data.pareto_front and plot_sample tables:
%             sel2_index/rpo_index are source row indices, tof_days is
%             the summed duration [days], and vinf_diff_m_s is the
%             speed-magnitude mismatch [m/s]. Data also records source
%             row counts and total/admissible pair counts; metrics
%             report frontier_rows and admissible_pairs.

% Resolve the selected encounter databases and retain their source paths.
s = cfg.sel2_to_rpo.proxy;
branch = 'reference'; if strcmp(database_source,'generated'), branch='generated'; end
relative_files = {fullfile('data',branch,'sel2_to_moon',s.sel2_database_file), ...
    fullfile('data',branch,'moon_to_rpo',s.rpo_database_file)};
for index=1:2
    if ~isfile(fullfile(repository_root,relative_files{index}))
        error('sel2rpo:MissingDatabase', ...
            ['Encounter database is missing:\n%s\nExtract the relevant database ZIP ', ...
            'or run the corresponding builder (see docs/DATA.md).'],relative_files{index});
    end
end
loaded = load(fullfile(repository_root,relative_files{1}),'result_explorer');
sel2 = loaded.result_explorer.SEL2MoonDatabase;
clear loaded;
loaded = load(fullfile(repository_root,relative_files{2}),'RPOMoonDatabase');
rpo_table = loaded.RPOMoonDatabase;
clear loaded;

% Reduce large source tables to encounter durations and speed magnitudes.
sel2_duration = double(sel2.t_enc_rel_days_em(:));
rpo_duration = abs(double(rpo_table.t_event_days(:)));
sel2_speed = vecnorm(double(sel2.vInfMoon_RTH_km_s),2,2);
rpo_speed = vecnorm(double(rpo_table.vInfMoon_RTH),2,2);
n_sel2 = height(sel2); n_rpo = height(rpo_table);
clear sel2 rpo_table;
rng(s.random_seed);
front = zeros(0,4); cloud = zeros(0,4); admissible = 0;
chunks = ceil(n_rpo/s.rpo_chunk_size);

% Accumulate the exact frontier while bounding the display-cloud size.
for chunk=1:chunks
    indices = (chunk-1)*s.rpo_chunk_size+1:min(chunk*s.rpo_chunk_size,n_rpo);
    duration = sel2_duration+rpo_duration(indices).';
    keep = duration<=s.maximum_duration_days;
    admissible = admissible+nnz(keep);
    if ~any(keep,'all'), continue; end
    metric = 1000*abs(sel2_speed-rpo_speed(indices).');
    [sel2_rows,rpo_rows] = find(keep);
    % Ordering is [duration, metric, SEL2 row, RPO row], including tie breaks.
    pairs = [duration(keep),metric(keep),sel2_rows,reshape(indices(rpo_rows),[],1)];
    front = reduce_front([front;reduce_front(pairs)]);
    count = min(size(pairs,1),max(1,ceil(0.10*s.plot_sample_size)));
    selection = randperm(size(pairs,1),count);
    cloud = [cloud;pairs(selection,:)]; %#ok<AGROW>
    if size(cloud,1)>s.plot_sample_size
        cloud = cloud(randperm(size(cloud,1),s.plot_sample_size),:);
    end
    fprintf('Database proxy %d/%d: %d admissible pairs, %d frontier rows\n', ...
        chunk,chunks,admissible,size(front,1));
end

% Restore source indices and units in the saved result tables.
names = {'sel2_index','rpo_index','tof_days','vinf_diff_m_s'};
data = struct('pareto_front',array2table(front(:,[3,4,1,2]),'VariableNames',names), ...
    'plot_sample',array2table(cloud(:,[3,4,1,2]),'VariableNames',names), ...
    'n_sel2',n_sel2,'n_rpo',n_rpo,'n_pairs_total',double(n_sel2)*n_rpo, ...
    'n_pairs_within_time_cap',admissible);
result = rpo.finalize_result(struct('data',data, ...
    'metrics',struct('frontier_rows',size(front,1),'admissible_pairs',admissible)), ...
    cfg,'run_04_sel2_to_rpo: database proxy',relative_files,'generated_result');
end

function front = reduce_front(points)
%REDUCE_FRONT Keep strict metric improvements in sorted [days,m/s,SEL2 row,RPO row] pairs.

points = sortrows(points,[1,2,3,4]);
% Strict improvement after sorting is the historical dominance rule.
keep = points(:,2)<[inf;cummin(points(1:end-1,2))];
front = points(keep,:);
end
