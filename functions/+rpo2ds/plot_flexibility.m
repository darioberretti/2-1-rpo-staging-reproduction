function figures = plot_flexibility(result, visible)
%PLOT_FLEXIBILITY Show the Sun-phase points and retuning choices.
%
% PROTOTYPE:
%   figures = rpo2ds.plot_flexibility(result, visible)
%
% DESCRIPTION:
%   Plot stored pre/post-retuning delta-v and post-retuning flight time.
%   The selected grid cases are marked on the untuned plots.
%
% INPUT:
%   result   Study result supplying venus/mars campaign structs and their tuning_log.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered names and figure handles, ready for export.
%              Plotting does not save result files or execute optimizations.

figures.flexibility_delta_v_before_tuning = local_plot(result, ...
    'total_delta_v_m_s', false, visible);
figures.flexibility_delta_v_after_tuning = local_plot(result, ...
    'total_delta_v_m_s', true, visible);
figures.flexibility_time_of_flight_after_tuning = local_plot(result, ...
    'time_of_flight_days', true, visible);
end

function handle = local_plot(result, value_field, tuned, visible)
%LOCAL_PLOT Draw Venus/Mars panels for one stored metric and tuning state.

handle = rpo.create_figure(visible,[1100,450]);
tiledlayout(handle, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
names = ["venus", "mars"];
labels = ["Venus", "Mars"];
for index = 1:2
    payload = result.(names(index));
    if tuned
        campaign = payload.tuned;
    else
        campaign = payload.untuned;
    end
    values = campaign.(value_field);
    ax = nexttile;
    plot(ax, campaign.phase_rad, values, 'ko', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 5, 'DisplayName', 'Computed transfer');
    hold(ax, 'on');
    if ~tuned
        % Several attempts at one phase still represent one transfer.
        selected = unique([payload.tuning_log.grid_index]);
        plot(ax, campaign.phase_rad(selected), values(selected), 'rd', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
            'DisplayName', 'Selected for retuning');
        title(ax, labels(index) + " before retuning");
        if index == 1
            legend(ax, 'Location', 'northeast');
        end
    else
        title(ax, labels(index));
    end
    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [0, 2 * pi]);
    xlabel(ax, '$\theta_{S,0}$ [rad]', 'Interpreter', 'latex');
    if strcmp(value_field, 'total_delta_v_m_s')
        ylabel(ax, '$\Delta V$ [m/s]', 'Interpreter', 'latex');
        if tuned
            ylim(ax, [0, 600]);
        else
            ylim(ax, [0, 1.05 * max(values)]);
        end
    else
        ylabel(ax, 'TOF [days]');
        ylim(ax, [10 + 10 * (index == 1), 100]);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 14);
end
end
