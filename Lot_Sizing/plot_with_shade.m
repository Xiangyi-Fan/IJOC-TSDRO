function p = plot_with_shade(x, y_avg, y_20,y_80, alpha, color, maker)
    y_avg = y_avg(~isnan(y_avg));
    y_20 = y_20(~isnan(y_20));
    y_80 = y_80(~isnan(y_80));
    length = size(y_avg,2);
    x = x(1:length);
    p = plot(x, y_avg, 'linewidth', 3,'color', color,'linestyle','--');
    p.Color(4) = 1;
    p.Marker = maker;
    p.MarkerEdgeColor = color;
    p.MarkerSize = 10;
    hold on
    x2 = [x, flip(x)];
    h = fill(x2,[y_20, flip(y_80)], color, 'LineStyle', 'none');
    set(h,'facealpha',alpha)
end