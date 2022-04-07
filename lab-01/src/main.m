function main()
    pkg load statistics

    function myhist(X, bins, counts, R, m)
        n = length(X);
        delta = R / m;
        middles = zeros(1, m);
        xx = zeros(1, m);

        for i = 1:m
            xx(i) = counts(i) / (n * delta);
        endfor

        for i = 1:m
            middles(i) = bins(i + 1) - (delta / 2);
        endfor

        fprintf("    высоты столбцов гистограммы:\n");

        for i = 1:m
            fprintf("    [%d] : %f\n", i, xx(i));
        endfor

        fprintf("[проверка] площадь гистограммы s = %f\n", sum(xx) * delta);

        set(gca, "xtick", bins);
        set(gca, "ytick", xx);
        set(gca, "xlim", [min(bins) - 1, max(bins) + 1]);
        bar(middles, xx, 1, "facecolor", "g", "edgecolor", "w");

        X_n = m_min:(sigma / 100):m_max;
        X_pdf = normpdf(X_n, mu, sigma);
        plot(X_n, X_pdf, "r");
    end

    function mycdf(X, bins, counts)
        n = length(X);
        xx = zeros(1, m + 3);

        bins = [(min(bins) - 0.5) bins (max(bins) + 1)];
        counts = [0 counts 0];

        m = m + 2;

        acc = 0;

        for i = 2:m
            acc = acc + counts(i);
            xx(i) = acc / n;
        end

        xx(m + 1) = 1;

        X_n = (min(X) - 0.5):(sigma / 100):(max(X) + 1.5);
        X_cdf = normcdf(X_n, mu, s_2);
        plot(X_n, X_cdf, "r");


        for i = 2:m
            fprintf("x = %f : F(x) = %f\n", bins(i), xx(i));
        end


        set(gca, "xtick", bins);
        set(gca, "ylim", [0, 1.1]);
        set(gca, "ytick", xx);
        stairs(bins, xx);
    end

    X_raw = [-0.68, 0.71, 2.27, 0.38, 0.14, 0.06, 1.21, -0.59, 0.44, 1.98, 1.00, ...
             - 0.88, -0.08, 1.87, -0.74, 0.83, -1.45, 0.58, 0.48, 3.26, 0.02, 0.26, ...
             2.96, 1.78, 0.58, 0.08, -1.60, 1.26, 1.28, -0.36, 0.15, -0.38, -1.04, ...
             0.95, -2.17, -0.30, 1.09, 0.39, 1.06, 0.98, -2.55, 2.62, -1.58, 3.75, ...
             -1.43, 0.92, 2.75, -0.55, 1.48, -0.96, 0.50, 2.67, -0.58, 0.41, -0.46, ...
             -0.48, 1.68, -0.08, 1.76, 0.08, -1.15, 0.66, 1.54, 0.17, -0.20, 1.34, ...
             1.08, 1.59, -0.05, 0.15, -0.35, 0.58, -0.87, 1.73, -0.27, 0.00, -0.67, ...
             0.13, 1.75, -0.59, 1.31, 1.20, 0.53, 0.14, -0.35, 1.00, -0.01, 0.21, ...
             1.58, -0.02, 1.28, 1.34, -1.66, 0.30, 0.08, 0.66, -0.26, 1.54, 1.22, ...
             1.24, 0.11, 0.79, -0.83, 1.41, 0.17, 0.55, 1.60, 1.26, 1.06, 0.39, ...
             -0.77, 1.49, 0.92, -1.58, 1.19, 0.13, 0.26, -2.14, 0.08, -1.75];

    % X_raw = [13.7, 14.7, 13.4, 12.3, 13.0, 9.7, 16.2, 11.9, 13.8, 13.4, 13.7, 14.7, 13.7, 12.5, 15.5];

    X = sort(X_raw);

    % (a) минимальное и максимальное значение
    m_max = max(X);
    m_min = min(X);

    fprintf("(a) Максимальное значение выборки (M_max) = %f\n", m_max);
    fprintf("    Минимальное значение выборки  (M_min) = %f\n", m_min);
    fprintf("----------------------------------------\n");

    % (б) размах выборки
    r = m_max - m_min;
    fprintf("(б) Размах выборки (R) = %f\n", r);
    fprintf("----------------------------------------\n");

    % (в) вычисление оценок MX DX
    n = length(X);
    mu = sum(X) / n;
    s_2 = sum((X - mu).^2) / (n - 1);
    sigma = sqrt(s_2);

    fprintf("(в) Оценка математического ожидания (mu) = %f\n", mu);
    fprintf("    Оценка дисперсии (s_2) = %f\n", s_2);
    fprintf("----------------------------------------\n");

    % (г) группировка значений выборки в m = [log_2 n] + 2 интервала
    m = floor(log2(n)) + 2;

    bins = [];
    cur = m_min;

    for i = 1:(m + 1)
        bins(i) = cur;
        cur = cur + r / m;
    end

    eps = 1e-6;
    counts = [];
    j = 1;

    for i = 1:(m - 1)
        cur_count = 0;

        for j = 1:n

            if (bins(i) < X(j) || abs(bins(i) - X(j)) < eps) && X(j) < bins(i + 1)
                cur_count = cur_count + 1;
            endif

        endfor

        counts(i) = cur_count;
    endfor

    cur_count = 0;

    for j = 1:n

        if (bins(m) < X(j) || abs(bins(m) - X(j)) < eps) && (X(j) < bins(m + 1) || abs(bins(m + 1) - X(j)) < eps)
            cur_count = cur_count + 1;
        endif

    endfor

    counts(m) = cur_count;

    fprintf("(г) группировка значений выборки в m = [log_2 n] + 2 интервала:\n");

    for i = 1:(m - 1)
        fprintf("    [%f : %f) - %d вхожд.\n", bins(i), bins(i + 1), counts(i));
    end

    fprintf("    [%f : %f] - %d вхожд.\n", bins(m), bins(m + 1), counts(m));

    fprintf("----------------------------------------\n");

    fprintf("(д) построение гистограммы и графика функции плотности\n");
    fprintf("    распределения вероятностей нормальной случайной величины\n");

    figure;
    hold on;
    grid on;
    myhist(X, bins, counts, r, m);
    xlabel('X')
    ylabel('P')
    print -djpg hist.jpg
    hold off;

    fprintf("----------------------------------------\n");

    fprintf("(е) построение графика эмпирической функции распределения\n");
    fprintf("    и функции распределения нормальной случайной величины\n");

    figure;
    hold on;
    grid on;
    mycdf(X, bins, counts);
    xlabel('X')
    ylabel('F')
    print -djpg cdf.jpg
    hold off;
end

