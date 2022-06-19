function [Sp,tp,A,b,tottime] = get_ellipsoids(S_supp,t_supp,points)
    % S_supp*x <= t_supp is the support set
    % points: centres of voronoi diagram; each row is a point
    % returns:
    %    Sp,tp: {x: Sp{j}*x <= tp{j} is boundary of j'th region
    %    A, b; j'th ellipsoid is {x : || A{j}*x + b{j} <= 1 ||}

    [M,dim] = size(points);
    N = size(S_supp,2);
    X = points;
    nm = sum(X.*X,2);

    A = cell(M,1);
    b = cell(M,1);
    Sp = cell(M,1);
    tp = cell(M,1);
    tottime = 0;
    for j = 1:M
        S1 = 2*(X-X(j,:));
        t1 = nm - nm(j);
        Sp{j} = [S_supp;[S1,zeros(M,N-dim)]];
        tp{j} = [t_supp;t1];
        [A{j},b{j},~,out] = polyhedron_copositive(Sp{j},tp{j},dim);
        tottime = tottime + out.solvertime;
    end

end