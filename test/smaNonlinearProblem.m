function smaNonlinearProblem()

    yCp = [5.8377002519964755e+01, 2.9352296732047269e-03, 1.5061023667222263e-02, 1.3523701213590386e-01];
    kA = [0.0, 35.5, 1.59, 7.7];
    kD = [0.0, 1000.0, 1000.0, 1000.0];
    nu = [0.0, 4.7, 5.29, 3.7];
    sigma = [0.0, 11.83, 10.6, 10.0];
    lambda = 1.2e3;

    qStart = [1.0485785488181000e+03, 1.1604726694141368e+01, 1.1469542586742687e+01, 9.7852311988018670e+00];

    opts = optimset('Jacobian', 'on', 'MaxIter', 400, 'MaxFunEvals', 400, ...
        'TolFun', 1e-12, 'TolX', 0, 'Display', 'iter', ...
        'Algorithm', {'levenberg-marquardt',.001});
    [sol, res] = fsolve(@residual, qStart, opts);

    function [res, J] = residual(x)
        res = zeros(size(x));

        q0bar = x(1) - sum(sigma .* x);
        res(1) = x(1) - lambda + sum(nu .* x);

        c0powNu = (yCp(1)).^nu;
        q0barPowNu = q0bar.^nu;

        res(2:end) = kD(2:end) .* x(2:end) .* c0powNu(2:end) - kA(2:end) .* yCp(2:end) .* q0barPowNu(2:end);

        q0barPowNuM1 = q0bar.^(nu - 1.0);

        J = zeros(length(x));
        J(1,:) = [1.0, nu(2:end)];
        for i = 2:length(x)
            J(i, :) = -kA(i) * yCp(i) * nu(i) * q0barPowNuM1(i) .* [1, -sigma(2:end)];
            J(i,i) = J(i,i) + kD(i) * c0powNu(i);
        end
    end


end
