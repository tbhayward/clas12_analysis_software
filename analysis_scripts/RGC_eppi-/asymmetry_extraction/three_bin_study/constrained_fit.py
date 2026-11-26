"""Fixed-physics-ratio fits with first-order dilution response of free amplitudes."""
import numpy as np
from iminuit import Minuit
from uu_model_study.likelihood import ConditionalLikelihood, NAMES


def fit(groups, dilution, fixed_llcos=None):
    df, sigma = dilution
    if not (0 < df < 1 and 0 <= sigma < 1):
        raise ValueError('Invalid dilution input')
    cost = ConditionalLikelihood(groups, target=True)
    fixed = {0: 0., 1: 0.}
    if fixed_llcos is not None:
        if not np.isfinite(fixed_llcos):
            raise ValueError('Nonfinite fixed LLcos ratio')
        fixed[6] = df*fixed_llcos
    candidates = []
    for amplitude in [0., .03, -.03]:
        seed = np.array([0., 0., amplitude, amplitude, -amplitude, amplitude, 0.])
        for i, value in fixed.items():
            seed[i] = value
        m = Minuit(cost, *seed, name=NAMES, grad=cost.gradient)
        m.errordef = .5
        for i in fixed:
            m.fixed[i] = True
        m.migrad(ncall=4000)
        m.hesse()
        if not (m.valid and m.accurate):
            m.simplex(ncall=2000)
            m.strategy = 2
            m.migrad(ncall=6000)
            m.hesse()
        candidates.append(m)
    valid = [m for m in candidates if m.valid and m.accurate and m.covariance is not None]
    if not valid:
        return dict(status='invalid', events=cost.events,
                    candidates=[dict(nll=float(m.fval), valid=bool(m.valid), edm=float(m.fmin.edm)) for m in candidates])
    best = min(valid, key=lambda m: m.fval)
    beta = np.array(best.values)
    cb = np.asarray(best.covariance)
    free = np.array([i not in fixed for i in range(7)])
    values = beta.copy()
    values[3:] /= df
    jac = np.ones(7)
    jac[3:] /= df
    cov = cb*jac[:, None]*jac[None, :]
    response = np.zeros(7)
    if fixed_llcos is not None:
        # Differentiate the free score equations while beta_LLcos = Df * fixed_ratio.
        n = 1+cost.observed@beta
        d = 1+cost.average@beta
        hessian = (cost.observed/n[:, None]).T@(cost.observed/n[:, None])-(cost.average/d[:, None]).T@(cost.average/d[:, None])
        response[free] = -np.linalg.solve(hessian[np.ix_(free, free)], hessian[free, 6]*fixed_llcos)
        response[6] = fixed_llcos
    derivative = response*jac
    derivative[3:] -= values[3:]/df
    if fixed_llcos is not None:
        derivative[6] = 0.
    cov += sigma**2*np.outer(derivative, derivative)
    errors = np.sqrt(np.maximum(np.diag(cov), 0))
    norm = np.outer(errors, errors)
    correlation = np.divide(cov, norm, out=np.zeros_like(cov), where=norm > 0)
    eigen = np.linalg.eigvalsh(cb[np.ix_(free, free)])
    minimum = float(np.min(1+cost.allowed@beta))
    warnings = []
    if minimum < 1e-5:
        warnings.append('near_positivity_boundary')
    if eigen[-1]/eigen[0] > 1e8:
        warnings.append('weak_parameter_identification')
    if max(m.fval for m in valid)-best.fval > .1:
        warnings.append('multistart_minima_disagree')
    return dict(status='valid', events=cost.events, names=NAMES, values=values.tolist(), errors=errors.tolist(),
                covariance=cov.tolist(), correlation=correlation.tolist(), diluted_values=beta.tolist(),
                diluted_covariance=cb.tolist(), dilution=df, dilution_stat=sigma,
                dilution_derivative=derivative.tolist(), diluted_response=response.tolist(),
                fixed_parameters={NAMES[i]:float(values[i]) for i in fixed},
                free_parameters=[NAMES[i] for i in range(7) if free[i]],
                nll=float(best.fval), edm=float(best.fmin.edm), minimum_state_factor=minimum,
                covariance_condition=float(eigen[-1]/eigen[0]), multistart_nll=[float(m.fval) for m in candidates], warnings=warnings)
