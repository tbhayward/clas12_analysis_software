"""Local Gaussian plus nonnegative quadratic background fits for provisional cuts."""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import xlogy

MASS_EDGES = np.linspace(.50, 1.50, 201)
MOM_EDGES = np.linspace(0, .50, 201)
ANGLE_EDGES = np.linspace(-.10, .10, 101)  # radians, matching the note
FIT_RANGES = {'mass': (.70, 1.30), 'momentum': (0, .25)}

def components(x, p, region):
    amp, mu, sigma, b0, b1, b2 = p
    u = (x-region[0])/(region[1]-region[0])
    return amp*np.exp(-.5*((x-mu)/sigma)**2), b0*(1-u)**2+2*b1*u*(1-u)+b2*u*u

def deviance_residual(obs, expected):
    expected = np.maximum(expected, 1e-12)
    return np.sign(expected-obs)*np.sqrt(np.maximum(0, 2*(expected-obs+xlogy(obs, obs/expected))))

def fit_histogram(hist, edges, kind):
    region = FIT_RANGES[kind]
    xall = (edges[:-1]+edges[1:])/2
    sel = (xall >= region[0]) & (xall < region[1])
    x, y = xall[sel], np.asarray(hist, float)[sel]
    result = dict(kind=kind, usable=False, warnings=[], fit_range=list(region), fit_entries=int(y.sum()))
    if y.sum() < 100:
        result['reason'] = 'fewer than 100 entries in fit range'
        return result
    upper_y = max(10., y.max()*10)
    if kind == 'mass':
        mu_bounds, sigma_bounds, seeds = (.78, 1.05), (.01, .20), [.86, .92, .98]
    else:
        mu_bounds, sigma_bounds, seeds = (.005, .18), (.005, .15), [.05, .08, .12]
    bounds = ([0, mu_bounds[0], sigma_bounds[0], 0, 0, 0],
              [upper_y, mu_bounds[1], sigma_bounds[1], upper_y, upper_y, upper_y])
    def residual(p):
        g, b = components(x, p, region)
        return deviance_residual(y, g+b)
    fits = []
    for seed in seeds:
        initial = [max(1., y.max()*.7), seed, .07 if kind == 'mass' else .035,
                   max(.1, y[0]), max(.1, np.percentile(y, 20)), max(.1, y[-1])]
        fits.append(least_squares(residual, initial, bounds=bounds, max_nfev=3000))
    opt = min(fits, key=lambda f: np.dot(f.fun, f.fun))
    pars = opt.x
    g, background = components(x, pars, region)
    expected = g+background
    ndf = len(x)-len(pars)
    pearson = float(np.sum((y-expected)**2/np.maximum(expected, 1e-12)))
    result.update(parameters=pars.tolist(), mu=float(pars[1]), sigma=float(pars[2]),
                  cut_low=float(max(0, pars[1]-2*pars[2]) if kind == 'momentum' else pars[1]-2*pars[2]),
                  cut_high=float(pars[1]+2*pars[2]), deviance=float(np.dot(opt.fun, opt.fun)),
                  chi2=pearson, ndf=ndf, chi2_ndf=pearson/ndf)
    # Local statistical errors conditional on the chosen Poisson fit model.
    errors=np.sqrt(np.maximum(0,np.diag(np.linalg.pinv(opt.jac.T@opt.jac))))
    result.update(mu_error=float(errors[1]),sigma_error=float(errors[2]))
    if not opt.success:
        result['reason'] = 'fit did not converge'
        return result
    for i in [1, 2]:
        relative_distance = min(pars[i]-bounds[0][i], bounds[1][i]-pars[i])/(bounds[1][i]-bounds[0][i])
        if relative_distance < .005:
            result['reason'] = 'peak mean or width reached a bound'
            return result
    # Require evidence for a peak over the same polynomial background.
    def bg_residual(b):
        return deviance_residual(y, components(x, [0, pars[1], pars[2], *b], region)[1])
    bgfit = least_squares(bg_residual, pars[3:], bounds=(0, upper_y), max_nfev=3000)
    improvement = float(np.dot(bgfit.fun, bgfit.fun)-result['deviance'])
    result['deviance_improvement'] = improvement
    if improvement < 9:
        result['reason'] = 'peak not sufficiently resolved above polynomial background'
        return result
    result['usable'] = True
    if result['chi2_ndf'] > 5:
        result['warnings'].append('poor Gaussian-plus-polynomial description')
    if result['cut_low'] < region[0] or result['cut_high'] > region[1]:
        result['warnings'].append('cut extends outside fitted range')
    result['reason'] = 'provisional cut, review warnings' if result['warnings'] else 'provisional cut'
    return result

def fit_angle(hist):
    """Gaussian visual guide only. Never supplies an exclusivity cut."""
    from scipy.optimize import curve_fit
    x=(ANGLE_EDGES[:-1]+ANGLE_EDGES[1:])/2
    select=np.abs(x)<.06
    x,y=x[select],hist[select].astype(float)
    if y.sum()<100:return dict(usable=False,reason='insufficient entries')
    def gaussian(x,a,mu,sigma):return a*np.exp(-.5*((x-mu)/sigma)**2)
    try:
        pars,_=curve_fit(gaussian,x,y,p0=[max(1,y.max()),0,.015],sigma=np.sqrt(np.maximum(1,y)),absolute_sigma=True,bounds=([0,-.04,.001],[max(10,y.max()*10),.04,.06]),maxfev=5000)
        return dict(usable=True,parameters=pars.tolist(),mu=float(pars[1]),sigma=float(pars[2]),role='visual guide only, no angular cut')
    except (RuntimeError,ValueError):return dict(usable=False,reason='Gaussian guide did not converge')
