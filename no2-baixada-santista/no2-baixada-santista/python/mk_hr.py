"""Mann-Kendall com correção de autocorrelação de Hamed e Rao (1998) e declive de Sen.
mk_hr(y, t) -> (sen_por_decada, p_original, p_hamed_rao)"""
import numpy as np
from scipy import stats

def mk_hr(y, t):
    y = np.asarray(y, float); t = np.asarray(t, float); ok = np.isfinite(y); y, t = y[ok], t[ok]; n = len(y)
    s = sum(np.sign(y[j] - y[i]) for i in range(n - 1) for j in range(i + 1, n))
    _, cont = np.unique(y, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(cont * (cont - 1) * (2 * cont + 5))) / 18
    z = lambda v: (s - np.sign(s)) / np.sqrt(v) if s != 0 else 0.0
    p0 = 2 * (1 - stats.norm.cdf(abs(z(var))))
    sen = stats.theilslopes(y, t)[0]
    r = stats.rankdata(y - sen * t); r = r - r.mean()
    soma = 0.0
    for k in range(1, n - 1):
        rho = np.sum(r[:-k] * r[k:]) / np.sum(r * r)
        if abs(rho) > stats.norm.ppf(0.975) / np.sqrt(n):          # só lags significativos
            soma += (n - k) * (n - k - 1) * (n - k - 2) * rho
    fator = max(1 + 2 / (n * (n - 1) * (n - 2)) * soma, 1e-6)
    p1 = 2 * (1 - stats.norm.cdf(abs(z(var * fator))))
    return sen * 10, p0, p1
