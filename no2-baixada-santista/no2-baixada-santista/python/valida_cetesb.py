"""Validação do NO2 do TROPOMI com as estações da CETESB (QUALAR, médias horárias).
Uso: python valida_cetesb.py pasta_csv_qualar pasta_pol"""
import sys, glob, json, numpy as np, pandas as pd
from scipy import stats
I, O = sys.argv[1].rstrip('/') + '/', sys.argv[2].rstrip('/') + '/'
reg = []
for f in glob.glob(I + '*.csv'):
    L = open(f, encoding='latin1').read().split('\n')
    est = L[4].split(';')[1].strip()
    d = pd.read_csv(f, sep=';', skiprows=8, header=None, names=['data', 'hora', 'v'], encoding='latin1', decimal=',')
    d['est'] = est; reg.append(d)
h = pd.concat(reg).dropna(subset=['v']).drop_duplicates(['est', 'data', 'hora'])
h['v'] = pd.to_numeric(h.v.astype(str).str.replace(',', '.'), errors='coerce'); h = h.dropna(subset=['v'])
h['hr'] = h.hora.str[:2].astype(int)
h['dia'] = pd.to_datetime(h.data, format='%d/%m/%Y')
h.loc[h.hr == 24, 'dia'] += pd.Timedelta(days=1); h.loc[h.hr == 24, 'hr'] = 0   # 24:00 = 0h do dia seguinte
# médias horárias marcadas no fim da hora: 13:00 e 14:00 = 12h–14h (passagem ~13h30)
p = h[h.hr.isin([13, 14])].groupby(['est', 'dia']).v.agg(['mean', 'count']); p = p[p['count'] == 2]['mean'].rename('st13')
dm = h.groupby(['est', 'dia']).v.agg(['mean', 'count']); dm = dm[dm['count'] >= 18]['mean'].rename('st24')
S = pd.concat([p, dm], axis=1).reset_index()
x = pd.read_csv(O + 'no2_diario_met.csv', parse_dates=['data'])
mapa = {'Santos-Ponta da Praia': 'Santos_Porto', 'Cubatão-Centro': 'Cubatao_polo'}
R = {}; pares = {}
for est, z in mapa.items():
    s = S[S.est == est].set_index('dia'); t = x[x.zona == z].set_index('data').no2
    m = s[['st13', 'st24']].join(t.rename('trop'), how='inner').dropna(subset=['st13', 'trop'])
    pares[est] = m
    r_d, _ = stats.pearsonr(m.trop, m.st13); rs_d = stats.spearmanr(m.trop, m.st13)[0]
    mm = m.resample('MS').mean().dropna(); mm = mm[m.resample('MS').size().reindex(mm.index) >= 5]
    r_m = stats.pearsonr(mm.trop, mm.st13)[0]
    an = mm.groupby(mm.index.month).transform(lambda v: v / v.mean())   # anomalias (sem sazonalidade)
    r_an = stats.pearsonr(an.trop, an.st13)[0]
    sv = s.st24.dropna(); ss = sv.groupby(sv.index.month).mean()
    st = s[['st13', 'st24']].copy()
    def dom(col):
        v = st[col].dropna(); rel = v / v.rolling(28, center=True, min_periods=14).mean(); g = rel.groupby(v.index.dayofweek).mean()
        return g[6] / g[:5].mean()
    def pand(col):
        v = s[col].dropna(); am = v[v.index.month.isin([4, 5])]
        return (am[am.index.year == 2020].mean() / am[am.index.year.isin([2019, 2021, 2022])].mean() - 1) * 100
    def pand_t():
        am = t[t.index.month.isin([4, 5])]; return (am[am.index.year == 2020].mean() / am[am.index.year.isin([2019, 2021, 2022])].mean() - 1) * 100
    anual = sv.groupby(sv.index.year).agg(['mean', 'count'])
    anual = anual[anual['count'] >= 0.75 * 365]['mean']
    R[est] = dict(n_dias=len(m), r_diario=r_d, rs_diario=rs_d, n_meses=len(mm), r_mensal=r_m, r_anomalias=r_an,
                  inv_ver_st=ss[[6, 7, 8]].mean() / ss[[12, 1, 2]].mean(),
                  inv_ver_trop=t[t.index.month.isin([6, 7, 8])].mean() / t[t.index.month.isin([12, 1, 2])].mean(),
                  dom_st13=dom('st13'), dom_st24=dom('st24'),
                  pand_st13=pand('st13'), pand_st24=pand('st24'), pand_trop=pand_t(),
                  media_st24=sv.mean(), media_anual=anual.round(1).to_dict(),
                  dias_acima_OMS25=float((sv > 25).mean() * 100), max_horaria=h[h.est == est].v.max(),
                  horas_acima_200=int((h[h.est == est].v > 200).sum()))
    pd.DataFrame({'trop': mm.trop, 'st13': mm.st13}).to_csv(O + f'cetesb_mensal_{z}.csv')
json.dump(R, open(O + 'resultados_cetesb.json', 'w'), indent=1, default=float, ensure_ascii=False)
for k, v in R.items(): print(k); [print('  ', a, round(b, 3) if isinstance(b, float) else b) for a, b in v.items()]
pd.to_pickle(pares, O + 'cetesb_pares.pkl')
