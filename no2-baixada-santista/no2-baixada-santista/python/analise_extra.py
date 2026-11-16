"""Análises adicionais: (1) validação cruzada do modelo de normalização meteorológica (anos sem pandemia);
(2) perfil de NO2 sobre a área de fundeio x faixa costeira; (3) teste de quebra (Pettitt) e salto nas mudanças de
versão do TROPOMI. Uso: python analise_extra.py pasta_GEE pasta_saida"""
import sys, json, numpy as np, pandas as pd
from PIL import Image
from scipy import stats
D, O = sys.argv[1].rstrip('/') + '/', sys.argv[2].rstrip('/') + '/'
x = pd.read_csv(O + 'no2_diario_met.csv', parse_dates=['data'])
Z = ['Sao_Paulo_centro', 'Cubatao_polo', 'Sao_Vicente', 'Santos_Porto', 'Praia_Grande', 'Guaruja', 'Mongagua', 'Itanhaem', 'Peruibe', 'Oceano_fundo']
R = {}
def X(d, anos):
    cols = [np.ones(len(d)), d.vel, d.vel ** 2, d.u10, d.v10, np.log1p(d.chuva_mm), d.t2m - 290]
    cols += [(d.mes == m).astype(float) for m in range(2, 13)] + [(d.dow == k).astype(float) for k in range(1, 7)]
    cols += [(d.ano == a).astype(float) for a in anos]
    return np.column_stack(cols)
CV = {}
for z in Z:
    g = x[x.zona == z].dropna(subset=['vel', 'no2']); g = g[g.no2 > 0].copy(); g['l'] = np.log(g.no2)
    err = {}
    for teste in (2019, 2021, 2022, 2023, 2024):
        anos_tr = [a for a in range(2019, 2026) if a not in (2020, teste)]
        tr = g[g.ano.isin(anos_tr)]; te = g[(g.ano == teste) & g.mes.isin([4, 5])]
        efe = anos_tr[1:]                                    # 1º ano é a referência
        b, *_ = np.linalg.lstsq(X(tr, efe), tr.l, rcond=None); s2 = (tr.l - X(tr, efe) @ b).var()
        Xt = X(te, efe)
        viz = [a for a in (teste - 1, teste + 1) if a in efe or a == anos_tr[0]]
        for j, a in enumerate(efe): Xt[:, -len(efe) + j] = (1.0 / len(viz)) if a in viz else 0.0
        pred = np.exp(Xt @ b + s2 / 2); err[teste] = (te.no2.mean() / pred.mean() - 1) * 100
    CV[z] = dict(**{str(k): v for k, v in err.items()}, dp=np.std(list(err.values()), ddof=1), media=np.mean(list(err.values())))
CV = pd.DataFrame(CV).T; CV.to_csv(O + 't_validacao_modelo.csv'); R['cv'] = CV.round(1).to_dict()

# (2) perfis a partir do mapa médio
im = Image.open(D + 'no2_media_2019_2025.tif'); A = np.array(im, float) * 1e6
sx, sy = im.tag_v2[33550][:2]; x0, y0 = im.tag_v2[33922][3:5]
def amostra(lon, lat):
    c = int((lon - x0) / sx); r = int((y0 - lat) / sy); return A[r - 1:r + 2, c - 1:c + 2].mean()
d_ = np.linspace(0, 45, 16)    # km para o sul a partir da costa
perf = {}
for nome, (lon0, lat0) in {'Barra de Santos (canal do porto)': (-46.30, -24.00), 'Praia Grande (costa sem porto)': (-46.50, -24.05),
                           'Itanhaém': (-46.78, -24.20)}.items():
    perf[nome] = [amostra(lon0 + 0.1 * k / 45 * 0, lat0 - k / 111.32) for k in d_]
P = pd.DataFrame(perf, index=d_); P.index.name = 'km_ao_sul'; P.to_csv(O + 'perfis_oceano.csv'); R['perfis'] = P.round(1).to_dict()
# área de fundeio: retângulo ao sul da barra (aprox. 24,05–24,15 S; 46,20–46,35 O) x faixa equivalente a oeste
def caixa(l0, l1, a0, a1):
    c0, c1 = int((l0 - x0) / sx), int((l1 - x0) / sx); r0, r1 = int((y0 - a1) / sy), int((y0 - a0) / sy); return A[r0:r1, c0:c1].mean()
R['fundeio'] = dict(area_fundeio=caixa(-46.35, -46.20, -24.15, -24.05), mesma_lat_oeste=caixa(-46.75, -46.60, -24.30, -24.20),
                    oceano_aberto=caixa(-46.40, -46.20, -24.42, -24.32))

# (3) quebra e salto nas mudanças de versão (anomalias mensais)
M = pd.read_csv(O + 'no2_mensal.csv', index_col=0, parse_dates=True)
an = M.groupby(M.index.month).transform(lambda s: s / s.mean() - 1) * 100
def pettitt(v):
    v = np.asarray(v); n = len(v); U = np.array([np.sign(v[:t + 1, None] - v[None, t + 1:]).sum() for t in range(n - 1)])
    k = np.abs(U).max(); return int(np.abs(U).argmax()), min(2 * np.exp(-6 * k ** 2 / (n ** 3 + n ** 2)), 1)
Q = {}
for z in Z:
    s = an[z].dropna(); t, p = pettitt(s.values)
    antes = s[s.index < '2021-07-01']; depois = s[s.index >= '2021-07-01']
    antes2 = antes[antes.index.year != 2020]
    Q[z] = dict(quebra=str(s.index[t + 1].date()) if p < 0.05 else '', p_pettitt=p, salto_jul2021=depois.mean() - antes2.mean(),
                p_salto=stats.mannwhitneyu(antes2, depois).pvalue)
Q = pd.DataFrame(Q).T; Q.to_csv(O + 't_quebras.csv'); R['quebras'] = Q.to_dict()
json.dump(R, open(O + 'resultados_extra.json', 'w'), indent=1, default=float)
print(CV.round(1).to_string()); print(P.round(1).to_string()); print(R['fundeio']); print(Q.round(3).to_string())
