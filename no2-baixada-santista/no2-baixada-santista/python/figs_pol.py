"""Figuras do artigo de poluição do ar na Baixada Santista. Uso: python figs_pol.py pasta_GEE pasta_geo pasta_saida"""
import sys, json, unicodedata, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle, Rectangle
from matplotlib.colors import LinearSegmentedColormap, LightSource, TwoSlopeNorm
from matplotlib.ticker import FuncFormatter
from PIL import Image
D, G, O = [a.rstrip('/') + '/' for a in sys.argv[1:4]]
CM = 1 / 2.54
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Carlito', 'Liberation Sans'], 'font.size': 8.5,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400, 'axes.linewidth': 0.6})
V = FuncFormatter(lambda v, p: ('%g' % round(v, 3)).replace('-', '−').replace('.', ','))
norm = lambda x: unicodedata.normalize('NFKD', str(x)).encode('ascii', 'ignore').decode().lower().strip()
Z = ['Sao_Paulo_centro', 'Cubatao_polo', 'Sao_Vicente', 'Santos_Porto', 'Praia_Grande', 'Guaruja', 'Mongagua', 'Itanhaem', 'Peruibe', 'Oceano_fundo']
NOME = {'Sao_Paulo_centro': 'São Paulo (centro)', 'Cubatao_polo': 'Cubatão (polo)', 'Sao_Vicente': 'São Vicente',
        'Santos_Porto': 'Santos/Porto', 'Praia_Grande': 'Praia Grande', 'Guaruja': 'Guarujá', 'Mongagua': 'Mongaguá',
        'Itanhaem': 'Itanhaém', 'Peruibe': 'Peruíbe', 'Oceano_fundo': 'Oceano (fundo)'}
SIG = {'Sao_Paulo_centro': 'SP', 'Cubatao_polo': 'CB', 'Sao_Vicente': 'SV', 'Santos_Porto': 'ST', 'Praia_Grande': 'PG',
       'Guaruja': 'GJ', 'Mongagua': 'MG', 'Itanhaem': 'IT', 'Peruibe': 'PE', 'Oceano_fundo': 'OC'}
XY = {'Santos_Porto': (-46.315, -23.950), 'Cubatao_polo': (-46.395, -23.865), 'Sao_Vicente': (-46.390, -23.960),
      'Guaruja': (-46.250, -23.990), 'Praia_Grande': (-46.420, -24.010), 'Mongagua': (-46.620, -24.090),
      'Itanhaem': (-46.790, -24.180), 'Peruibe': (-47.000, -24.310), 'Oceano_fundo': (-46.300, -24.350), 'Sao_Paulo_centro': (-46.633, -23.550)}
COR = {'Sao_Paulo_centro': '#444444', 'Cubatao_polo': '#c0392b', 'Santos_Porto': '#e08e0b', 'Sao_Vicente': '#8e44ad',
       'Praia_Grande': '#2a6ebb', 'Guaruja': '#16a085', 'Peruibe': '#27ae60', 'Oceano_fundo': '#7fb3d5', 'Mongagua': '#95a5a6', 'Itanhaem': '#bdc3c7'}
ASP = 1 / np.cos(np.radians(23.9))

def aneis(g):
    t = g['type']
    if t == 'Polygon': yield g['coordinates']
    elif t == 'MultiPolygon': yield from g['coordinates']
    elif t == 'GeometryCollection':
        for gg in g['geometries']: yield from aneis(gg)
def desenha(ax, feats, **kw):
    for f in feats:
        vs, cs = [], []
        for poly in aneis(f['geometry']):
            for r in poly:
                r = np.asarray(r)[:, :2]
                if len(r) < 3: continue
                vs += r.tolist(); cs += [Path.MOVETO] + [Path.LINETO] * (len(r) - 2) + [Path.CLOSEPOLY]
        if vs: ax.add_patch(PathPatch(Path(vs, cs), **kw))
def tif(f):
    im = Image.open(f); A = np.array(im, dtype=float); sx, sy = im.tag_v2[33550][:2]; x0, y0 = im.tag_v2[33922][3:5]
    return A, (x0, x0 + A.shape[1] * sx, y0 - A.shape[0] * sy, y0)
def grade(ax, xs, ys, fs=6.8):
    ax.set_xticks(xs); ax.set_yticks(ys)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{abs(v):.1f}°O".replace('.', ',')))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{abs(v):.1f}°S".replace('.', ',')))
    ax.tick_params(labelsize=fs, length=2, direction='in', top=True, right=True); ax.grid(color='#666', lw=0.25, ls=':', alpha=0.5)
    for s in ('top', 'right'): ax.spines[s].set_visible(True)
def escala(ax, x, y, km=20):
    d = km / (111.32 * np.cos(np.radians(24)))
    for i in range(2): ax.add_patch(Rectangle((x + i * d / 2, y), d / 2, 0.012, fc='k' if i == 0 else 'white', ec='k', lw=0.5, zorder=9))
    for v_, t in ((0, '0'), (d / 2, str(km // 2)), (d, f'{km} km')): ax.text(x + v_, y + 0.022, t, ha='center', fontsize=6.3, zorder=9)
def norte(ax, x, y):
    ax.annotate('', xy=(x, y + 0.07), xytext=(x, y), arrowprops=dict(arrowstyle='-|>,head_width=0.35,head_length=0.7', lw=1.1, color='k'), zorder=9)
    ax.text(x, y + 0.085, 'N', ha='center', fontsize=8, fontweight='bold', zorder=9)

mun = json.load(open(G + 'municipios_area.geojson'))['features']
msp = json.load(open(G + 'municipios_sp.geojson'))['features']
pai = json.load(open(G + 'paises.geojson'))['features']
alvo = [f for f in mun if norm(f['properties']['ADM2_NAME']) in ('santos', 'sao vicente', 'guaruja', 'praia grande', 'mongagua', 'itanhaem', 'peruibe', 'cubatao')]
im = Image.open(G + 'srtm.tif'); Zs = np.array(im, dtype=float); sx, sy = im.tag_v2[33550][:2]; x0, y0 = im.tag_v2[33922][3:5]
zext = (x0, x0 + Zs.shape[1] * sx, y0 - Zs.shape[0] * sy, y0); mar = Zs <= 0
mar[:3] = False                                   # primeiras linhas do SRTM sem dado (borda), não são oceano

EXT = (-47.1, -45.9, -24.45, -23.4); EXT1 = (-47.1, -45.9, -24.45, -23.5)

def base_mapa(ax, cor_lim='#555'):
    desenha(ax, mun, fc='none', ec=cor_lim, lw=0.3, zorder=3)
    for z, (x, y) in XY.items():
        ax.add_patch(Circle((x, y), 4 / 111.32 / np.cos(np.radians(24)), fc='none', ec='k', lw=0.7, zorder=5))
        ax.text(x, y, SIG[z], ha='center', va='center', fontsize=6, fontweight='bold', zorder=6,
                bbox=dict(boxstyle='round,pad=0.05', fc='white', ec='none', alpha=0.65))
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_aspect(ASP)

# ======== Fig 1: localização
hips = LinearSegmentedColormap.from_list('h', ['#7fb27a', '#b9d39a', '#e8e1a9', '#d6b583', '#b08b67', '#8a6d5a', '#f2efe9'])
rgb = LightSource(315, 40).shade(np.clip(Zs, 0, 1300), cmap=hips, vmin=0, vmax=1300, blend_mode='soft', vert_exag=0.08, dx=90, dy=90)
rgb[mar] = matplotlib.colors.to_rgba('#d6e8f4')
fig = plt.figure(figsize=(16 * CM, 11.2 * CM)); ax = fig.add_axes([0.06, 0.07, 0.68, 0.91])
ax.set_facecolor('#d6e8f4'); ax.imshow(rgb, extent=zext, zorder=0, interpolation='bilinear')
desenha(ax, alvo, fc='none', ec='k', lw=0.9, zorder=4); base_mapa(ax)
for t, x, y in (('Porto de\nSantos', -46.27, -23.905), ('Polo industrial\nde Cubatão', -46.43, -23.81), ('OCEANO ATLÂNTICO', -46.75, -24.40),
                ('SERRA DO MAR', -46.82, -23.83), ('Região Metropolitana\nde São Paulo', -46.40, -23.58)):
    ax.text(x, y, t, ha='center', fontsize=6.5, style='italic', color='#222' if 'OCEANO' not in t else '#4a7aa5', zorder=7,
            rotation=-25 if 'SERRA' in t else 0, fontweight='bold' if 'SERRA' in t else 'normal',
            bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.6) if 'OCEANO' not in t and 'SERRA' not in t else None)
ax.set_ylim(EXT1[2], EXT1[3]); grade(ax, [-47.0, -46.6, -46.2], [-24.4, -24.0, -23.6]); escala(ax, -47.05, -24.43); norte(ax, -45.97, -23.66)
cax = ax.inset_axes([0.40, 0.10, 0.22, 0.025])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=hips, norm=plt.Normalize(0, 1300)), cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=5.8, pad=1); cb.set_label('Altitude (m)', fontsize=6, labelpad=1); cax.xaxis.set_label_position('top')
a1 = fig.add_axes([0.76, 0.52, 0.23, 0.46])
for f in pai:
    vs = f['properties']['country_na'] == 'Brazil'; desenha(a1, [f], fc='#e9e3d4' if vs else '#f4f4f4', ec='#777', lw=0.3)
desenha(a1, msp, fc='#c0392b', ec='none'); a1.set_xlim(-82, -34); a1.set_ylim(-56, 13); a1.set_aspect(1 / np.cos(np.radians(20)))
a1.set_facecolor('#d6e8f4'); a1.set_xticks([]); a1.set_yticks([])
a2 = fig.add_axes([0.76, 0.07, 0.23, 0.42]); desenha(a2, msp, fc='#efe9da', ec='#bbb', lw=0.08)
a2.add_patch(Rectangle((EXT[0], EXT[2]), EXT[1] - EXT[0], EXT[3] - EXT[2], fc='#c0392b', ec='k', lw=0.6, alpha=0.6))
a2.set_xlim(-53.3, -44.0); a2.set_ylim(-25.5, -19.6); a2.set_aspect(1 / np.cos(np.radians(22.5))); a2.set_facecolor('#d6e8f4')
a2.set_xticks([]); a2.set_yticks([]); a2.text(-49, -22.2, 'SÃO PAULO', fontsize=6, ha='center')
for a_, t in ((ax, '(a)'), (a1, '(b)'), (a2, '(c)')): a_.text(0.03, 0.97, t, transform=a_.transAxes, fontsize=7.5, fontweight='bold', va='top')
fig.savefig(O + 'fig1_localizacao.png'); plt.close()

# ======== Fig 2: mapa médio de NO2 + barras por zona
A, ext = tif(D + 'no2_media_2019_2025.tif'); A = A * 1e6
cmap = LinearSegmentedColormap.from_list('no2', ['#fff7ec', '#fdd49e', '#fc8d59', '#d7301f', '#7f0000'])
fig = plt.figure(figsize=(16 * CM, 7.6 * CM)); ax = fig.add_axes([0.04, 0.1, 0.46, 0.86])
h = ax.imshow(A, extent=ext, cmap=cmap, vmin=15, vmax=180, interpolation='bilinear', zorder=0)
ax.contour(np.where(mar, 1, 0), levels=[0.5], extent=zext, origin='upper', colors='#1f4e79', linewidths=0.6, zorder=2)
base_mapa(ax); grade(ax, [-47.0, -46.5, -46.0], [-24.4, -24.0, -23.6], 6.3)
cb = fig.colorbar(h, ax=ax, fraction=0.04, pad=0.02); cb.set_label('NO$_2$ (µmol m$^{-2}$)', fontsize=7); cb.ax.tick_params(labelsize=6.5)
ax.text(0.03, 0.97, '(a)', transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')
t = pd.read_csv(O + 't_medias.csv', index_col=0).loc[Z]
b = fig.add_axes([0.75, 0.15, 0.23, 0.81])
yy = np.arange(len(Z))[::-1]
b.barh(yy, t.media, color=[COR[z] for z in Z], height=0.7)
b.errorbar(t.mediana, yy, fmt='|', color='k', ms=7, mew=1.2, label='mediana')
for yv, z in zip(yy, Z): b.text(t.loc[z, 'media'] + 3, yv, f"{t.loc[z, 'razao_fundo']:.1f}×".replace('.', ','), va='center', fontsize=6.5)
b.set_yticks(yy); b.set_yticklabels([NOME[z] for z in Z], fontsize=7); b.set_xlabel('NO$_2$ médio (µmol m$^{-2}$)', fontsize=7.5)
b.set_xlim(0, 205); b.legend(frameon=False, fontsize=6.5, loc='lower right'); b.xaxis.set_major_formatter(V)
b.text(-0.62, 1.0, '(b)', transform=b.transAxes, fontsize=8, fontweight='bold', va='top')
fig.savefig(O + 'fig2_no2_mapa.png'); plt.close()

# ======== Fig 3: sazonalidade e série mensal + carga do porto
S = pd.read_csv(O + 'sazonal.csv', index_col=0); M = pd.read_csv(O + 'no2_mensal.csv', index_col=0, parse_dates=True)
SEL = ['Sao_Paulo_centro', 'Cubatao_polo', 'Santos_Porto', 'Praia_Grande', 'Peruibe', 'Oceano_fundo']
fig, ax = plt.subplots(1, 2, figsize=(16 * CM, 6.2 * CM), gridspec_kw=dict(width_ratios=[1, 1.6]))
for z in SEL: ax[0].plot(S.index, S[z], '-o', ms=2.5, lw=1.1, color=COR[z], label=NOME[z])
ax[0].set_xticks(range(1, 13)); ax[0].set_xticklabels('JFMAMJJASOND'); ax[0].set_ylabel('NO$_2$ (µmol m$^{-2}$)'); ax[0].yaxis.set_major_formatter(V)
ax[0].set_title('(a)', loc='left', fontsize=8.5, fontweight='bold')
for z in SEL: ax[1].plot(M.index, M[z], lw=0.9, color=COR[z])
ax[1].axvspan(pd.Timestamp('2020-04-01'), pd.Timestamp('2020-05-31'), color='#ccc', alpha=0.6, lw=0)
ax[1].text(pd.Timestamp('2020-05-01'), M[SEL].max().max() * 0.98, 'abr–mai\n2020', ha='center', va='top', fontsize=6)
ax[1].set_ylabel('NO$_2$ mensal (µmol m$^{-2}$)'); ax[1].yaxis.set_major_formatter(V); ax[1].set_title('(b)', loc='left', fontsize=8.5, fontweight='bold')
fig.legend(*ax[0].get_legend_handles_labels(), loc='lower center', ncol=6, frameon=False, fontsize=7)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(O + 'fig3_series.png'); plt.close()

# ======== Fig 4: ciclo semanal
W = pd.read_csv(O + 't_semanal.csv', index_col=0)
fig, ax = plt.subplots(1, 2, figsize=(16 * CM, 6 * CM), gridspec_kw=dict(width_ratios=[1.3, 1]))
for z in ['Sao_Paulo_centro', 'Cubatao_polo', 'Santos_Porto', 'Peruibe', 'Oceano_fundo']:
    ax[0].plot(range(7), [W.loc[z, f'd{k}'] for k in range(7)], '-o', ms=3, lw=1.1, color=COR[z], label=NOME[z])
ax[0].axhline(1, color='k', lw=0.5, ls=':'); ax[0].set_xticks(range(7)); ax[0].set_xticklabels(['seg', 'ter', 'qua', 'qui', 'sex', 'sáb', 'dom'])
ax[0].set_ylabel('NO$_2$ relativo à média dos dias úteis'); ax[0].yaxis.set_major_formatter(V); ax[0].legend(frameon=False, fontsize=6.5, loc='lower left')
ax[0].set_title('(a)', loc='left', fontsize=8.5, fontweight='bold')
yy = np.arange(len(Z))[::-1]
ax[1].errorbar(W.loc[Z, 'dom_util'], yy, xerr=[W.loc[Z, 'dom_util'] - W.loc[Z, 'ic_lo'], W.loc[Z, 'ic_hi'] - W.loc[Z, 'dom_util']],
               fmt='o', ms=4, color='k', ecolor='#888', capsize=2, lw=0.8)
for yv, z in zip(yy, Z): ax[1].plot(W.loc[z, 'dom_util'], yv, 'o', ms=4.5, color=COR[z], zorder=5)
ax[1].axvline(1, color='k', lw=0.5, ls=':'); ax[1].set_yticks(yy); ax[1].set_yticklabels([NOME[z] for z in Z], fontsize=7)
ax[1].set_xlabel('Razão domingo / dias úteis (IC 95%)'); ax[1].xaxis.set_major_formatter(V); ax[1].set_xlim(0.4, 1.05)
ax[1].set_title('(b)', loc='left', fontsize=8.5, fontweight='bold')
fig.tight_layout(); fig.savefig(O + 'fig4_semanal.png'); plt.close()

# ======== Fig 5: vento (velocidade e direção) em Santos/Porto e Cubatão
x = pd.read_csv(O + 'no2_diario_met.csv', parse_dates=['data']); x = x[(x.ano >= 2019) & (x.ano <= 2025)].dropna(subset=['vel'])
fig = plt.figure(figsize=(16 * CM, 6.4 * CM))
a0 = fig.add_axes([0.07, 0.17, 0.32, 0.75])
cls = [0, 1, 2, 3, 4, 10]; lab = ['<1', '1–2', '2–3', '3–4', '>4']
for z in ['Cubatao_polo', 'Santos_Porto', 'Praia_Grande', 'Peruibe']:
    g = x[x.zona == z]; m_ = g.groupby(pd.cut(g.vel, cls), observed=True).no2.mean()
    a0.plot(range(5), m_.values, '-o', ms=3, color=COR[z], label=NOME[z])
a0.set_xticks(range(5)); a0.set_xticklabels(lab); a0.set_xlabel('Velocidade do vento a 10 m (m s$^{-1}$)'); a0.set_ylabel('NO$_2$ (µmol m$^{-2}$)')
a0.legend(frameon=False, fontsize=6.5); a0.yaxis.set_major_formatter(V); a0.set_title('(a)', loc='left', fontsize=8.5, fontweight='bold')
setores = np.radians(np.arange(0, 360, 45))
for k, z in enumerate(['Santos_Porto', 'Cubatao_polo']):
    ap = fig.add_axes([0.45 + k * 0.28, 0.1, 0.25, 0.78], projection='polar'); g = x[x.zona == z]
    s_ = (np.round(g.dir / 45) % 8).astype(int); m_ = g.groupby(s_).no2.mean().reindex(range(8)); f_ = g.groupby(s_).size().reindex(range(8)) / len(g)
    ap.set_theta_zero_location('N'); ap.set_theta_direction(-1)
    bars = ap.bar(setores, m_.values, width=np.radians(40), color=plt.cm.YlOrRd(f_.values / f_.max() * 0.8 + 0.15), edgecolor='k', lw=0.4)
    ap.set_xticks(setores); ap.set_xticklabels(['N', 'NE', 'L', 'SE', 'S', 'SO', 'O', 'NO'], fontsize=6.5)
    rmax = np.nanmax(m_.values); tk = [50, 100, 150, 200]; ap.set_rticks([t_ for t_ in tk if t_ < rmax * 1.05]); ap.set_rlim(0, rmax * 1.08)
    ap.tick_params(axis='y', labelsize=5.5); ap.set_rlabel_position(200)
    ap.set_title(f"({'bc'[k]}) {NOME[z]}", fontsize=8, fontweight='bold', pad=10)
fig.text(0.72, 0.02, 'Raio: NO$_2$ médio (µmol m$^{-2}$) por direção de origem do vento; cor: frequência do setor', ha='center', fontsize=6.3)
fig.savefig(O + 'fig5_vento.png'); plt.close()

# ======== Fig 6: pandemia (mapa da diferença + barras bruto x normalizado)
A, ext = tif(D + 'no2_abrmai2020_menos_referencia.tif'); A = A * 1e6
P = pd.read_csv(O + 't_pandemia.csv', index_col=0).loc[Z]
fig = plt.figure(figsize=(16 * CM, 7.6 * CM)); ax = fig.add_axes([0.04, 0.1, 0.46, 0.86])
h = ax.imshow(A, extent=ext, cmap='RdBu_r', norm=TwoSlopeNorm(0, -120, 40), interpolation='bilinear')
ax.contour(np.where(mar, 1, 0), levels=[0.5], extent=zext, origin='upper', colors='#333', linewidths=0.5, zorder=2)
base_mapa(ax); grade(ax, [-47.0, -46.5, -46.0], [-24.4, -24.0, -23.6], 6.3)
cb = fig.colorbar(h, ax=ax, fraction=0.04, pad=0.02); cb.set_label('ΔNO$_2$ (µmol m$^{-2}$)', fontsize=7); cb.ax.tick_params(labelsize=6.5)
ax.text(0.03, 0.97, '(a)', transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')
b = fig.add_axes([0.75, 0.15, 0.23, 0.81]); yy = np.arange(len(Z))[::-1]
b.barh(yy + 0.18, P.bruto_pct, height=0.34, color='#bbb', label='Bruta')
b.barh(yy - 0.18, P.normalizado_pct, height=0.34, color=[COR[z] for z in Z], label='Normalizada')
b.axvline(0, color='k', lw=0.6); b.set_yticks(yy); b.set_yticklabels([NOME[z] for z in Z], fontsize=7)
b.set_xlabel('Variação do NO$_2$ (%)', fontsize=7.5); b.legend(frameon=False, fontsize=6.5, loc='lower left'); b.xaxis.set_major_formatter(V)
b.text(-0.62, 1.0, '(b)', transform=b.transAxes, fontsize=8, fontweight='bold', va='top')
fig.savefig(O + 'fig6_pandemia.png'); plt.close()

# ======== Fig 7: NO2 anual x carga do porto
CARGA = {2018: 133.2, 2019: 134.0, 2020: 146.6, 2021: 147.0, 2022: 162.4, 2023: 173.3, 2024: 179.8, 2025: 186.4}
An = pd.read_csv(O + 'no2_anual.csv', index_col=0)
fig, ax = plt.subplots(figsize=(12 * CM, 7.2 * CM))
for z in ['Cubatao_polo', 'Santos_Porto', 'Sao_Paulo_centro', 'Oceano_fundo']:
    ax.plot(An.index, An[z] / An.loc[2019, z] * 100, '-o', ms=3, color=COR[z], label=f'NO$_2$ {NOME[z]}')
ax.plot(list(CARGA)[1:], [CARGA[k] / CARGA[2019] * 100 for k in list(CARGA)[1:]], 's--', ms=3.5, color='k', label='Carga do porto')
ax.axhline(100, color='k', lw=0.4, ls=':'); ax.set_ylabel('Índice (2019 = 100)'); ax.yaxis.set_major_formatter(V)
ax.legend(frameon=False, fontsize=6.3, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.12)); fig.tight_layout(); fig.savefig(O + 'fig7_porto.png'); plt.close()

# ======== Fig 8: aerossóis (mapa médio + climatologia mensal + série anual)
A, ext = tif(D + 'aod_media_2019_2025.tif')
# máscara de terra: reamostra o SRTM para a grade do AOD
xs_ = np.linspace(ext[0], ext[1], A.shape[1], endpoint=False) + (ext[1] - ext[0]) / A.shape[1] / 2
ys_ = np.linspace(ext[3], ext[2], A.shape[0], endpoint=False) - (ext[3] - ext[2]) / A.shape[0] / 2
ci = np.clip(((xs_ - zext[0]) / sx).astype(int), 0, Zs.shape[1] - 1); ri = ((zext[3] - ys_) / sy).astype(int)
terra = np.zeros_like(A, bool); ok_ = (ri >= 0) & (ri < Zs.shape[0])
terra[ok_] = ~mar[ri[ok_]][:, ci]; terra[~ok_] = True
A = np.where(terra, A, np.nan); ax0_bg = '#e8eef3'
AM = pd.read_csv(O + 'aod_mensal.csv', index_col=0, parse_dates=True)
fig = plt.figure(figsize=(16 * CM, 6.6 * CM)); ax = fig.add_axes([0.04, 0.1, 0.42, 0.86])
ax.set_facecolor('#e3eef6'); h = ax.imshow(A, extent=ext, cmap='YlOrBr', vmin=0.08, vmax=0.2, interpolation='nearest')
ax.contour(np.where(mar, 1, 0), levels=[0.5], extent=zext, origin='upper', colors='#333', linewidths=0.5, zorder=2)
base_mapa(ax); grade(ax, [-47.0, -46.5, -46.0], [-24.4, -24.0, -23.6], 6)
cb = fig.colorbar(h, ax=ax, fraction=0.04, pad=0.02); cb.set_label('AOD 550 nm', fontsize=7); cb.ax.tick_params(labelsize=6.5)
ax.text(0.03, 0.97, '(a)', transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')
b = fig.add_axes([0.60, 0.58, 0.38, 0.38]); cl = AM.groupby(AM.index.month).mean()
for z in ['Cubatao_polo', 'Santos_Porto', 'Peruibe', 'Sao_Paulo_centro']: b.plot(cl.index, cl[z], '-o', ms=2.5, lw=1, color=COR[z], label=NOME[z])
b.set_xticks(range(1, 13)); b.set_xticklabels('JFMAMJJASOND', fontsize=6.5); b.set_ylabel('AOD', fontsize=7); b.yaxis.set_major_formatter(V)
b.legend(frameon=False, fontsize=5.8, ncol=2, loc='upper left'); b.text(-0.2, 1.05, '(b)', transform=b.transAxes, fontsize=8, fontweight='bold')
c = fig.add_axes([0.60, 0.1, 0.38, 0.36]); an = AM.groupby(AM.index.year).mean()
for z in ['Cubatao_polo', 'Santos_Porto', 'Peruibe', 'Sao_Paulo_centro']: c.plot(an.index, an[z], lw=1, color=COR[z])
c.set_ylabel('AOD anual', fontsize=7); c.yaxis.set_major_formatter(V); c.tick_params(labelsize=6.5)
c.text(-0.2, 1.05, '(c)', transform=c.transAxes, fontsize=8, fontweight='bold')
fig.savefig(O + 'fig8_aod.png'); plt.close()
print('ok')
