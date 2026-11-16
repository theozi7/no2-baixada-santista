import numpy as np, pandas as pd, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Carlito', 'Liberation Sans'], 'font.size': 8.5,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400, 'axes.linewidth': 0.6})
im = Image.open('../gee/no2_media_2019_2025.tif'); A = np.array(im, float) * 1e6
sx, sy = im.tag_v2[33550][:2]; x0, y0 = im.tag_v2[33922][3:5]; ny, nx = A.shape
ext = [x0, x0 + nx * sx, y0 - ny * sy, y0]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={'width_ratios': [1.15, 1]})
a = ax[0]; h = a.imshow(A, extent=ext, cmap='YlOrRd', vmin=20, vmax=160, origin='upper')
gj = json.load(open('../mapas/municipios_area.geojson'))
for f in gj['features']:
    g = f['geometry']
    if 'coordinates' not in g: continue
    polys = g['coordinates'] if g['type'] == 'MultiPolygon' else [g['coordinates']]
    for p in polys:
        r = np.array(p[0]); a.plot(r[:, 0], r[:, 1], color='0.25', lw=0.35)
for (lo, la, nm, c) in ((-46.30, -24.00, 'Barra de Santos', '#1f4e79'), (-46.50, -24.05, 'Praia Grande', '#2e8b57')):
    a.plot([lo, lo], [la, max(la - 45 / 111.32, -24.435)], color=c, lw=1.6); a.plot(lo, la, 'o', color=c, ms=3)
for (l0, l1, nm, c) in ((-46.35, -46.20, 'Fundeio', '#1f4e79'), (-46.55, -46.40, 'Mar em frente\na Praia Grande', '#2e8b57'), (-46.15, -46.00, 'Mar a leste\nde Guarujá', '#7f7f7f')):
    a.add_patch(Rectangle((l0, -24.15), l1 - l0, 0.10, fill=False, ec=c, lw=1.1, ls='--'))
    a.text((l0 + l1) / 2 + 0.035, -24.165, {'Fundeio': '50', 'Mar em frente\na Praia Grande': '64'}.get(nm, '34'), ha='center', va='top', fontsize=7.5, color=c, fontweight='bold')
a.set_xlim(-46.85, -45.9); a.set_ylim(-24.445, -23.75)
a.set_xticks([-46.8, -46.4, -46.0]); a.set_xticklabels(['46,8°O', '46,4°O', '46,0°O'])
a.set_yticks([-24.4, -24.0, -23.8]); a.set_yticklabels(['24,4°S', '24,0°S', '23,8°S'])
cb = fig.colorbar(h, ax=a, fraction=0.045, pad=0.02); cb.set_label('NO$_2$ (µmol m$^{-2}$)', fontsize=7.5); cb.ax.tick_params(labelsize=7)
a.set_title('(a)', loc='left', fontweight='bold')
P = pd.read_csv('perfis_oceano.csv', index_col=0)
b = ax[1]
b.plot(P.index, P.iloc[:, 0], '-o', ms=2.5, color='#1f4e79', label='Barra de Santos (canal do porto)')
b.plot(P.index, P.iloc[:, 1], '-s', ms=2.5, color='#2e8b57', label='Praia Grande (costa sem porto)')
b.set_xlabel('Distância da costa para o sul (km)'); b.set_ylabel('NO$_2$ (µmol m$^{-2}$)')
b.axhline(27.3, color='0.5', lw=0.6, ls=':'); b.text(1, 28.3, 'fundo oceânico', ha='left', fontsize=7, color='0.4')
b.legend(frameon=False, fontsize=7); b.set_title('(b)', loc='left', fontweight='bold')
fig.tight_layout(); fig.savefig('fig_fundeio.png')
