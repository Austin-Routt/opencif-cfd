"""pubstyle.py: one figure style for the dissertation, in the manner of peer-reviewed journals.

White background, a closed black frame with inward ticks on all four sides (minor ticks on), Liberation Sans (metric
equal to Arial; axis labels 8 pt, ticks and legends 7 to 7.5 pt, no text below 7 pt), thin lines, no titles inside the axes (the caption carries them), bold panel letters outside
the top-left corner, frameless legends, and a perceptually uniform colormap (viridis) for every continuous field.
Line colors are the Okabe-Ito set, which stays distinct for the common forms of color blindness and is told apart in
grayscale by line style and marker as well. Figures are sized to the 6.0 in text width (full; the UH college checklist's 1.5 in left margin on letter paper) or
2.9 in (half), and saved as vector PDF, with no creation date so that a rebuild is byte-identical, and as PNG at 600 dpi.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402

TEXT_W = 6.0           # in, the text width with the UH college checklist's 1.5 in left margin and 1 in right margin
HALF_W = 2.9
BLACK = '#000000'
OI = dict(blue='#0072B2', vermillion='#D55E00', green='#009E73', orange='#E69F00', sky='#56B4E9',
          purple='#CC79A7', yellow='#F0E442', grey='#7F7F7F')
SERIES = [OI['blue'], OI['vermillion'], OI['green'], OI['orange']]
STYLES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'D']
CMAP = 'viridis'
# One color per entity across every figure of the dissertation (validated with the Okabe-Ito order: worst adjacent
# colour-blind separation dE 11.0, normal vision 24.2; all three pairs of blue, vermillion and green: dE 11.0 and 18.7;
# the orange is kept off thin lines, 2.2:1 on white). As in the figures of Chapters 1 and 2, a series takes blue,
# vermillion or green, and black is kept for a single curve, an analytic reference or a frame (9 November 2026: V1 and
# run C, black before, are now blue and green).
FCOL = {0.3: OI['blue'], 0.5: OI['vermillion'], 0.7: OI['green']}            # filtrate fraction F
FSTY = {0.3: ('-', 'o'), 0.5: ('--', 's'), 0.7: ('-.', '^')}
ACOL = {0.80: OI['green'], 1.40: OI['blue'], 2.69: OI['vermillion']}        # aspect ratio of the CFD runs C, B, A
ASTY = {0.80: ('-', 'o'), 1.40: ('--', 's'), 2.69: ('-.', '^')}
WCOL = {100: OI['blue'], 125: OI['vermillion'], 150: OI['green']}           # Gifford et al.'s calibration widths (µm)
VCOL = {'V1': OI['blue'], 'V3': OI['vermillion']}                           # Dinh et al.'s devices


def apply():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
        'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'legend.fontsize': 7.5, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans',
        'mathtext.it': 'Liberation Sans:italic', 'mathtext.bf': 'Liberation Sans:bold',
        'axes.linewidth': 0.6, 'axes.edgecolor': BLACK, 'axes.labelcolor': BLACK, 'axes.facecolor': 'white',
        'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'text.color': BLACK,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'xtick.minor.visible': True, 'ytick.minor.visible': True,
        'xtick.major.size': 3.0, 'ytick.major.size': 3.0, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
        'xtick.color': BLACK, 'ytick.color': BLACK,
        'lines.linewidth': 1.0, 'lines.markersize': 4.0, 'axes.grid': False,
        'legend.frameon': False, 'legend.handlelength': 2.2, 'legend.borderaxespad': 0.3,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 600, 'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02, 'image.cmap': CMAP})


def letter(ax, s, dx=-0.02, dy=1.02):
    """Write the bold panel letter just outside the axes' top-left corner."""
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')


def save(fig, path_stem):
    fig.savefig(path_stem + '.pdf', metadata={'CreationDate': None})
    fig.savefig(path_stem + '.png', dpi=600)
    plt.close(fig)
