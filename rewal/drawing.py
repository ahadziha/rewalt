"""
Drawing backends.
"""

from abc import ABC

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

matplotlib.use('TkCairo')


DEFAULT = {
            'bgcolor': 'white',
            'fgcolor': 'black',
            'orientation': 'bt',
            }


class DrawBackend(ABC):
    def __init__(self, **params):
        self.bgcolor = params.get(
                'bgcolor', DEFAULT['bgcolor'])
        self.fgcolor = params.get(
                'fgcolor', DEFAULT['fgcolor'])
        self.orientation = params.get(
                'orientation', DEFAULT['orientation'])
        self.name = params.get('name', None)

    def rotate(self, xy):
        if self.orientation == 'tb':
            return (xy[0], 1-xy[1])
        if self.orientation == 'lr':
            return (xy[1], 1-xy[0])
        if self.orientation == 'rl':
            return (1-xy[1], 1-xy[0])
        return xy


class TikZBackend(DrawBackend):
    """
    TikZ drawing backend.
    """
    def __init__(self, **params):
        super().__init__(**params)
        self.bg = '\\path[fill, color={}] (0, 0) rectangle (1, 1)'.format(
                self.bgcolor)
        self.wirelayer = []
        self.nodelayer = []
        self.labellayer = []

    def draw_wire(self, wire_xy, node_xy, **params):
        """
        Draws a wire from a wire vertex to a node vertex.
        """
        color = params.get('color', self.fgcolor)
        alpha = params.get('alpha', 1)

        width = .02

        def to_cubic(p0, p1, p2):
            control1 = (p0[0]/3 + 2*p1[0]/3, p0[1]/3 + 2*p1[1]/3)
            control2 = (2*p1[0]/3 + p2[0]/3, 2*p1[1]/3 + p2[1]/3)
            return p0, control1, control2, p2

        contour = '\\path[fill, color={}] {} .. controls {} and {} .. {} '\
            'to {} .. controls {} and {} .. {};\n'.format(
                    self.bgcolor,
                    *[self.rotate(p) for p in to_cubic(
                        node_xy,
                        (wire_xy[0] - (width/2), node_xy[1]),
                        (wire_xy[0] - (width/2), wire_xy[1])
                        )],
                    *[self.rotate(p) for p in to_cubic(
                        (wire_xy[0] + (width/2), wire_xy[1]),
                        (wire_xy[0] + (width/2), node_xy[1]),
                        node_xy)]
                   )
        wire = '\\draw[color={}, opacity={}] {} .. controls {} and {} .. '\
            '{};\n'.format(
                    color,
                    alpha,
                    *[self.rotate(p) for p in to_cubic(
                        node_xy,
                        (wire_xy[0], node_xy[1]),
                        wire_xy)]
                   )
        self.wirelayer.append(contour)
        self.wirelayer.append(wire)

    def draw_label(self, label, xy, offset, **params):
        color = params.get('color', self.fgcolor)

        xy = self.rotate(xy)
        label = '\\node[text={}, font={{\\scriptsize \\sffamily}}, '\
            'xshift={}pt, yshift={}pt] at {} {{{}}};\n'.format(
                color,
                offset[0],
                offset[1],
                xy,
                label)
        self.labellayer.append(label)

    def draw_node(self, xy, **params):
        color = params.get('color', self.fgcolor)
        stroke = params.get('stroke', color)

        xy = self.rotate(xy)
        node = '\\node[circle, fill={}, draw={}, inner sep=1pt] '\
            'at {} {{}};\n'.format(
                    color, stroke, xy)
        self.nodelayer.append(node)

    def output(self, path=None, show=True):
        lines = [
                '\\begin{tikzpicture}[scale=3]\n',
                '\\path[fill={}] (0, 0) rectangle (1, 1);\n'.format(
                    self.bgcolor),
                *self.wirelayer,
                *self.nodelayer,
                *self.labellayer,
                '\\end{tikzpicture}']
        if path is None and show:
            print(''.join(lines))
        if path is not None:
            with open(path, 'w+') as file:
                file.writelines(lines)


class MatBackend(DrawBackend):
    """
    Matplotlib drawing backend.
    """
    def __init__(self, **params):
        super().__init__(**params)

        self.fig, self.axes = plt.subplots()
        self.axes.set_facecolor(self.bgcolor)
        self.axes.set_xlim(0, 1)
        self.axes.set_ylim(0, 1)
        for side in ('top', 'right', 'bottom', 'left'):
            self.axes.spines[side].set_visible(False)

    def draw_wire(self, wire_xy, node_xy, **params):
        """
        Draws a wire from a wire vertex to a node vertex.
        """
        color = params.get('color', self.fgcolor)
        alpha = params.get('alpha', 1)

        width = .02
        contour = Path(
                [
                    self.rotate(node_xy),
                    self.rotate(
                        (wire_xy[0] - (width/2), node_xy[1])),
                    self.rotate(
                        (wire_xy[0] - (width/2), wire_xy[1])),
                    self.rotate(
                        (wire_xy[0] + (width/2), wire_xy[1])),
                    self.rotate(
                        (wire_xy[0] + (width/2), node_xy[1])),
                    self.rotate(node_xy)
                ], [
                    Path.MOVETO,
                    Path.CURVE3,
                    Path.CURVE3,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.CURVE3
                ])
        wire = Path(
                [
                    self.rotate(wire_xy),
                    self.rotate(
                        (wire_xy[0], node_xy[1])),
                    self.rotate(node_xy)
                ], [
                    Path.MOVETO,
                    Path.CURVE3,
                    Path.CURVE3
                ])
        p_contour = PathPatch(
                contour,
                facecolor=self.bgcolor,
                edgecolor='none')
        p_wire = PathPatch(
                wire,
                facecolor='none',
                edgecolor=color,
                alpha=alpha,
                lw=1)
        self.axes.add_patch(p_contour)
        self.axes.add_patch(p_wire)

    def draw_label(self, label, xy, offset, **params):
        color = params.get('color', self.fgcolor)
        ha = params.get('ha', 'left')
        va = params.get('va', 'baseline')

        xy = self.rotate(xy)
        xytext = (xy[0] + offset[0], xy[1] + offset[1])
        self.axes.annotate(
                label,
                xy,
                xytext=xytext,
                textcoords='offset pixels',
                color=color,
                ha=ha,
                va=va)

    def draw_node(self, xy, **params):
        color = params.get('color', self.fgcolor)
        stroke = params.get('stroke', color)

        xy = self.rotate(xy)
        self.axes.scatter(
                xy[0],
                xy[1],
                s=40,
                c=color,
                edgecolors=stroke)

    def draw_arrow(self, xy0, xy1, **params):
        color = params.get('color', self.fgcolor)
        shorten = params.get('shorten', 1)

        xy0 = self.rotate(xy0)
        xy1 = self.rotate(xy1)
        dxy = (xy1[0] - xy0[0], xy1[1] - xy0[1])
        self.axes.arrow(
                xy0[0] + ((1 - shorten)/2)*dxy[0],
                xy0[1] + ((1 - shorten)/2)*dxy[1],
                dxy[0]*shorten,
                dxy[1]*shorten,
                fc=color,
                ec=color,
                overhang=0.8,
                lw=0.5,
                head_width=0.01,
                head_length=0.01,
                length_includes_head=True)

    def output(self, path=None, show=True):
        self.fig.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        self.fig.canvas.manager.set_window_title(self.name)
        if path is None and show:
            self.fig.show()
        if path is not None:
            self.fig.savefig(path)
            plt.close(self.fig)
