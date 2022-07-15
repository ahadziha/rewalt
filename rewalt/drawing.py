"""
Drawing backends.
"""

from abc import ABC

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch


DEFAULT = {
            'bgcolor': 'white',
            'fgcolor': 'black',
            'orientation': 'bt',
            }


class DrawBackend(ABC):
    """
    Abstract drawing backend for placing nodes, wires, arrows,
    and labels on a canvas.

    The purpose of this class is simply to describe the signature
    of methods that subclasses have to implement.

    Keyword arguments
    -----------------
    bgcolor : multiple types
        The background colour (default is :code:`'white'`).
    fgcolor : multiple types
        The foreground colour (default is :code:`'black'`).
    orientation : :class:`str`
        Orientation: one of :code:`'bt'` (bottom-to-top), :code:`'lr'`
        (left-to-right), :code:`'tb'` (top-to-bottom), :code:`'rl'`
        (right-to-left) (default is :code:`'bt'`).

    Notes
    -----
    All coordinates should be passed to the backend *as if* the
    orientation was bottom-to-top; the backend will then make rotations
    and adjustments according to the chosen orientation.
    """
    def __init__(self, **params):
        self.bgcolor = params.get(
                'bgcolor', DEFAULT['bgcolor'])
        self.fgcolor = params.get(
                'fgcolor', DEFAULT['fgcolor'])
        self.orientation = params.get(
                'orientation', DEFAULT['orientation'])
        self.name = params.get('name', None)

    def draw_wire(self, wire_xy, node_xy, **params):
        """
        Draws a wire from a wire vertex to a node vertex on the canvas.

        Arguments
        ---------
        wire_xy : :class:`tuple[float]`
            The coordinates of the wire vertex.
        node_xy : :class:`tuple[float]`
            The coordinates of the node vertex.

        Keyword arguments
        -----------------
        color : multiple types
            The colour of the wire (default is :code:`self.fgcolor`).
        alpha : :class:`float`
            Alpha factor of the wire (default is :code:`1`).
        depth : :class:`bool`
            Whether to draw the wire with a contour, to simulate "crossing
            over" objects that are already on the canvas (default is
            :code:`True`).
        """
        pass

    def draw_label(self, label, xy, offset, **params):
        """
        Draws a label next to a location on the canvas.

        Arguments
        ---------
        label : :class:`str`
            The label.
        xy : :class:`tuple[float]`
            The coordinates of the object to be labelled.
        offset : :class:`tuple[float]`
            Point offset of the label relative to the object.

        Keyword arguments
        -----------------
        color : multiple types
            The colour of the label (default is :code:`self.fgcolor`).
        """
        pass

    def draw_node(self, xy, **params):
        """
        Draws a node on the canvas.

        Arguments
        ---------
        xy : :class:`tuple[float]`
            The coordinates of the node.

        Keyword arguments
        -----------------
        color : multiple types
            Fill colour of the node (default is :code:`self.fgcolor`).
        stroke : multiple types
            Stroke colour of the node (default is same as `color`).
        """
        pass

    def draw_arrow(self, xy0, xy1, **params):
        """
        Draws an arrow on the canvas.

        Arguments
        ---------
        xy0 : :class:`tuple[float]`
            The coordinates of the starting point.
        xy1 : :class:`tuple[float]`
            The coordinates of the ending point.

        Keyword arguments
        -----------------
        color : multiple types
            Colour of the arrow (default is :code:`self.fgcolor`).
        shorten : :class:`float`
            Factor by which to scale the length (default is :code:`1`).
        """
        pass

    def output(self, **params):
        """
        Output the picture.

        Keyword arguments
        -----------------
        show : :class:`bool`
            Whether to show the output (default is :code:`True`).
        path : :class:`str`
            Path where to save the output (default is :code:`None`).
        scale : :class:`float`
            (TikZ only) Scale factor to apply to output (default is
            :code:`3`).
        xscale : :class:`float`
            (TikZ only) Scale factor to apply to x axis in output
            (default is same as `scale`)
        yscale : :class:`float`
            (TikZ only) Scale factor to apply to y axis in output
            (default is same as `scale`)
        """
        pass

    def rotate(self, xy):
        """
        Returns coordinates rotated according to the orientation
        of the picture.

        Arguments
        ---------
        xy : :class:`tuple[float]`
            The coordinates to rotate.

        Returns
        -------
        rotate : :class:`tuple[float]`
            The rotated coordinates.
        """
        if self.orientation == 'tb':
            return (xy[0], 1-xy[1])
        if self.orientation == 'lr':
            return (xy[1], 1-xy[0])
        if self.orientation == 'rl':
            return (1-xy[1], 1-xy[0])
        return xy


class TikZBackend(DrawBackend):
    """
    Drawing backend outputting TikZ code that can be embedded in a
    LaTeX document.
    """
    def __init__(self, **params):
        super().__init__(**params)
        self.bg = '\\path[fill, color={}] (0, 0) rectangle (1, 1)'.format(
                self.bgcolor)
        self.wirelayer = []
        self.nodelayer = []
        self.arrowlayer = []
        self.labellayer = []

    def draw_wire(self, wire_xy, node_xy, **params):
        super().draw_wire(wire_xy, node_xy, **params)
        color = params.get('color', self.fgcolor)
        alpha = params.get('alpha', 1)
        depth = params.get('depth', True)

        def to_cubic(p0, p1, p2):
            control1 = (p0[0]/3 + 2*p1[0]/3, p0[1]/3 + 2*p1[1]/3)
            control2 = (2*p1[0]/3 + p2[0]/3, 2*p1[1]/3 + p2[1]/3)
            return p0, control1, control2, p2

        if depth:
            width = .02
            contour = '\\path[fill, color={}] {} .. controls {} '\
                'and {} .. {} to {} .. controls {} and {} .. {};\n'.format(
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
            self.wirelayer.append(contour)

        wire = '\\draw[color={}, opacity={}] {} .. controls {} and {} .. '\
            '{};\n'.format(
                    color,
                    alpha,
                    *[self.rotate(p) for p in to_cubic(
                        node_xy,
                        (wire_xy[0], node_xy[1]),
                        wire_xy)]
                   )
        self.wirelayer.append(wire)

    def draw_label(self, label, xy, offset, **params):
        super().draw_label(label, xy, offset, **params)
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
        super().draw_node(xy, **params)
        color = params.get('color', self.fgcolor)
        stroke = params.get('stroke', color)

        xy = self.rotate(xy)
        node = '\\node[circle, fill={}, draw={}, inner sep=1pt] '\
            'at {} {{}};\n'.format(
                    color, stroke, xy)
        self.nodelayer.append(node)

    def draw_arrow(self, xy0, xy1, **params):
        super().draw_arrow(xy0, xy1, **params)
        color = params.get('color', self.fgcolor)
        shorten = params.get('shorten', 1)

        xy0 = self.rotate(xy0)
        xy1 = self.rotate(xy1)
        dxy = (xy1[0] - xy0[0], xy1[1] - xy0[1])

        xy0_off = (
                xy0[0] + 0.5*(1 - shorten)*dxy[0],
                xy0[1] + 0.5*(1 - shorten)*dxy[1])
        xy1_off = (
                xy1[0] - 0.5*(1 - shorten)*dxy[0],
                xy1[1] - 0.5*(1 - shorten)*dxy[1])

        arrow = '\\draw[->, draw={}] {} -- {};\n'.format(
                    color,
                    xy0_off,
                    xy1_off)
        self.arrowlayer.append(arrow)

    def output(self, **params):
        super().output(**params)
        path = params.get('path', None)
        show = params.get('show', True)

        scale = params.get('scale', 1)
        xscale = params.get('xscale', scale)
        yscale = params.get('yscale', scale)

        baseline = '{([yshift=-.5ex]current bounding box.center)}'
        lines = [
                '\\begin{{tikzpicture}}[xscale={}, yscale={}, '
                'baseline={}]\n'.format(
                    xscale,
                    yscale,
                    baseline),
                '\\path[fill={}] (0, 0) rectangle (1, 1);\n'.format(
                    self.bgcolor),
                *self.wirelayer,
                *self.nodelayer,
                *self.arrowlayer,
                *self.labellayer,
                '\\end{tikzpicture}']
        if path is None and show:
            print(''.join(lines))
        if path is not None:
            with open(path, 'w+') as file:
                file.writelines(lines)


class MatBackend(DrawBackend):
    """
    Drawing backend outputting Matplotlib figures.
    """
    def __init__(self, **params):
        super().__init__(**params)

        self.fig, self.axes = plt.subplots()
        self.axes.set_facecolor(self.bgcolor)
        self.axes.set_xlim(0, 1)
        self.axes.set_ylim(0, 1)

        self.axes.xaxis.set_visible(False)
        self.axes.yaxis.set_visible(False)
        for side in ('top', 'right', 'bottom', 'left'):
            self.axes.spines[side].set_visible(False)

    def draw_wire(self, wire_xy, node_xy, **params):
        super().draw_wire(wire_xy, node_xy, **params)
        color = params.get('color', self.fgcolor)
        alpha = params.get('alpha', 1)
        depth = params.get('depth', True)

        if depth:
            width = .02
            contour = Path(
                [
                    self.rotate(node_xy),
                    self.rotate(
                        (wire_xy[0] - 0.5*width, node_xy[1])),
                    self.rotate(
                        (wire_xy[0] - 0.5*width, wire_xy[1])),
                    self.rotate(
                        (wire_xy[0] + 0.5*width, wire_xy[1])),
                    self.rotate(
                        (wire_xy[0] + 0.5*width, node_xy[1])),
                    self.rotate(node_xy)
                ], [
                    Path.MOVETO,
                    Path.CURVE3,
                    Path.CURVE3,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.CURVE3
                ])
            p_contour = PathPatch(
                contour,
                facecolor=self.bgcolor,
                edgecolor='none')
            self.axes.add_patch(p_contour)

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
        p_wire = PathPatch(
                wire,
                facecolor='none',
                edgecolor=color,
                alpha=alpha,
                lw=1)
        self.axes.add_patch(p_wire)

    def draw_label(self, label, xy, offset, **params):
        super().draw_label(label, xy, offset, **params)
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
        super().draw_node(xy, **params)
        color = params.get('color', self.fgcolor)
        stroke = params.get('stroke', color)

        xy = self.rotate(xy)
        self.axes.scatter(
                xy[0],
                xy[1],
                s=40,
                c=color,
                edgecolors=stroke,
                zorder=2)

    def draw_arrow(self, xy0, xy1, **params):
        super().draw_arrow(xy0, xy1, **params)
        color = params.get('color', self.fgcolor)
        shorten = params.get('shorten', 1)

        xy0 = self.rotate(xy0)
        xy1 = self.rotate(xy1)
        dxy = (xy1[0] - xy0[0], xy1[1] - xy0[1])

        xy0_off = (
                xy0[0] + 0.5*(1 - shorten)*dxy[0],
                xy0[1] + 0.5*(1 - shorten)*dxy[1])
        xy1_off = (
                xy1[0] - 0.5*(1 - shorten)*dxy[0],
                xy1[1] - 0.5*(1 - shorten)*dxy[1])
        self.axes.annotate(
                '',
                xy=xy1_off,
                xytext=xy0_off,
                arrowprops=dict(
                    arrowstyle='->',
                    color=color,
                    shrinkA=0,
                    shrinkB=0))

    def output(self, **params):
        super().output(**params)
        path = params.get('path', None)
        show = params.get('show', True)

        self.fig.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        if path is None and show:
            plt.show()
        if path is not None:
            self.fig.savefig(path)
            plt.close(self.fig)
