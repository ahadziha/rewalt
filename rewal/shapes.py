"""
Implements shapes of cells and diagrams.
"""

import networkx as nx

from rewal import utils
from rewal.ogposets import (El, OgPoset, GrSet, GrSubset, Closed,
                            OgMap, OgMapPair)


class Shape(OgPoset):
    """
    Class for shapes of cells and diagrams.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    @property
    def isatom(self):
        """
        Returns whether the shape is an atom (has a greatest element).
        """
        return len(self.maximal()) == 1

    @property
    def isround(self):
        """
        Returns whether the shape is round (has spherical boundary).
        """
        return self.all().isround

    @property
    def layers(self):
        return self.id().layers

    @property
    def rewrite_steps(self):
        return self.id().rewrite_steps

    # Main constructors
    @staticmethod
    def atom(fst, snd, **params):
        """
        Given two shapes with identical round boundaries, returns a new
        atomic shape whose input boundary is the first one and output
        boundary the second one.
        """
        cospan = params.get('cospan', False)

        for u in fst, snd:
            utils.typecheck(u, {
                'type': Shape,
                'st': lambda v: v.isround,
                'why': 'expecting a round Shape'})
        if fst.dim != snd.dim:
            raise ValueError(utils.value_err(
                snd, 'dimension does not match dimension of {}'.format(
                    repr(fst))))
        dim = fst.dim

        if dim == -1:  # Avoid more work in this simple case
            if cospan:
                return OgMapPair(
                        Shape.point().initial(),
                        Shape.point().initial()
                    )
            return Shape.point()

        in_span = OgMapPair(
                fst.boundary('-'), snd.boundary('-'))
        out_span = OgMapPair(
                fst.boundary('+'), snd.boundary('+'))
        if not in_span.isspan:
            raise ValueError(utils.value_err(
                snd, 'input boundary does not match '
                'input boundary of {}'.format(repr(fst))))
        if not out_span.isspan:
            raise ValueError(utils.value_err(
                snd, 'output boundary does not match '
                'output boundary of {}'.format(repr(fst))))

        glue_in = in_span.pushout(wfcheck=False)
        glue_out = out_span.then(glue_in).coequaliser(wfcheck=False)
        inclusion = glue_in.then(glue_out)

        sphere = inclusion.target
        # Add a greatest element
        face_data = [
                *sphere.face_data,
                [
                    {'-':
                        {inclusion.fst[x].pos for x in fst[dim]},
                     '+':
                        {inclusion.snd[x].pos for x in snd[dim]}}
                ]]
        coface_data = [
                *sphere.coface_data,
                [{'-': set(), '+': set()}]]
        for x in fst[dim]:
            coface_data[dim][inclusion.fst[x].pos]['-'].add(0)
        for x in snd[dim]:
            coface_data[dim][inclusion.snd[x].pos]['+'].add(0)

        ogatom = OgPoset(face_data, coface_data,
                         wfcheck=False, matchcheck=False)

        def inheritance():
            if isinstance(fst, OpetopeTree) and isinstance(snd, Opetope):
                if fst.dim == 0:
                    return Arrow
                if isinstance(fst, Globe):
                    return Globe
                return Opetope
            return Shape

        if cospan:
            boundary_in = OgMap(
                fst, ogatom, inclusion.fst.mapping,
                wfcheck=False)
            boundary_out = OgMap(
                snd, ogatom, inclusion.snd.mapping,
                wfcheck=False)
            atom_cospan = OgMapPair(boundary_in, boundary_out).then(
                    Shape._reorder(ogatom).inv())
            return inheritance()._upgrademaptgt(atom_cospan)

        atom = Shape._reorder(ogatom).source
        return inheritance()._upgrade(atom)

    def _atom(self, other, **params):
        return Shape.atom(self, other, **params)

    @staticmethod
    def paste(fst, snd, dim=None, **params):
        """
        Returns the pasting of two shapes along the output k-boundary
        of the first and the input k-boundary of the second.
        """
        cospan = params.get('cospan', False)

        for u in fst, snd:
            utils.typecheck(u, {'type': Shape})
        if dim is None:  # default is principal composition
            dim = min(fst.dim, snd.dim) - 1
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})

        span = OgMapPair(
                fst.boundary('+', dim),
                snd.boundary('-', dim))
        if not span.isspan:
            raise ValueError(utils.value_err(
                    snd,
                    'input {}-boundary does not match '
                    'output {}-boundary of {}'.format(
                        str(dim), str(dim), repr(fst))))

        if dim >= fst.dim:
            return OgMapPair(span.snd, snd.id()) if cospan else snd
        if dim >= snd.dim:
            return OgMapPair(fst.id(), span.fst) if cospan else fst

        pushout = span.pushout(wfcheck=False)
        paste_cospan = pushout.then(Shape._reorder(pushout.target).inv())

        def inheritance():
            if isinstance(fst, Theta) and isinstance(snd, Theta):
                if isinstance(fst, GlobeString) and \
                        isinstance(snd, GlobeString) \
                        and fst.dim == snd.dim == dim+1:
                    return GlobeString
                return Theta
            if isinstance(fst, OpetopeTree) and isinstance(snd, GlobeString) \
                    and fst.dim == snd.dim == dim+1:
                return OpetopeTree
            return Shape

        paste_cospan = inheritance()._upgrademaptgt(paste_cospan)

        if fst.dim == snd.dim == dim + 1:  # add layering
            if hasattr(fst, '_layering') and len(fst[fst.dim]) > 1:
                layering_fst = fst._layering
            else:
                layering_fst = [fst.id()]

            if hasattr(snd, '_layering') and len(snd[snd.dim]) > 1:
                layering_snd = snd._layering
            else:
                layering_snd = [snd.id()]

            layering = [
                *[f.then(paste_cospan.fst) for f in layering_fst],
                *[f.then(paste_cospan.snd) for f in layering_snd]]
            paste_cospan.target._layering = layering

        if cospan:
            return paste_cospan
        return paste_cospan.target

    def _paste(self, other, dim=None, **params):
        return Shape.paste(self, other, dim, **params)

    # Other constructors
    @staticmethod
    def paste_along(fst, snd, **params):
        """
        Returns the pasting of two shapes along the entire input
        (output) k-boundary of one, and a subshape of the output
        (input) k-boundary of the other.
        """
        wfcheck = params.get('wfcheck', True)
        cospan = params.get('cospan', False)

        span = OgMapPair(fst, snd)
        utils.typecheck(span, {
            'type': OgMapPair,
            'st': lambda x: x.isspan and x.isinjective,
            'why': 'expecting a span of injective maps'
            }, {'type': ShapeMap})

        dim = span.source.dim
        fst_image = fst.source.maximal().image(fst)
        snd_image = snd.source.maximal().image(snd)
        fst_output = fst.target.all().boundary_max('+', dim)
        snd_input = snd.target.all().boundary_max('-', dim)

        if fst_image == fst_output and snd_image == snd_input:
            return Shape.paste(fst.target, snd.target, dim,
                               cospan=cospan)

        if wfcheck:
            def condition():
                t1 = fst_image.issubset(fst_output)
                t2 = snd_image.issubset(snd_input)
                t3 = fst_image == fst_output or snd_image == snd_input
                return t1 and t2 and t3

            if not condition():
                raise ValueError(utils.value_err(
                    span, 'not a well-formed span for pasting'))

            if fst_image == fst_output:
                if not snd.target._ispastable(
                        snd_image.support, snd_input.support):
                    raise ValueError(utils.value_err(
                        snd, 'cannot paste along this map'))
            else:
                if not fst.target._ispastable(
                        fst_image.support, fst_output.support):
                    raise ValueError(utils.value_err(
                        fst, 'cannot paste along this map'))

        pushout = span.pushout(wfcheck=False)

        def inheritance():
            if isinstance(fst.target, OpetopeTree) and \
                    isinstance(snd.target, OpetopeTree) \
                    and fst.target.dim == snd.target.dim == dim+1:
                return OpetopeTree
            return Shape

        if cospan:
            paste_cospan = pushout.then(
                    Shape._reorder(pushout.target).inv())
            return inheritance()._upgrademaptgt(paste_cospan)

        paste = Shape._reorder(pushout.target).source
        return inheritance()._upgrade(paste)

    def to_outputs(self, positions, other, dim=None, **params):
        """
        Paste along the inclusion of several outputs.
        """
        if isinstance(positions, int):
            positions = [positions]
        if dim is None:
            dim = self.dim-1

        fst = GrSet(*[El(dim, pos) for pos in positions])
        snd = self.all().boundary_max('+', dim).support
        if not self._ispastable(fst, snd):
            raise ValueError(utils.value_err(
                positions, 'cannot paste to these outputs'))

        oginclusion = GrSubset(fst, self).closure().as_map
        inclusion = ShapeMap(Shape._reorder(
                oginclusion.source).then(oginclusion),
                wfcheck=False)

        other_boundary = other.boundary('-')
        if inclusion.source != other_boundary.source:
            raise ValueError(utils.value_err(
                positions, 'does not match input boundary of {}'.format(
                    repr(other))))
        return Shape.paste_along(
            inclusion,
            other_boundary, wfcheck=False, **params)

    def to_inputs(self, positions, other, dim=None, **params):
        """
        Paste along the inclusion of several outputs.
        """
        if isinstance(positions, int):
            positions = [positions]
        if dim is None:
            dim = self.dim-1

        fst = GrSet(*[El(dim, pos) for pos in positions])
        snd = self.all().boundary_max('-', dim).support
        if not self._ispastable(fst, snd):
            raise ValueError(utils.value_err(
                positions, 'cannot paste to these inputs'))

        oginclusion = GrSubset(fst, self).closure().as_map
        inclusion = ShapeMap(Shape._reorder(
                oginclusion.source).then(oginclusion),
                wfcheck=False)

        other_boundary = other.boundary('+')
        if inclusion.source != other_boundary.source:
            raise ValueError(utils.value_err(
                positions, 'does not match output boundary of {}'.format(
                    repr(other))))
        return Shape.paste_along(
            other_boundary,
            inclusion,
            wfcheck=False, **params)

    @staticmethod
    def suspend(shape, n=1):
        if n == 0:
            return shape

        if not isinstance(shape, Shape) or isinstance(shape, Empty):
            return OgPoset.suspend(shape, n)
        if isinstance(shape, Point) and n == 1:
            return Arrow()

        suspension = Shape._reorder(
                OgPoset.suspend(shape, n)).source

        def inheritance():
            if isinstance(shape, Theta):
                if isinstance(shape, GlobeString):
                    if isinstance(shape, Globe):
                        return Globe
                    return GlobeString
                return Theta
            return Shape

        return inheritance()._upgrade(suspension)

    @staticmethod
    def gray(*shapes):
        """ Gray product of shapes. """
        for x in shapes:
            if not isinstance(x, Shape):
                return OgPoset.gray(*shapes)
        if len(shapes) == 0:
            return Point()
        if len(shapes) == 1:
            return shapes[0]

        oggray = OgPoset.gray(*shapes)
        if oggray in shapes:
            return oggray

        gray = Shape._reorder(oggray).source

        def inheritance(l):
            if all([isinstance(x, Cube) for x in l]):
                return Cube
            return Shape

        return inheritance(shapes)._upgrade(gray)

    @staticmethod
    def join(*shapes):
        for x in shapes:
            if not isinstance(x, Shape):
                return OgPoset.join(*shapes)
        if len(shapes) == 0:
            return Empty()
        if len(shapes) == 1:
            return shapes[0]

        ogjoin = OgPoset.join(*shapes)
        if len(ogjoin) == 3:  # check explicitly if it's arrow
            return Arrow()
        if ogjoin in shapes:
            return ogjoin

        join = Shape._reorder(ogjoin).source

        def inheritance(l):
            if all([isinstance(x, Simplex) for x in l]):
                return Simplex
            return Shape

        return inheritance(shapes)._upgrade(join)

    @staticmethod
    def dual(shape, *dims, **params):
        reordering = params.get('reordering', False)

        reordermap = Shape._reorder(OgPoset.dual(shape, *dims))
        dual = reordermap.source
        if shape == dual:
            if reordering:
                return shape.id()
            return shape

        def inheritance():
            if isinstance(shape, Theta):
                return Theta
            return Shape

        if reordering:
            return inheritance()._upgrademapsrc(reordermap)
        return inheritance()._upgrade(dual)

    def merge(self):
        """
        Returns an atom with the same boundary as the shape (if round).
        """
        if self.isatom:
            return self
        if not self.isround:
            raise ValueError(utils.value_err(
                self, 'not a round shape'))

        merged = Shape.atom(
                self.boundary('-').source,
                self.boundary('+').source)

        def inheritance():
            if self.dim == 1:
                return Arrow
            if isinstance(self, GlobeString):
                return Globe
            if isinstance(self, OpetopeTree):
                return Opetope
            return Shape

        return inheritance()._upgrade(merged)

    # Named shapes
    @staticmethod
    def empty():
        return Empty()

    @staticmethod
    def point():
        return Point()

    @staticmethod
    def arrow():
        return Arrow()

    @staticmethod
    def simplex(dim=-1):
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= -1,
            'why': 'expecting integer >= -1'})
        point = Point()
        return Shape.join(*[point for _ in range(dim+1)])

    @staticmethod
    def cube(dim=0):
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})
        arrow = Arrow()
        return Shape.gray(*[arrow for _ in range(dim)])

    @staticmethod
    def globe(dim=0):
        """ Globes. """
        return Shape.suspend(Point(), dim)

    @staticmethod
    def theta(*thetas):
        """ Batanin cells. """
        if len(thetas) > 0:
            theta = thetas[0]
            utils.typecheck(theta, {'type': Theta})
            thetas = thetas[1:]
            if len(thetas) > 0:
                return Shape.paste(
                        Shape.suspend(theta),
                        Shape.theta(*thetas), 0)
            return Shape.suspend(theta)
        return Point()

    # Special maps
    def id(self):
        return ShapeMap(super().id(),
                        wfcheck=False)

    def boundary(self, sign=None, dim=None):
        """
        Input and output boundaries of Shapes are Shapes.
        """
        if isinstance(dim, int) and dim >= self.dim:
            return self.id()

        boundary_ogmap = super().boundary(sign, dim)
        if sign is None:
            return boundary_ogmap
        reordering = Shape._reorder(boundary_ogmap.source)
        boundary = reordering.then(boundary_ogmap)

        def inheritance():
            if dim == -1:
                return Empty
            if dim == 0:
                return Point
            if isinstance(self, OpetopeTree):
                if isinstance(self, GlobeString):
                    if dim == 1:
                        return Arrow
                    return Globe
                if utils.mksign(sign) == '+':
                    return Opetope
                return OpetopeTree
            return Shape

        return inheritance()._upgrademapsrc(boundary)

    def atom_inclusion(self, element):
        """
        Returns the inclusion of an atom as the closure of an element
        in the shape.
        """
        oginclusion = self.underset(element).as_map
        reordering = Shape._reorder(oginclusion.source)
        inclusion = reordering.then(oginclusion)

        def inheritance():
            if element.dim == 0:
                return Point
            if element.dim == 1:
                return Arrow
            if isinstance(self, Theta):
                return Globe
            if isinstance(self, OpetopeTree):
                return Opetope
            if isinstance(self, Simplex):
                return Simplex
            if isinstance(self, Cube):
                return Cube
            return Shape

        return inheritance()._upgrademapsrc(inclusion)

    def initial(self):
        """
        Returns the unique map from the initial (empty) shape.
        """
        return ShapeMap(
                OgMap(Empty(), self,
                      wfcheck=False),
                wfcheck=False)

    def terminal(self):
        """
        Returns the unique map to the point.
        """
        mapping = [
                [El(0, 0) for _ in n_data]
                for n_data in self.face_data]
        return ShapeMap(
                OgMap(self, Point(), mapping,
                      wfcheck=False),
                wfcheck=False)

    def inflate(self, collapsed=None):
        """
        Forms a cylinder on a shape with some "collapsed" sides, specified
        by a closed subset of the boundary of the shape, and returns its
        projection on the original shape.
        Used in constructing units and unitors on diagrams.
        """
        if self.dim == -1:  # Some simple cases
            return self.id()
        if self.dim == 0:
            return Shape.arrow().terminal()

        boundary_set = self.all().boundary()
        if collapsed is not None:
            utils.typecheck(collapsed, {
                'type': Closed,
                'st': lambda x: x.issubset(boundary_set),
                'why': "expecting a closed subset of the shape's boundary"})
        else:
            collapsed = boundary_set  # Default is whole boundary.

        asmap = collapsed.as_map
        arrow = Shape.arrow()
        map1 = OgMap.gray(arrow.id(), asmap)
        map2 = OgMap.gray(arrow.terminal(), asmap.source.id())
        pushout = OgMapPair(map1, map2).pushout(wfcheck=False)

        # We use the cylinder projection map and the first leg of the pushout
        # to define the projection map.
        cyl_proj = OgMap.gray(arrow.terminal(), self.id())
        collapse = pushout.fst
        mapping = [
                [None for _ in n_data]
                for n_data in pushout.target.face_data]
        for x in collapse.source:
            y = collapse[x]
            mapping[y.dim][y.pos] = cyl_proj[x]

        ogproj = OgMap(collapse.target, self, mapping,
                       wfcheck=False)
        proj = Shape._reorder(collapse.target).then(ogproj)

        def inheritance():
            if collapsed == boundary_set and isinstance(self, Opetope):
                if isinstance(self, Globe):
                    return Globe
                return Opetope
            return Shape

        return inheritance()._upgrademapsrc(proj)

    def all_layerings(self):
        """
        Returns an iterator on all layerings of a shape.
        """
        dim = self.dim
        maximal = self.maximal().support
        topdim = maximal[dim]
        flow = self._flowgraph(topdim)
        all_sorts = nx.all_topological_sorts(flow)

        def test(sort):
            return self._islayering(
                    list(sort), maximal, layerlist=True)

        def layering(layers):
            layering = []
            for layer in layers:
                oginclusion = GrSubset(
                        layer, self, wfcheck=False).closure().as_map
                inclusion = Shape._reorder(oginclusion.source).then(
                        oginclusion)
                layering.append(ShapeMap(inclusion, wfcheck=False))
            return layering

        return (
                layering(test(sort)[1])
                for sort in all_sorts if test(sort)[0]
                )

    def generate_layering(self):
        """
        Iterates through layerings.
        """
        if not hasattr(self, '_layering_gen'):
            self._layering_gen = self.all_layerings()
        try:
            self._layering = next(self._layering_gen)
        except StopIteration:
            self._layering_gen = self.all_layerings()
            self._layering = next(self._layering_gen)

        return self.layers

    def draw(self, **params):
        from rewal.strdiags import draw
        return draw(self, **params)

    def draw_boundaries(self, **params):
        from rewal.strdiags import draw_boundaries
        return draw_boundaries(self, **params)

    # Private methods
    @staticmethod
    def _reorder(shape):
        """
        Traverses all elements and returns an isomorphism
        from the shape with elements reordered in traversal order.
        """
        mapping = [[] for _ in range(shape.dim + 1)]
        marked = GrSet()

        focus_stack = [shape.maximal().support]  # traversal begins
        while len(focus_stack) > 0:
            focus = focus_stack[-1]
            if focus.issubset(marked):
                del focus_stack[-1]
            else:
                dim = focus.dim
                top_dim = focus[dim]
                if dim == 0:
                    for x in top_dim:
                        mapping[dim].append(x)
                        marked.add(x)
                    del focus_stack[-1]
                else:
                    focus_in_faces = GrSet().union(
                            *[shape.faces(x, '-') for x in top_dim])
                    focus_input = focus[:dim]
                    for y in focus_in_faces:
                        if shape.cofaces(y, '+').isdisjoint(focus):
                            focus_input.add(y)
                    if not focus_input.issubset(marked):
                        focus_stack.append(focus_input)
                    else:
                        if len(focus) == 1:
                            for x in top_dim:
                                mapping[dim].append(x)
                                marked.add(x)
                            del focus_stack[-1]
                            focus_out_faces = GrSet().union(
                                *[shape.faces(x, '+') for x in top_dim])
                            focus_output = focus[:dim]
                            for y in focus_out_faces:
                                if shape.cofaces(y, '-').isdisjoint(focus):
                                    focus_output.add(y)
                            if not focus_output.issubset(marked):
                                focus_stack.append(focus_output)
                        else:
                            def candidates(x):
                                return [y for y in shape.cofaces(x, '-')
                                        if y in focus and y not in marked]
                            x = next(
                                    x for x in mapping[dim-1]
                                    if x in focus_in_faces
                                    and len(candidates(x)) > 0)
                            focus_stack.append(
                                GrSet(candidates(x)[0]))

        def reordered_faces(x, sign):
            return {k for y in shape.faces(x, sign)
                    for k in range(shape.size[x.dim - 1])
                    if y == mapping[x.dim - 1][k]}

        def reordered_cofaces(x, sign):
            return {k for y in shape.cofaces(x, sign)
                    for k in range(shape.size[x.dim + 1])
                    if y == mapping[x.dim + 1][k]}

        face_data = [
                [
                    {sign: reordered_faces(x, sign)
                     for sign in ('-', '+')}
                    for x in n_data
                ]
                for n_data in mapping]
        coface_data = [
                [
                    {sign: reordered_cofaces(x, sign)
                     for sign in ('-', '+')}
                    for x in n_data
                ]
                for n_data in mapping]
        reordered = Shape._upgrade(
                OgPoset(face_data, coface_data,
                        wfcheck=False, matchcheck=False))

        return OgMap(reordered, shape, mapping,
                     wfcheck=False)

    def _flowgraph(self, grset):
        """
        The 'flow graph' of the set in the shape.
        """
        flowgraph = nx.DiGraph()
        flowgraph.add_nodes_from(grset)
        if grset.dim > 0:
            for x in grset:
                for y in grset:
                    if not self.faces(x, '+').isdisjoint(
                            self.faces(y, '-')):
                        flowgraph.add_edge(x, y)
        return flowgraph

    def _ispastable(self, fst, snd, **params):
        """
        Returns whether fst is a 'pastable', round region of snd
        (both given just as GrSets of maximal elements)
        """
        wfcheck = params.get('wfcheck', True)

        if wfcheck:
            if not fst.issubset(snd):
                return False
            if not fst.dim == snd.dim:
                return False
        if len(fst) == 1:
            return True
        if wfcheck:
            fst_closed = GrSubset(fst, self).closure()
            if not fst_closed.isround:
                return False
        dim = fst.dim
        fst_flow = self._flowgraph(fst)
        snd_flow = self._flowgraph(snd[dim])
        mapping = {
            x: x for x in fst}
        remaining = set(fst_flow.nodes)
        for e in fst_flow.edges:
            src = mapping[e[0]]
            tgt = mapping[e[1]]
            if src != tgt:
                snd_flow = nx.contracted_edge(
                        snd_flow, (src, tgt), self_loops=False)
                for x in mapping:
                    if mapping[x] == tgt:
                        mapping[x] = src
                remaining.remove(tgt)

        if len(remaining) > 1:  # fst_flow not connected
            return False
        if not nx.is_directed_acyclic_graph(snd_flow):  # fst_flow not convex
            return False
        if fst.dim < 3:  # nothing else needs to be checked in dim <= 2
            return True

        fst_sort = next(
                (
                    list(sort)
                    for sort in nx.all_topological_sorts(fst_flow)
                    if self._islayering(list(sort), fst)
                ), None)
        if fst_sort is None:  # cannot layer fst
            return False

        for x in remaining:
            fst_el = x
        for sort in nx.all_topological_sorts(snd_flow):
            snd_sort = list(sort)
            fst_index = snd_sort.index(fst_el)
            amended = [
                    *snd_sort[:fst_index],
                    *fst_sort,
                    *snd_sort[fst_index+1:]
                    ]
            if self._islayering(amended, snd):
                return True
        return False

    def _islayering(self, ellist, grset, **params):
        """
        Returns whether a list of top-dimensional elements is a valid
        layering of a molecule (given as its set of maximal elements)
        """
        layerlist = params.get('layerlist', False)

        dim = grset.dim
        top_dim = grset[dim]

        focus = grset[:dim]  # initialise to input boundary
        in_faces = GrSet().union(
            *[self.faces(x, '-') for x in top_dim])
        for y in in_faces:
            if self.cofaces(y, '+').isdisjoint(grset):
                focus.add(y)

        if layerlist:
            layers = []
        for x in ellist:
            in_x = self.faces(x, '-')
            if not self._ispastable(
                    in_x, focus):
                if layerlist:
                    return False, []
                return False

            if layerlist:
                layers.append(
                    focus.difference(in_x).union(GrSet(x)))
            out_x = self.faces(x, '+')
            focus = focus.difference(in_x).union(out_x)
        if layerlist:
            return True, layers
        return True

    @classmethod
    def _upgrade(cls, ogp):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        shape = OgPoset.__new__(cls)
        shape.__init__(ogp.face_data, ogp.coface_data,
                       wfcheck=False, matchcheck=False)
        return shape

    @classmethod
    def _upgrademapsrc(cls, ogmap):
        """
        Upgrades the source of an OgMap to a shape class, and declares
        it a ShapeMap.
        """
        if isinstance(ogmap, OgMapPair):
            return OgMapPair(
                    cls._upgrademapsrc(ogmap.fst),
                    cls._upgrademapsrc(ogmap.snd))
        return ShapeMap(OgMap(
                cls._upgrade(ogmap.source),
                ogmap.target,
                ogmap.mapping,
                wfcheck=False), wfcheck=False)

    @classmethod
    def _upgrademaptgt(cls, ogmap):
        """
        Upgrades the target of an OgMap to a shape class, and declares
        it a ShapeMap.
        """
        if isinstance(ogmap, OgMapPair):
            return OgMapPair(
                    cls._upgrademaptgt(ogmap.fst),
                    cls._upgrademaptgt(ogmap.snd))
        return ShapeMap(OgMap(
            ogmap.source,
            cls._upgrade(ogmap.target),
            ogmap.mapping,
            wfcheck=False), wfcheck=False)


class Simplex(Shape):
    """
    Class for oriented simplices.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    def simplex_face(self, k):
        """ Simplicial face maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim + 1),
            'why': 'out of bounds'})
        pointid = Point().id()
        maps = [
                *[pointid for _ in range(k)],
                Empty().terminal(),
                *[pointid for _ in range(self.dim - k)]
            ]
        return ShapeMap.join(*maps)

    def simplex_degeneracy(self, k):
        """ Simplicial degeneracy maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim + 1),
            'why': 'out of bounds'})
        pointid = Point().id()
        maps = [
                *[pointid for _ in range(k)],
                Arrow().terminal(),
                *[pointid for _ in range(self.dim - k)]
            ]
        return ShapeMap.join(*maps)


class Empty(Simplex):
    """
    Class for the empty shape.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    def __init__(self):
        super().__init__([], [],
                         wfcheck=False, matchcheck=False)


class Cube(Shape):
    """
    Class for cubes.
    """
    def __new__(self):
        return OgPoset.__new__(Point)

    def cube_face(self, k, sign):
        """ Cubical face maps. """
        sign = utils.mksign(sign)
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim),
            'why': 'out of bounds'})
        basic_faces = {
            '-': ShapeMap(OgMap(
                Point(), Arrow(),
                [[El(0, 0)]],
                wfcheck=False), wfcheck=False),
            '+': ShapeMap(OgMap(
                Point(), Arrow(),
                [[El(0, 1)]],
                wfcheck=False), wfcheck=False)}
        arrowid = Arrow().id()
        maps = [
                *[arrowid for _ in range(k)],
                basic_faces[sign],
                *[arrowid for _ in range(self.dim - k - 1)]
            ]
        return ShapeMap.gray(*maps)

    def cube_degeneracy(self, k):
        """ Cubical degeneracy maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim + 1),
            'why': 'out of bounds'})
        arrowid = Arrow().id()
        maps = [
                *[arrowid for _ in range(k)],
                Arrow().terminal(),
                *[arrowid for _ in range(self.dim - k)]
            ]
        return ShapeMap.gray(*maps)

    def cube_connection(self, k, sign):
        """ Cubical connection maps. """
        sign = utils.mksign(sign)
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim),
            'why': 'out of bounds'})
        basic_connections = {
                '-': ShapeMap(OgMap(
                        Shape.cube(2), Arrow(),
                        [[El(0, 0), El(0, 0), El(0, 1), El(0, 0)],
                         [El(0, 0), El(1, 0), El(0, 0), El(1, 0)],
                         [El(1, 0)]],
                        wfcheck=False), wfcheck=False),
                '+': ShapeMap(OgMap(
                        Shape.cube(2), Arrow(),
                        [[El(0, 0), El(0, 1), El(0, 1), El(0, 1)],
                         [El(1, 0), El(0, 1), El(1, 0), El(0, 1)],
                         [El(1, 0)]],
                        wfcheck=False), wfcheck=False)}
        arrowid = Arrow().id()
        maps = [
                *[arrowid for _ in range(k)],
                basic_connections[sign],
                *[arrowid for _ in range(self.dim - k - 1)]
            ]
        return ShapeMap.gray(*maps)


class Theta(Shape):
    """
    Class for the objects of the Theta category (Batanin cells).
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class OpetopeTree(Shape):
    """
    Class for shapes that can appear as input boundaries of positive
    opetopes.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class GlobeString(Theta, OpetopeTree):
    """
    Class for 'strings of globes' pasted in the top dimension, the
    intersection of opetope trees and Batanin cells.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Opetope(OpetopeTree):
    """
    Class for positive opetopes.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Globe(GlobeString, Opetope):
    """ Class for globes. """
    def __new__(self):
        return OgPoset.__new__(Point)


class Point(Globe, Simplex, Cube):
    """ Class for the point. """
    def __new__(self):
        return OgPoset.__new__(Point)

    def __init__(self, face_data=None, coface_data=None,
                 wfcheck=False, matchcheck=False):
        super().__init__(
                [[{'-': set(), '+': set()}]],
                [[{'-': set(), '+': set()}]],
                wfcheck=False, matchcheck=False)


class Arrow(Globe, Simplex, Cube):
    """ Class for the arrow. """
    def __new__(self):
        return OgPoset.__new__(Arrow)

    def __init__(self, face_data=None, coface_data=None,
                 wfcheck=False, matchcheck=False):
        super().__init__(
                [
                    [{'-': set(), '+': set()}, {'-': set(), '+': set()}],
                    [{'-': {0}, '+': {1}}]
                ], [
                    [{'-': {0}, '+': set()}, {'-': set(), '+': {0}}],
                    [{'-': set(), '+': set()}]
                ],
                wfcheck=False, matchcheck=False)


class ShapeMap(OgMap):
    """
    An OgMap overlay for total maps between Shape objects.
    """

    def __init__(self, ogmap, **params):
        wfcheck = params.get('wfcheck', True)
        if wfcheck:
            utils.typecheck(ogmap, {
                'type': OgMap,
                'st': lambda f: f.istotal,
                'why': 'a ShapeMap must be total'})
            for x in ogmap.source, ogmap.target:
                utils.typecheck(x, {'type': Shape})

        super().__init__(ogmap.source, ogmap.target, ogmap.mapping,
                         wfcheck=False)

    def then(self, other, *others):
        for f in [other, *others]:
            if not isinstance(f, ShapeMap):
                return super().then(other, *others)
        return ShapeMap(
                super().then(other, *others),
                wfcheck=False)

    @property
    def layers(self):
        if not hasattr(self.source, '_layering'):
            return [self]
        return [
                f.then(self)
                for f in self.source._layering
            ]

    @property
    def rewrite_steps(self):
        rewrite_steps = [
                *[layer.input for layer in self.layers],
                self.layers[-1].output
                ]
        return rewrite_steps

    @staticmethod
    def gray(*maps):
        for f in maps:
            if not isinstance(f, ShapeMap):
                return OgMap.gray(*maps)
        if len(maps) == 0:
            return Point().id()
        if len(maps) == 1:
            return maps[0]

        gray = OgMap.gray(*maps)
        if gray.source in [f.source for f in maps]:
            if gray.target in [f.target for f in maps]:
                return gray

        if gray.source not in [f.source for f in maps]:
            gray = Shape._reorder(gray.source).then(gray)
        if gray.target not in [f.target for f in maps]:
            gray = gray.then(Shape._reorder(gray.target).inv())

        def inheritance(l):
            if all([isinstance(x, Cube) for x in l]):
                return Cube
            return Shape

        return inheritance([f.source for f in maps])._upgrademapsrc(
                inheritance([f.target for f in maps])._upgrademaptgt(
                    gray))

    @staticmethod
    def join(*maps):
        for f in maps:
            if not isinstance(f, ShapeMap):
                return OgMap.join(*maps)
        if len(maps) == 0:
            return Empty().id()
        if len(maps) == 1:
            return maps[0]

        join = OgMap.join(*maps)
        if join.source in [f.source for f in maps]:
            if join.target in [f.target for f in maps]:
                return join

        if join.source not in [f.source for f in maps]:
            join = Shape._reorder(join.source).then(join)
        if join.target not in [f.target for f in maps]:
            join = join.then(Shape._reorder(join.target).inv())

        def inheritance(l):
            if all([isinstance(x, Simplex) for x in l]):
                if sum([x.dim+1 for x in l]) == 2:
                    return Arrow
                return Simplex
            return Shape

        return inheritance([f.source for f in maps])._upgrademapsrc(
                inheritance([f.target for f in maps])._upgrademaptgt(
                    join))

    def dual(shapemap, *dims):
        if not isinstance(shapemap, ShapeMap):
            return OgMap.dual(shapemap, *dims)

        ogdual = OgMap.dual(shapemap, *dims)
        dual = Shape._reorder(ogdual.source).then(ogdual).then(
                Shape._reorder(ogdual.target).inv())

        def inheritance(x, y):
            if x == y:
                return x.__class__
            if isinstance(x, Theta):
                return Theta
            return Shape

        return inheritance(shapemap.source, dual.source)._upgrademapsrc(
                inheritance(shapemap.target, dual.target)._upgrademaptgt(
                    dual))

    def generate_layering(self):
        self.source.generate_layering()

    def draw(self, **params):
        from rewal.strdiags import draw
        return draw(self, **params)

    def draw_boundaries(self, **params):
        from rewal.strdiags import draw_boundaries
        return draw_boundaries(self, **params)
