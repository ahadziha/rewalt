"""
Implements shapes of cells and diagrams.
"""

import networkx as nx

from rewalt import utils
from rewalt.ogposets import (El, OgPoset, GrSet, GrSubset, Closed,
                             OgMap, OgMapPair)


class Shape(OgPoset):
    """
    Inductive subclass of :class:`ogposets.OgPoset` for shapes of
    cells and diagrams.

    Properly formed objects of the class are unique encodings of the
    *regular molecules* from the theory of diagrammatic sets (plus the
    empty shape, which is not considered a regular molecule).

    To create shapes, we start from basic constructors such as
    :meth:`empty`, :meth:`point`, or one of the named shape constructors,
    such as :meth:`globe`, :meth:`simplex`, :meth:`cube`.

    Then we generate new shapes by gluing basic shapes together with
    :meth:`paste`, :meth:`to_inputs`, :meth:`to_outputs`, or by
    producing new higher-dimensional shapes with operations such as
    :meth:`atom`, :meth:`gray`, :meth:`join`.

    When possible, the constructors place the shapes in appropriate
    subclasses of separate interest, which include the *globes*,
    the *oriented simplices*, the *oriented cubes*, and the
    *positive opetopes*. This is to enable the specification of special
    methods for subclasses of shapes.

    The following diagram summarises the hierarchy of subclasses of
    shapes:
    ::

         Simplex  Cube  OpetopeTree  Theta
         |    |\  |\     |      |      |
         |    | \ | \ Opetope  GlobeString
         |    |  \|  \   |    /
         |    |   \   \  Globe
         |    |   |\   \/    |
        Empty |   | \  /\    |
              |   |  \/  \   |
              |   |  /\   \  |
              |   | /  \   \ |
              |   |/    \   \|
              Point      Arrow

    Currently only the :class:`Cube` and :class:`Simplex` classes have
    special methods implemented.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    @property
    def isatom(self):
        """
        Returns whether the shape is an atom (has a greatest element).

        Returns
        -------
        isatom : :class:`bool`
            :code:`True` if and only if the shape has a greatest element.

        Examples
        --------
        >>> arrow = Shape.arrow()
        >>> assert arrow.isatom
        >>> assert not arrow.paste(arrow).isatom
        """
        return len(self.maximal()) == 1

    @property
    def isround(self):
        """
        Shorthand for :code:`all().isround`.
        """
        return self.all().isround

    @property
    def layers(self):
        """
        Returns the current layering of the shape.

        Returns
        -------
        layers : :class:`list[ShapeMap]`
            The current layering.

        Examples
        --------
        >>> arrow = Shape.arrow()
        >>> globe = Shape.globe(2)
        >>> cospan = globe.paste(arrow).paste(
        ...         arrow.paste(globe), cospan=True)
        >>> shape = cospan.target
        >>> assert shape.layers == [cospan.fst, cospan.snd]
        """
        return self.id().layers

    @property
    def rewrite_steps(self):
        """
        Returns the sequence of rewrite steps associated to the current
        layering of the shape.

        The :code:`0`-th rewrite step is the input boundary of the shape.
        For :code:`n > 0`, the :code:`n`-th rewrite step is the output
        boundary of the :code:`(n-1)`-th layer.

        Returns
        -------
        rewrite_steps : :class:`list[ShapeMap]`
            The current sequence of rewrite steps.

        Examples
        --------
        >>> arrow = Shape.arrow()
        >>> globe = Shape.globe(2)
        >>> cospan = globe.paste(arrow).paste(
        ...         arrow.paste(globe), cospan=True)
        >>> shape = cospan.target
        >>> assert shape.rewrite_steps == [
        ...         cospan.fst.input,
        ...         cospan.fst.output,
        ...         cospan.snd.output]
        """
        return self.id().rewrite_steps

    # Main constructors
    @staticmethod
    def atom(fst, snd, **params):
        """
        Given two shapes with identical round boundaries, returns a new
        atomic shape whose input boundary is the first one and output
        boundary the second one.

        Arguments
        ---------
        fst : :class:`Shape`
            The input boundary shape.
        snd : :class:`Shape`
            The output boundary shape.

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the input and
            output boundaries (default is :code:`False`).

        Returns
        -------
        atom : :class:`Shape` | :class:`ogposets.OgMapPair`
            The new atomic shape (optionally with the cospan of
            inclusions of its boundaries).

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match, or are not round.

        Examples
        --------
        We create a 2-dimensional cell shape with two input 1-cells
        and one output 2-cell.

        >>> arrow = Shape.arrow()
        >>> binary = arrow.paste(arrow).atom(arrow)
        >>> binary.draw(path='docs/_static/img/Shape_atom.png')

        .. image:: ../_static/img/Shape_atom.png
            :width: 400
            :align: center
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
        Given two shapes and :code:`k` such that the output
        :code:`k`-boundary of the first is equal to the input
        :code:`k`-boundary of the second, returns their pasting along
        the matching boundaries.

        Arguments
        ---------
        fst : :class:`Shape`
            The first shape.
        snd : :class:`Shape`
            The second shape.
        dim : :class:`int`, optional
            The dimension of the boundary along which they will be pasted
            (default is :code:`min(fst.dim, snd.dim) - 1`).

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the two shapes
            into the pasting (default is :code:`False`).

        Returns
        -------
        paste : :class:`Shape` | :class:`ogposets.OgMapPair`
            The pasted shape (optionally with the cospan of
            inclusions of its components).

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match.

        Examples
        --------
        We can paste two 2-dimensional globes either "vertically" along
        their 1-dimensional boundary or "horizontally" along their
        0-dimensional boundary.

        >>> globe = Shape.globe(2)
        >>> vert = globe.paste(globe)
        >>> horiz = globe.paste(globe, 0)
        >>> vert.draw(path='docs/_static/img/Shape_paste_vert.png')

        .. image:: ../_static/img/Shape_paste_vert.png
            :width: 400
            :align: center

        >>> horiz.draw(path='docs/_static/img/Shape_paste_horiz.png')

        .. image:: ../_static/img/Shape_paste_horiz.png
            :width: 400
            :align: center

        We can also check that the interchange equation holds.

        >>> assert vert.paste(vert, 0) == horiz.paste(horiz)
        >>> horiz.paste(horiz).draw(
        ...     path='docs/_static/img/Shape_paste_interchange.png')

        .. image:: ../_static/img/Shape_paste_interchange.png
            :width: 400
            :align: center
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
        Given a span of shape maps, where one is the inclusion of the
        input (resp output) :code:`k`-boundary of a shape,
        and the other the inclusion of a round subshape of the
        output (resp input) :code:`k`-boundary of another shape,
        returns the pasting (pushout) of the two shapes along the span.

        In practice, it is convenient to use :meth:`to_inputs` and
        :meth:`to_outputs` instead, where the data of the span is specified
        by :code:`k` and the positions of the :code:`k`-dimensional
        elements in the round subshape along which the pasting occurs.

        Arguments
        ---------
        fst : :class:`ShapeMap`
            The first inclusion.
        snd : :class:`ShapeMap`
            The second inclusion.

        Keyword arguments
        -----------------
        wfcheck : :class:`bool`
            Check if the span gives rise to a well-formed pasting
            (default is :code:`True`).
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the two shapes
            into the pasting (default is :code:`False`).

        Returns
        -------
        paste_along : :class:`Shape` | :class:`ogposets.OgMapPair`
            The pasted shape (optionally with the cospan of
            inclusions of its components).

        Raises
        ------
        :class:`ValueError`
            If the pair of maps is not an injective span.
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
        Returns the pasting of another shape along a round subshape of
        the output :code:`k`-boundary, specified by the positions of its
        :code:`k`-dimensional elements.

        Arguments
        ---------
        positions : :class:`list[int]` | :class:`int`
            The positions of the outputs along which to paste. If given
            an integer :code:`n`, interprets it as the list :code:`[n]`.
        other : :class:`Shape`
            The other shape to paste.
        dim : :class:`int`, optional
            The dimension of the boundary along which to paste
            (default is :code:`self.dim - 1`)

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the two shapes
            into the pasting (default is :code:`False`).

        Returns
        -------
        to_outputs : :class:`Shape` | :class:`ogposets.OgMapPair`
            The pasted shape (optionally with the cospan of
            inclusions of its components).

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match, or the pasting does not produce
            a well-formed shape.

        Examples
        --------
        We create a 2-simplex and visualise it as a string diagram with the
        :code:`positions` parameter enabled.

        >>> simplex = Shape.simplex(2)
        >>> simplex.draw(
        ...     positions=True, path='docs/_static/img/Shape_to_outputs1.png')

        .. image:: ../_static/img/Shape_to_outputs1.png
            :width: 400
            :align: center

        We paste another 2-simplex to the output in position :code:`2`.

        >>> paste1 = simplex.to_outputs(2, simplex)
        >>> paste1.draw(
        ...     positions=True, path='docs/_static/img/Shape_to_outputs2.png')

        .. image:: ../_static/img/Shape_to_outputs2.png
            :width: 400
            :align: center

        Finally, we paste the *dual* of a 2-simplex to the outputs in
        positions :code:`2, 3`.

        >>> paste2 = paste1.to_outputs([1, 3], simplex.dual())
        >>> paste2.draw(
        ...     positions=True, path='docs/_static/img/Shape_to_outputs3.png')

        .. image:: ../_static/img/Shape_to_outputs3.png
            :width: 400
            :align: center
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
        Returns the pasting of another shape along a round subshape
        of the input :code:`k`-boundary, specified by the positions of its
        :code:`k`-dimensional elements.

        Arguments
        ---------
        positions : :class:`list[int]` | :class:`int`
            The positions of the inputs along which to paste. If given
            an integer :code:`n`, interprets it as the list :code:`[n]`.
        other : :class:`Shape`
            The other shape to paste.
        dim : :class:`int`, optional
            The dimension of the boundary along which to paste
            (default is :code:`self.dim - 1`)

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the two shapes
            into the pasting (default is :code:`False`).

        Returns
        -------
        to_inputs : :class:`Shape` | :class:`ogposets.OgMapPair`
            The pasted shape (optionally with the cospan of
            inclusions of its components).

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match, or the pasting does not produce
            a well-formed shape.

        Examples
        --------
        We work dually to the example for :meth:`to_outputs`.

        >>> binary = Shape.simplex(2).dual()
        >>> binary.draw(
        ...     positions=True, path='docs/_static/img/Shape_to_inputs1.png')

        .. image:: ../_static/img/Shape_to_inputs1.png
            :width: 400
            :align: center

        >>> paste1 = binary.to_inputs(1, binary)
        >>> paste1.draw(
        ...     positions=True, path='docs/_static/img/Shape_to_inputs2.png')

        .. image:: ../_static/img/Shape_to_inputs2.png
            :width: 400
            :align: center

        >>> paste2 = paste1.to_inputs([0, 1], binary.dual())
        >>> paste2.draw(
        ...     positions=True, path='docs/_static/img/Shape_to_inputs3.png')

        .. image:: ../_static/img/Shape_to_inputs3.png
            :width: 400
            :align: center
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
        """
        Returns the n-fold suspension of a shape.

        This static method can be also used as a bound method after
        an object is initialised, that is, :code:`shape.suspend(n)` is
        equivalent to :code:`suspend(shape, n)`.

        Arguments
        ---------
        shape : :class:`Shape`
            The object to suspend.
        n : :class:`int`, optional
            The number of iterations of the suspension (default is 1).

        Returns
        -------
        suspension : :class:`Shape`
            The suspended shape.

        Examples
        --------
        The suspension of the point is the arrow, and the suspension of
        an arrow is the 2-globe.

        >>> assert Shape.point().suspend() == Shape.arrow()
        >>> assert Shape.arrow().suspend() == Shape.globe(2)

        In general, the suspension of the n-globe is the (n+1)-globe.
        """
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
        """
        Returns the Gray product of any number of shapes.

        This method can be called with the math operator :code:`*`, that is,
        :code:`fst * snd` is equivalent to :code:`gray(fst, snd)`.

        This static method can also be used as a bound method after an object
        is initialised, that is, :code:`fst.gray(*shapes)` is equivalent to
        :code:`gray(fst, *shapes)`.

        Arguments
        ---------
        *shapes : :class:`Shape`
            Any number of shapes.

        Returns
        -------
        gray : :class:`Shape`
            The Gray product of the arguments.

        Example
        -------
        The point is a unit for the Gray product.

        >>> point = Shape.point()
        >>> arrow = Shape.arrow()
        >>> assert point*arrow == arrow*point == arrow

        The Gray product of two arrows is the oriented square (2-cube).

        >>> arrow = Shape.arrow()
        >>> assert arrow*arrow == Shape.cube(2)

        In general, the Gray product of the n-cube with the k-cube
        is the (n+k)-cube.
        """
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
        """
        Returns the join of any number of shapes.

        This method can be called with the shift operators :code:`>>`
        and :code:`<<`, that is, :code:`fst >> snd` is equivalent to
        :code:`join(fst, snd)` and :code:`fst << snd` is equivalent to
        :code:`join(snd, fst)`.

        This static method can also be used as a bound method after an
        object is initialised, that is, :code:`fst.join(*shapes)` is
        equivalent to :code:`join(fst, *shapes)`.

        Arguments
        ---------
        *shapes : :class:`Shape`
            Any number of shapes.

        Returns
        -------
        join : :class:`Shape`
            The join of the arguments.

        Examples
        --------
        The empty shape is a unit for the join.

        >>> empty = Shape.empty()
        >>> point = Shape.point()
        >>> assert empty >> point == point >> empty == point

        The join of two points is the arrow, and the join of an arrow
        and a point is the 2-simplex.

        >>> arrow = Shape.arrow()
        >>> assert point >> point == Shape.arrow()
        >>> assert arrow >> point == Shape.simplex(2)

        In general, the join of an n-simplex with a k-simplex is
        the (n+k+1)-simplex.
        """
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
        """
        Returns the shape with orientations reversed in given dimensions.

        The dual in all dimensions can also be called with the bit negation
        operator :code:`~`, that is, :code:`~shape` is equivalent to
        :code:`shape.dual()`.

        This static method can be also used as a bound method after an object
        is initialised, that is, :code:`shape.dual(*dims)` is equivalent to
        :code:`dual(shape, *dims)`.

        Arguments
        ---------
        shape : :class:`Shape`
            A shape.
        *dims : :class:`int`
            Any number of dimensions; if none, defaults to *all* dimensions.

        Returns
        -------
        dual : :class:`Shape`
            The shape, dualised in the given dimensions.

        Examples
        --------
        >>> arrow = Shape.arrow()
        >>> simplex = Shape.simplex(2)
        >>> binary = arrow.paste(arrow).atom(arrow)
        >>> assert binary == simplex.dual()

        >>> assoc_l = binary.to_inputs(0, binary)
        >>> assoc_r = binary.to_inputs(1, binary)
        >>> assert assoc_r == assoc_l.dual(1)
        """
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
        Returns the unique atomic shape with the same boundary,
        if the shape is round.

        Returns
        -------
        merge : :class:`Shape`
            The unique atomic shape with the same boundary.

        Raises
        ------
        :class:`ValueError`
            If the shape is not round.

        Examples
        --------
        We create a 2-dimensional shape with two input 1-cells and
        one output 1-cell, and paste it to itself along one of the
        inputs.

        >>> arrow = Shape.arrow()
        >>> binary = arrow.paste(arrow).atom(arrow)
        >>> to_merge = binary.to_inputs(1, binary)
        >>> to_merge.draw(path='docs/_static/img/Shape_merge1.png')

        .. image:: ../_static/img/Shape_merge1.png
            :width: 400
            :align: center

        The "merged" shape is the 2-dimensional atom with three input
        2-cells and one output 1-cell.

        >>> merged = to_merge.merge()
        >>> merged.draw(path='docs/_static/img/Shape_merge2.png')

        .. image:: ../_static/img/Shape_merge2.png
            :width: 400
            :align: center
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
        """
        Constructs the initial, empty shape.

        Returns
        -------
        empty : :class:`Empty`
            The empty shape.
        """
        return Empty()

    @staticmethod
    def point():
        """
        Constructs the terminal shape, consisting of a single point.

        Returns
        -------
        point : :class:`Point`
            The point.
        """
        return Point()

    @staticmethod
    def arrow():
        """
        Constructs the arrow, the unique 1-dimensional atomic shape.

        Returns
        -------
        arrow : :class:`Arrow`
            The arrow.
        """
        return Arrow()

    @staticmethod
    def simplex(dim=-1):
        """
        Constructs the oriented simplex of a given dimension.

        Arguments
        ---------
        dim : :class:`int`
            The dimension of the simplex (default is :code:`-1`).

        Returns
        -------
        simplex : :class:`Simplex`
            The simplex of the requested dimension.
        """
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= -1,
            'why': 'expecting integer >= -1'})
        point = Point()
        return Shape.join(*[point for _ in range(dim+1)])

    @staticmethod
    def cube(dim=0):
        """
        Constructs the oriented cube of a given dimension.

        Arguments
        ---------
        dim : :class:`int`
            The dimension of the cube (default is :code:`0`).

        Returns
        -------
        cube : :class:`Cube`
            The cube of the requested dimension.
        """
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})
        arrow = Arrow()
        return Shape.gray(*[arrow for _ in range(dim)])

    @staticmethod
    def globe(dim=0):
        """
        Constructs the globe of a given dimension.

        Arguments
        ---------
        dim : :class:`int`
            The dimension of the globe (default is :code:`0`).

        Returns
        -------
        globe : :class:`Globe`
            The globe of the requested dimension.
        """
        return Shape.suspend(Point(), dim)

    @staticmethod
    def theta(*thetas):
        """
        Inductive constructor for the objects of the Theta category,
        sometimes known as Batanin cells.

        Batanin cells are in 1-to-1 correspondence with finite plane trees.
        The constructor is based on this correspondence, using the
        well-known inductive definition of plane trees: given any number
        :code:`k` of Batanin cells, it returns the Batanin cell encoded by
        a root with :code:`k` children, to which the :code:`k` plane trees
        encoding the arguments are attached.

        Arguments
        ---------
        thetas : :class:`Theta`
            Any number of Batanin cells.

        Returns
        -------
        theta : :class:`Theta`
            The resulting Batanin cell.

        Examples
        --------
        Every globe is a Batanin cell, encoded by the linear tree of length
        equal to its dimension.

        >>> assert Shape.theta() == Shape.globe(0)
        >>> assert Shape.theta(Shape.theta()) == Shape.globe(1)
        >>> assert Shape.theta(Shape.theta(Shape.theta())) == Shape.globe(2)

        The tree with one root with n children corresponds to a string
        of n arrows.

        >>> point = Shape.theta()
        >>> arrow = Shape.arrow()
        >>> assert Shape.theta(point, point) == arrow.paste(arrow)
        """
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
        """
        Returns the identity map on the shape.

        Returns
        -------
        id : :class:`ShapeMap`
            The identity map on the object.
        """
        return ShapeMap(super().id(),
                        wfcheck=False)

    def boundary(self, sign=None, dim=None):
        """
        Returns the inclusion of the boundary of a given orientation
        and dimension into the shape.

        Note that input and output boundaries of shapes are shapes,
        so they are returned as shape maps; however, the entire (input
        + output) boundary of a shape is not a shape, so it is returned
        simply as a map of oriented graded posets.

        Arguments
        ---------
        sign : :class:`str`, optional
            Orientation: :code:`'-'` for input, :code:`'+'` for output,
            :code:`None` (default) for both.
        dim : :class:`int`, optional
            Dimension of the boundary (default is :code:`self.dim - 1`).

        Returns
        -------
        boundary : :class:`ShapeMap` | :class:`OgMap`
            The inclusion of the requested boundary into the object.

        Examples
        --------
        >>> point = Shape.point()
        >>> arrow = Shape.arrow()
        >>> binary = arrow.paste(arrow).atom(arrow)
        >>> assert binary.boundary('-').source == arrow.paste(arrow)
        >>> assert binary.boundary('+').source == arrow
        >>> assert binary.boundary('-', 0).source == point
        >>> assert binary.boundary('-').target == binary
        """
        dim = self.dim-1 if dim is None else dim
        if dim >= self.dim:
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
                    if dim == 1:
                        return Arrow
                    return Opetope
                return OpetopeTree
            return Shape

        return inheritance()._upgrademapsrc(boundary)

    def atom_inclusion(self, element):
        """
        Returns the inclusion of the closure of an element, which
        is an atomic shape, in the shape.

        Arguments
        ---------
        element : :class:`El`
            An element of the shape.

        Returns
        -------
        atom_inclusion : :class:`ShapeMap`
            The inclusion of the closure of the element.

        Examples
        --------
        >>> arrow = Shape.arrow()
        >>> globe = Shape.globe(2)
        >>> whisker_l = arrow.paste(globe)
        >>> assert whisker_l.atom_inclusion(El(2, 0)).source == globe
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
        Returns the unique map from the initial, empty shape.

        Returns
        -------
        initial : :class:`ShapeMap`
            The unique map from the empty shape.

        Examples
        --------
        >>> point = Shape.point()
        >>> empty = Shape.empty()
        >>> assert point.initial() == empty.terminal()
        >>> assert empty.initial() == empty.id()
        """
        return ShapeMap(
                OgMap(Empty(), self,
                      wfcheck=False),
                wfcheck=False)

    def terminal(self):
        """
        Returns the unique map to the point, the terminal shape.

        Returns
        -------
        terminal : :class:`ShapeMap`
            The unique map to the point.

        Examples
        --------
        >>> point = Shape.point()
        >>> assert point.terminal() == point.id()
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
        Given a closed subset of the boundary of the shape, forms a
        cylinder on the shape, with the sides incident to the closed subset
        collapsed, and returns its projection map onto the original shape.

        This is mainly used in constructing units and unitors on diagrams;
        see :meth:`diagrams.Diagram.unit`, :meth:`diagrams.Diagram.lunitor`,
        :meth:`diagrams.Diagram.runitor`.

        Arguments
        ---------
        collapsed : :class:`Closed`, optional
            A closed subset of the boundary of the shape (default is
            the entire boundary).

        Returns
        -------
        inflate : :class:`Closed`
            The projection map of the "partially collapsed cylinder" onto
            the shape.

        Raises
        ------
        :class:`ValueError`
            If `collapsed` is not a subset of the boundary.
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
        Returns an iterator on all *layerings* of a shape of dimension
        :code:`n` into shapes with a single :code:`n`-dimensional element,
        pasted along their :code:`(n-1)`-dimensional boundary.

        Returns
        -------
        all_layerings : :class:`Iterable`
            The iterator on all layerings of the shape.
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
        Assigns a layering to the shape, iterating through all
        the layerings, and returns it.

        Returns
        -------
        layers : :class:`list[ShapeMap]`
            The generated layering.

        Examples
        --------
        >>> arrow = Shape.arrow()
        >>> globe = Shape.globe(2)
        >>> chain = globe.paste(globe, 0)
        >>> chain.generate_layering()
        >>> assert chain.layers[0].source == arrow.paste(globe)
        >>> assert chain.layers[1].source == globe.paste(arrow)
        >>> chain.generate_layering()
        >>> assert chain.layers[0].source == globe.paste(arrow)
        >>> assert chain.layers[1].source == arrow.paste(globe)
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
        """
        Bound version of :meth:`strdiags.draw`.

        Calling :code:`x.draw(**params)` is equivalent to calling
        :code:`strdiags.draw(x, **params)`.
        """
        from rewalt.strdiags import draw
        return draw(self, **params)

    def draw_boundaries(self, **params):
        """
        Bound version of :meth:`strdiags.draw_boundaries`.

        Calling :code:`x.draw_boundaries(**params)` is equivalent to
        calling :code:`strdiags.draw_boundaries(x, **params)`.
        """
        from rewalt.strdiags import draw_boundaries
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
            dim = focus.dim
            top_dim = focus[dim]
            top_unmarked = GrSet()
            for x in top_dim:
                if x not in marked:
                    top_unmarked.add(x)
            if len(top_unmarked) == 0:
                del focus_stack[-1]
            else:
                if dim == 0:
                    for x in top_dim:
                        mapping[dim].append(x)
                        marked.add(x)
                    del focus_stack[-1]
                else:
                    focus_in_faces = GrSet()
                    focus_out_faces = GrSet()
                    candidates = {}
                    for x in top_dim:
                        isunmarked = x in top_unmarked
                        for y in shape.faces(x, '-'):
                            focus_in_faces.add(y)
                            if isunmarked:
                                candidates[y] = x
                        for y in shape.faces(x, '+'):
                            focus_out_faces.add(y)
                    focus_input = focus_in_faces.difference(
                            focus_out_faces)
                    if not focus_input.issubset(marked):
                        if len(focus[:dim]) > 0:
                            focus_input = focus_input.union(
                                focus[:dim])
                        focus_stack.append(focus_input)
                    else:
                        if len(focus) == 1:
                            for x in top_dim:
                                mapping[dim].append(x)
                                marked.add(x)
                            del focus_stack[-1]
                            focus_output = focus_out_faces.difference(
                                    focus_in_faces)
                            if not focus_output.issubset(marked):
                                if len(focus[:dim]) > 0:
                                    focus_output = focus_output.union(
                                            focus[:dim])
                                focus_stack.append(focus_output)
                        else:
                            x = next(
                                    x for x in mapping[dim-1]
                                    if x in candidates.keys())
                            focus_stack.append(
                                GrSet(candidates[x]))

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
    Subclass of :class:`Shape` for oriented simplices.

    The methods of this class provide a full implementation of the
    category of simplices, which is generated by the face and
    degeneracy maps between simplices one dimension apart.

    Use :meth:`Shape.simplex` to construct.

    Examples
    --------
    We create a 1-simplex (arrow), a 2-simplex (triangle),
    and a 3-simplex (tetrahedron).

    >>> arrow = Shape.simplex(1)
    >>> triangle = Shape.simplex(2)
    >>> tetra = Shape.simplex(3)

    We can then check some of the simplicial relations between
    degeneracy and face maps.

    >>> map1 = triangle.simplex_degeneracy(2).then(
    ...     arrow.simplex_degeneracy(1))
    >>> map2 = triangle.simplex_degeneracy(1).then(
    ...     arrow.simplex_degeneracy(1))
    >>> assert map1 == map2

    >>> map3 = tetra.simplex_face(2).then(
    ...     triangle.simplex_degeneracy(2))
    >>> assert map3 == triangle.id()

    >>> map4 = tetra.simplex_face(0).then(
    ...     triangle.simplex_degeneracy(2))
    >>> map5 = arrow.simplex_degeneracy(1).then(
    ...     triangle.simplex_face(0))
    >>> assert map4 == map5
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    def simplex_face(self, k):
        """
        Returns one of the face inclusion maps of the simplex.

        Arguments
        ---------
        k : :class:`int`
            The index of the face map, ranging from :code:`0` to
            :code:`self.dim`.

        Returns
        -------
        simplex_face : :class:`ShapeMap`
            The face map.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
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
        """
        Returns one of the collapse (degeneracy) maps of the simplex
        one dimension higher.

        Arguments
        ---------
        k : :class:`int`
            The index of the degeneracy map, ranging from :code:`0` to
            :code:`self.dim`.

        Returns
        -------
        simplex_degeneracy : :class:`ShapeMap`
            The degeneracy map.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
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
    Subclass of :class:`Shape` for the empty shape.

    Use :meth:`Shape.empty` to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    def __init__(self, face_data=None, coface_data=None,
                 **params):
        super().__init__([], [],
                         wfcheck=False, matchcheck=False)


class Cube(Shape):
    """
    Subclass of :class:`Shape` for oriented cubes.

    The methods of this class provide a full implementation of the
    category of cubes with connections, which is generated by the
    face, degeneracy, and connection maps between cubes one
    dimension apart.

    Use :meth:`Shape.cube` to construct.

    Examples
    --------
    We create a 1-cube (arrow), 2-cube (square), and 3-cube (cube).

    >>> arrow = Shape.cube(1)
    >>> square = Shape.cube(2)
    >>> cube = Shape.cube(3)

    We can then check some of the relations between cubical face,
    connection, and degeneracy maps.

    >>> map1 = square.cube_degeneracy(2).then(
    ...     arrow.cube_degeneracy(1))
    >>> map2 = square.cube_degeneracy(1).then(
    ...     arrow.cube_degeneracy(1))
    >>> assert map1 == map2

    >>> map3 = square.cube_face(0, '+').then(
    ...     cube.cube_face(2, '-'))
    >>> map4 = square.cube_face(1, '-').then(
    ...     cube.cube_face(0, '+'))
    >>> assert map3 == map4

    >>> map5 = square.cube_connection(1, '-').then(
    ...     arrow.cube_connection(0, '-'))
    >>> map6 = square.cube_connection(0, '-').then(
    ...     arrow.cube_connection(0, '-'))
    >>> assert map5 == map6
    """
    def __new__(self):
        return OgPoset.__new__(Point)

    def cube_face(self, k, sign):
        """
        Returns one of the face inclusion maps of the cube.

        Arguments
        ---------
        k : :class:`int`
            Index of the face map, ranging from :code:`0` to
            :code:`self.dim - 1`.
        sign : :class:`str`
            Side: :code:`'-'` or :code:`'+'`.

        Returns
        -------
        cube_face : :class:`ShapeMap`
            The face map.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
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
        """
        Returns one of the "degeneracy" collapse maps of the cube
        one dimension higher.

        Arguments
        ---------
        k : :class:`int`
            The index of the degeneracy map, ranging from :code:`0` to
            :code:`self.dim`.

        Returns
        -------
        cube_degeneracy : :class:`ShapeMap`
            The degeneracy map.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
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
        """
        Returns one of the "connection" collapse maps of the cube
        one dimension higher.

        Arguments
        ---------
        k : :class:`int`
            Index of the connection map, ranging from :code:`0` to
            :code:`self.dim - 1`.
        sign : :class:`str`
            Side: :code:`'-'` or :code:`'+'`.

        Returns
        -------
        cube_face : :class:`ShapeMap`
            The connection map.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
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
    Subclass of :class:`Shape` for Batanin cells.

    Use :meth:`Shape.theta` to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class OpetopeTree(Shape):
    """
    Subclass of :class:`Shape` for shapes that can appear as input
    boundaries of opetopes.

    Use :class:`Shape` methods to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class GlobeString(Theta, OpetopeTree):
    """
    Subclass of :class:`Shape` for "strings of globes" pasted in the top
    dimension.

    This is the "intersection" of :class:`OpetopeTree` and :class:`Theta`.
    Use :class:`Shape` methods to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Opetope(OpetopeTree):
    """
    Subclass of :class:`Shape` for (positive) opetopes.

    Use :class:`Shape` methods to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Globe(GlobeString, Opetope):
    """
    Subclass of :class:`Shape` for globes.

    Use :meth:`Shape.globe` to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Point(Globe, Simplex, Cube):
    """
    Subclass of :class:`Shape` for the point.

    Use :meth:`Shape.point` to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Point)

    def __init__(self, face_data=None, coface_data=None,
                 **params):
        super().__init__(
                [[{'-': set(), '+': set()}]],
                [[{'-': set(), '+': set()}]],
                wfcheck=False, matchcheck=False)


class Arrow(Globe, Simplex, Cube):
    """
    Subclass of :class:`Shape` for the arrow shape.

    Use :meth:`Shape.arrow` to construct.
    """
    def __new__(self):
        return OgPoset.__new__(Arrow)

    def __init__(self, face_data=None, coface_data=None,
                 **params):
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
    An overlay of :class:`ogposets.OgMap` for total maps between
    :class:`Shape` objects.

    It is used to extend constructions of shapes functorially to their
    maps, in a way that is compatible with the unique representation
    of shapes by their underlying :class:`ogposets.OgPoset` objects.

    The most common :class:`ShapeMap` objects are created by methods of
    :class:`Shape` such as :meth:`Shape.boundary` and :meth:`Shape.inflate`,
    or of its subclasses, such as :meth:`Simplex.simplex_degeneracy` or
    :meth:`Cube.cube_connection`.

    Nevertheless, occasionally we may need to define a map explicitly,
    in which case we first define an object :code:`f` of class
    :class:`ogposets.OgMap`, then upgrade it to a :class:`ShapeMap`
    with the constructor :code:`ShapeMap(f)`.

    Arguments
    ---------
    ogmap : :class:`ogposets.OgMap`
        A total map between shapes.

    Keyword arguments
    -----------------
    wfcheck : :class:`bool`
        Check whether the given map is a total map between shapes
        (default is :code:`True`).
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
        """
        Returns the current layering of the map's source, composed
        with the map.

        Returns
        -------
        layers : :class:`list[ShapeMap]`
            The source's current layering, composed with the map.
        """
        if not hasattr(self.source, '_layering'):
            return [self]
        return [
                f.then(self)
                for f in self.source._layering
            ]

    @property
    def rewrite_steps(self):
        """
        Returns the sequence of rewrite steps associated to the current
        layering of the map's source, composed with the map.

        Returns
        -------
        rewrite_steps : :class:`list[ShapeMap]`
            The source's current sequence of rewrite steps, composed
            with the map.
        """
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
            return Shape.point().id()
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
            return Shape.empty().id()
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
        """
        Shorthand for :code:`source.generate_layering()`.
        """
        self.source.generate_layering()

    def draw(self, **params):
        """
        Bound version of :meth:`strdiags.draw`.

        Calling :code:`f.draw(**params)` is equivalent to calling
        :code:`strdiags.draw(f, **params)`.
        """
        from rewalt.strdiags import draw
        return draw(self, **params)

    def draw_boundaries(self, **params):
        """
        Bound version of :meth:`strdiags.draw_boundaries`.

        Calling :code:`f.draw_boundaries(**params)` is equivalent to calling
        :code:`strdiags.draw_boundaries(f, **params)`.
        """
        from rewalt.strdiags import draw_boundaries
        return draw_boundaries(self, **params)
