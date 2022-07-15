"""
Implements diagrammatic sets and diagrams.
"""

from rewalt import utils, shapes
from rewalt.ogposets import (El, GrSet, GrSubset)
from rewalt.shapes import (Shape, ShapeMap)


class DiagSet:
    """
    Class for diagrammatic sets, a model of higher-dimensional rewrite
    systems and/or directed cell complexes.

    A diagrammatic set is constructed by creating an empty object, then
    adding named *generators* of different dimensions. The addition of a
    generator models the gluing of an atomic :class:`shapes.Shape` object
    along its boundary.

    This operation produces a *diagram*, that is, a map from a shape
    to the diagrammatic set, modelled as a :class:`Diagram` object.
    From these "basic" diagrams, we can construct "derived" diagrams
    either by pasting, or by pulling back along shape maps (this is
    used to produce "unit" or "degenerate" diagrams).

    To add a 0-dimensional generator (a point), we just give it a name.
    In the main constructor :meth:`add`, the gluing of an
    :code:`n`-dimensional generator is specified by a pair of round,
    :code:`(n-1)`-dimensional :class:`Diagram` objects, describing
    the gluing maps for the input and output boundaries of a shape.

    Simplicial sets, cubical sets with connections, and reflexive globular
    sets are all special cases of diagrammatic sets, where the generators
    have simplicial, cubical, or globular shapes.
    There are special constructors :meth:`add_simplex` and
    :meth:`add_cube` for adding simplicial and cubical generators by
    listing all their faces.

    The generators of a diagrammatic set are, by default, "directed" and
    not invertible. The class supports a model of weak or pseudo-
    invertibility, where two generators being each other's "weak inverse"
    is witnessed by a pair of higher-dimensional generators (*invertors*).
    This is produced by the methods :meth:`invert` (creates an inverse) and
    :meth:`make_inverses` (makes an existing generator the inverse).

    Diagrammatic sets do not have an intrinsic notion of composition
    of diagrams, so they are not by themselves a model of higher categories.
    However, the class supports a model of higher categories in which one
    generator being the composite of a diagram is witnessed by a
    higher-dimensional generator (a *compositor*). This is produced
    by the methods :meth:`compose` (creates a composite) and
    :meth:`make_composite` (makes an existing generator the composite).

    Notes
    -----
    There is an alternative constructor :meth:`yoneda` which turns
    a :class:`shapes.Shape` object into a diagrammatic set with one
    generator for every face of the shape.
    """
    _PRIVATE = (
            'shape',
            'mapping',
            'faces',
            'cofaces',
            'inverse',
            'linvertor',
            'rinvertor',
            'inverts',
            'composite')

    def __init__(self):
        self._generators = dict()
        self._by_dim = dict()
        self._compositors = dict()

    def __str__(self):
        return '{} with {} generators'.format(
                type(self).__name__, str(len(self)))

    def __getitem__(self, key):
        if key in self.generators:
            return Diagram._new(
                    self.generators[key]['shape'],
                    self,
                    self.generators[key]['mapping'],
                    key)
        raise KeyError(str(key))

    def __contains__(self, item):
        return item in self.generators.keys()

    def __len__(self):
        return len(self.generators)

    def __iter__(self):
        return iter(self.generators)

    @property
    def generators(self):
        """
        Returns the object's internal representation of the set of
        generators and related data.

        This is a dictionary whose keys are the generators' names.
        For each generator, the object stores another dictionary,
        which contains at least

        - the generator's shape (:code:`shape`, :class:`shapes.Shape`),
        - the mapping of the shape (:code:`mapping`,
          :class:`list[list[hashable]]`),
        - the generator's set of "faces", that is, other generators
          appearing as codimension-1 faces of the generator
          (:code:`faces`, :class:`set[hashable]`),
        - the generator's set of "cofaces", that is, other generators
          that have the generator as a face (:code:`cofaces`,
          :class:`set[hashable]`).

        If the generator has been inverted, it will also contain

        - its inverse's name (:code:`inverse`, :class:`hashable`),
        - the left invertor's name (:code:`linvertor`, :class:`hashable`),
        - the right invertor's name (:code:`rinvertor`, :class:`hashable`).

        If the generator happens to be a compositor, it will also
        contain the name of the composite it is exhibiting
        (:code:`composite`, :class:`hashable`).

        This also stores any additional keyword arguments passed when
        adding the generator.

        Returns
        -------
        generators : :class:`dict[dict]`
            The generators data.
        """
        return self._generators

    @property
    def by_dim(self):
        """
        Returns the set of generators indexed by dimension.

        Returns
        -------
        by_dim : :class:`dict[hashable]`
            The set of generators indexed by dimension.
        """
        return self._by_dim

    @property
    def compositors(self):
        """
        Returns a dictionary of diagrams that have a non-trivial
        composite, indexed by their compositor's name.

        More precisely, rather than :class:`Diagram` objects,
        the dictionary stores the :code:`shape` and :code:`mapping`
        data that allows to reconstruct them.

        Returns
        -------
        compositors : :class:`dict[dict]`
            The dictionary of composed diagrams.
        """
        return self._compositors

    @property
    def dim(self):
        """
        Returns the maximal dimension of a generator.

        Returns
        -------
        dim : :class:`int`
            The maximal dimension of a generator, or :code:`-1` if empty.
        """
        return max(self.by_dim, default=-1)

    @property
    def issimplicial(self):
        """
        Returns whether the diagrammatic sets is simplicial, that is,
        all its generators are simplex-shaped.

        Returns
        -------
        issimplicial : :class:`bool`
            :code:`True` if and only if the shape of every generator is
            a :class:`shapes.Simplex` object.
        """
        for x in self:
            if not isinstance(self.generators[x]['shape'], shapes.Simplex):
                return False
        return True

    @property
    def iscubical(self):
        """
        Returns whether the diagrammatic sets is cubical, that is,
        all its generators are cube-shaped.

        Returns
        -------
        iscubical : :class:`bool`
            :code:`True` if and only if the shape of every generator is
            a :class:`shapes.Cube` object.
        """
        for x in self:
            if not isinstance(self.generators[x]['shape'], shapes.Cube):
                return False
        return True

    def add(self, name, input=None, output=None, **kwargs):
        """
        Adds a generator and returns the diagram that maps the new
        generator into the diagrammatic set.

        The gluing of the generator is specified by a pair of round
        diagrams with identical boundaries, corresponding to the input
        and output diagrams of the new generator. If none are given,
        adds a point (0-dimensional generator).

        Arguments
        ---------
        name : :class:`hashable`
            Name to assign to the new generator.
        input : :class:`Diagram`, optional
            The input diagram of the new generator (default is :code:`None`)
        output : :class:`Diagram`, optional
            The output diagram of the new generator (default is :code:`None`)

        Keyword arguments
        -----------------
        color : multiple types
            Fill color when pictured as a node in string diagrams.
            If :code:`stroke` is not specified, this is
            also the color when pictured as a wire.
        stroke : multiple types
            Stroke color when pictured as a node, and color when pictured
            as a wire.
        draw_node : :class:`bool`
            If :code:`False`, no node is drawn when picturing the
            generator in string diagrams.
        draw_label : :class:`bool`
            If :code:`False`, no label is drawn when picturing the
            generator in string diagrams.

        Returns
        -------
        generator : :class:`Diagram`
            The diagram picking the new generator.

        Raises
        ------
        :class:`ValueError`
            If the name is already in use, or the input and output diagrams
            do not have round and matching boundaries.
        """
        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'name already in use'))
        for key in self._PRIVATE:
            kwargs.pop(key, None)

        if input is None:
            input = Diagram(self)
        if output is None:
            output = Diagram(self)
        for x in (input, output):
            utils.typecheck(x, {
                'type': Diagram,
                'st': lambda x: x.ambient == self,
                'why': 'not a diagram in {}'.format(repr(self))})

        boundary = Shape.atom(
                input.shape, output.shape, cospan=True)
        shape = boundary.target
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for x in input.shape:
            y = boundary.fst[x]
            mapping[y.dim][y.pos] = input.mapping[x.dim][x.pos]
        for x in output.shape:
            y = boundary.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = output.mapping[x.dim][x.pos]
            else:
                if mapping[y.dim][y.pos] != output.mapping[x.dim][x.pos]:
                    raise ValueError(utils.value_err(
                            output, 'boundary does not match '
                            'boundary of {}'.format(repr(input))))
        mapping[-1][0] = name
        self._update_generators(name, shape, mapping, **kwargs)

        return self[name]

    def add_simplex(self, name, *faces, **kwargs):
        """
        Variant of :meth:`add` for simplex-shaped generators.

        The gluing of the generator is specified by a number of
        :class:`SimplexDiagram` objects, corresponding to the faces
        of the new generator as listed by
        :class:`SimplexDiagram.simplex_face`.

        Arguments
        ---------
        name : :class:`hashable`
            Name to assign to the new generator.
        *faces : :class:`SimplexDiagram`
            The simplicial faces of the new generator.

        Keyword arguments
        -----------------
        **kwargs
            Same as :meth:`add`.

        Returns
        -------
        generator : :class:`SimplexDiagram`
            The diagram picking the new generator.

        Raises
        ------
        :class:`ValueError`
            If the name is already in use, or the faces do not have
            matching boundaries.
        """
        if len(faces) <= 1:
            return self.add(name, **kwargs)

        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'name already in use'))
        for key in self._PRIVATE:
            kwargs.pop(key, None)

        dim = len(faces) - 1
        for x in faces:
            utils.typecheck(x, {
                'type': SimplexDiagram,
                'st': lambda x: x.ambient == self and x.dim == dim-1,
                'why': 'expecting a {}-simplex in {}'.format(
                    str(dim - 1), repr(self))})

        shape = Shape.simplex(dim)
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for k, face in enumerate(faces):
            face_map = shape.simplex_face(k)
            for x in face:
                y = face_map[x]
                if mapping[y.dim][y.pos] is None:
                    mapping[y.dim][y.pos] = face.mapping[x.dim][x.pos]
                else:
                    if mapping[y.dim][y.pos] != face.mapping[x.dim][x.pos]:
                        raise ValueError(utils.value_err(
                            face, 'boundary of face does not '
                            'match other faces'))

        mapping[-1][0] = name
        self._update_generators(name, shape, mapping, **kwargs)

        return self[name]

    def add_cube(self, name, *faces, **kwargs):
        """
        Variant of :meth:`add` for cube-shaped generators.

        The gluing of the generator is specified by a number of
        :class:`CubeDiagram` objects, corresponding to the faces
        of the new generator as listed by
        :class:`CubeDiagram.cube_face`, in the order
        :code:`(0, '-')`, :code:`(0, '+')`, :code:`(1, '-')`,
        :code:`(1, '+')`, etc.

        Arguments
        ---------
        name : :class:`hashable`
            Name to assign to the new generator.
        *faces : :class:`CubeDiagram`
            The cubical faces of the new generator.

        Keyword arguments
        -----------------
        **kwargs
            Same as :meth:`add`.

        Returns
        -------
        generator : :class:`CubeDiagram`
            The diagram picking the new generator.

        Raises
        ------
        :class:`ValueError`
            If the name is already in use, or the faces do not have
            matching boundaries.
        """
        if len(faces) % 2 == 1:
            raise ValueError(utils.value_err(
                faces, 'expecting an even number of faces'))
        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'name already in use'))
        for key in self._PRIVATE:
            kwargs.pop(key, None)

        dim = int(len(faces) / 2)
        for x in faces:
            utils.typecheck(x, {
                'type': CubeDiagram,
                'st': lambda x: x.ambient == self and x.dim == dim-1,
                'why': 'expecting a {}-cube in {}'.format(
                    str(dim - 1), repr(self))})

        shape = Shape.cube(dim)
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for n, face in enumerate(faces):
            k = int(n/2)
            sign = utils.mksign(n % 2)
            face_map = shape.cube_face(k, sign)
            for x in face:
                y = face_map[x]
                if mapping[y.dim][y.pos] is None:
                    mapping[y.dim][y.pos] = face.mapping[x.dim][x.pos]
                else:
                    if mapping[y.dim][y.pos] != face.mapping[x.dim][x.pos]:
                        raise ValueError(utils.value_err(
                            face, 'boundary of face does not '
                            'match other faces'))

        mapping[-1][0] = name
        self._update_generators(name, shape, mapping, **kwargs)

        return self[name]

    def invert(self, generatorname,
               inversename=None,
               rinvertorname=None,
               linvertorname=None,
               **kwargs):
        """
        Adds a weak inverse for a generator, together
        with left and right invertors that witness the
        inversion, and returns them as diagrams.

        Both the inverse and the invertors can be given custom names.
        If the generator to be inverted is named :code:`'a'`, the
        default names are

        - :code:`'a⁻¹'` for the inverse,
        - :code:`'inv(a, a⁻¹)'` for the right invertor,
        - :code:`'inv(a⁻¹, a)'` for the left invertor.

        In the theory of diagrammatic sets, weak invertibility would
        correspond to the situation where the invertors themselves
        are weakly invertible, coinductively.
        In the implementation, we take an "invert when necessary"
        approach, where invertors are not invertible by default, and
        should be inverted when needed.

        Notes
        -----
        The right invertor for the generator is the left invertor
        for its inverse, and the left invertor for the generator is the
        right invertor for its inverse.

        Arguments
        ---------
        generatorname : :class:`hashable`
            Name of the generator to invert.
        inversename : :class:`hashable`, optional
            Name assigned to the inverse.
        rinvertorname : :class:`hashable`, optional
            Name assigned to the right invertor.
        linvertorname : :class:`hashable`, optional
            Name assigned to the left invertor.

        Keyword arguments
        -----------------
        **kwargs
            Passed to :meth:`add` when adding the inverse.

        Returns
        -------
        inverse : :class:`Diagram`
            The diagram picking the inverse.
        rinvertor : :class:`Diagram`
            The diagram picking the right invertor.
        linvertor : :class:`Diagram`
            The diagram picking the left invertor.

        Raises
        ------
        :class:`ValueError`
            If the generator is already inverted, or 0-dimensional.
        """
        if isinstance(generatorname, Diagram):
            generatorname = generatorname.name
        generator = self[generatorname]
        if generator.dim == 0:
            raise ValueError(utils.value_err(
                generatorname, 'cannot invert 0-cell'))
        if generator.isinvertiblecell:
            raise ValueError(utils.value_err(
                generatorname, 'already inverted'))

        if inversename is None:
            inversename = '{}⁻¹'.format(str(generatorname))
        if rinvertorname is None:
            rinvertorname = 'inv({}, {})'.format(
                    str(generatorname), str(inversename))
        if linvertorname is None:
            linvertorname = 'inv({}, {})'.format(
                    str(inversename), str(generatorname))
        for x in (inversename, rinvertorname, linvertorname):
            if x in self.generators:
                raise ValueError(utils.value_err(
                    generatorname, 'name already in use'))

        inverse = self.add(
                inversename,
                generator.output,
                generator.input,
                **kwargs)
        rinvertor = self.add(
                rinvertorname,
                generator.paste(inverse),
                generator.input.unit())
        linvertor = self.add(
                linvertorname,
                inverse.paste(generator),
                generator.output.unit())

        self._generators[generatorname].update({
            'inverse': inversename,
            'rinvertor': rinvertorname,
            'linvertor': linvertorname})
        self._generators[inversename].update({
            'inverse': generatorname,
            'rinvertor': linvertorname,
            'linvertor': rinvertorname})
        self._generators[rinvertorname].update({
            'inverts': (generatorname, inversename)})
        self._generators[linvertorname].update({
            'inverts': (inversename, generatorname)})

        return inverse, rinvertor, linvertor

    def make_inverses(self, generatorname1, generatorname2,
                      rinvertorname=None,
                      linvertorname=None):
        """
        Makes two generators each other's weak inverse by adding
        invertors, and returns the invertors.

        In what follows, "right/left" invertors are relative to
        the first generator.
        Both invertors can be given custom names.
        If the generators are named :code:`'a'`, :code:`'b'`, the
        default names for the invertors are

        - :code:`'inv(a, b)'` for the right invertor,
        - :code:`'inv(b, a)'` for the left invertor.

        In the theory of diagrammatic sets, weak invertibility would
        correspond to the situation where the invertors themselves
        are weakly invertible, coinductively.
        In the implementation, we take an "invert when necessary"
        approach, where invertors are not invertible by default, and
        should be inverted when needed.

        Arguments
        ---------
        generatorname1 : :class:`hashable`
            Name of the first generator.
        generatorname2 : :class:`hashable`, optional
            Name of the second generator.
        rinvertorname : :class:`hashable`, optional
            Name assigned to the right invertor.
        linvertorname : :class:`hashable`, optional
            Name assigned to the left invertor.

        Returns
        -------
        rinvertor : :class:`Diagram`
            The diagram picking the right invertor.
        linvertor : :class:`Diagram`
            The diagram picking the left invertor.

        Raises
        ------
        :class:`ValueError`
            If the generators are already inverted, or 0-dimensional,
            or do not have compatible boundaries.
        """
        for x in (generatorname1, generatorname2):
            if isinstance(x, Diagram):
                x = x.name
        generator1 = self[generatorname1]
        generator2 = self[generatorname2]

        selfinverse = generatorname1 == generatorname2
        for x in (generator1, generator2):
            if x.isinvertiblecell:
                raise ValueError(utils.value_err(
                    x.name, 'already inverted'))

        rpaste = generator1.paste(generator2)
        if rinvertorname is None:
            rinvertorname = 'inv({}, {})'.format(
                    str(generatorname1), str(generatorname2))

        if selfinverse:
            lpaste = rpaste
            linvertorname = rinvertorname
        else:
            lpaste = generator2.paste(generator1)
            if linvertorname is None:
                linvertorname = 'inv({}, {})'.format(
                        str(generatorname2), str(generatorname1))

        for x in (rinvertorname, linvertorname):
            if x in self.generators:
                raise ValueError(utils.value_err(
                    x, 'name already in use'))

        rinvertor = self.add(
                rinvertorname,
                rpaste,
                generator1.input.unit())
        if selfinverse:
            linvertor = rinvertor
        else:
            linvertor = self.add(
                    linvertorname,
                    lpaste,
                    generator2.input.unit())

        self._generators[generatorname1].update({
            'inverse': generatorname2,
            'rinvertor': rinvertorname,
            'linvertor': linvertorname})
        self._generators[rinvertorname].update({
            'inverts': (generatorname1, generatorname2)})
        if not selfinverse:
            self._generators[generatorname2].update({
                'inverse': generatorname1,
                'rinvertor': linvertorname,
                'linvertor': rinvertorname})
            self._generators[linvertorname].update({
                'inverts': (generatorname2, generatorname1)})

        return rinvertor, linvertor

    def compose(self, diagram,
                name=None, compositorname=None,
                **kwargs):
        """
        Given a round diagram, adds a weak composite for it,
        together with a compositor witnessing the composition, and
        returns them as diagrams.

        Both the composite and the compositor can be given custom names.
        If the diagram to be composed is named :code:`'a'`, the
        default names are

        - :code:`'⟨a⟩'` for the composite,
        - :code:`'comp(a)'` for the compositor.

        In the theory of diagrammatic sets, a weak composite is
        witnessed by a weakly invertible compositor.
        In the implementation, we take an "invert when necessary"
        approach, where compositors are not invertible by default, and
        should be inverted when needed.

        Notes
        -----
        A cell (a diagram whose shape is an atom) is treated as already
        having itself as a composite, witnessed by a unit cell; this
        method can only be used on non-atomic diagrams.

        Arguments
        ---------
        diagram : :class:`Diagram`
            The diagram to compose.
        name : :class:`hashable`, optional
            Name of the weak composite.
        compositorname : :class:`hashable`, optional
            Name of the compositor.

        Keyword arguments
        -----------------
        **kwargs
            Passed to :meth:`add` when adding the composite.

        Returns
        -------
        composite : :class:`Diagram`
            The diagram picking the composite.
        compositor : :class:`Diagram`
            The diagram picking the compositor.

        Raises
        ------
        :class:`ValueError`
            If the diagram is not round, or already has a composite.
        """
        utils.typecheck(diagram, {
            'type': Diagram,
            'st': lambda x: x.shape.isround,
            'why': 'composable diagrams must have round shape'})

        if diagram.hascomposite:
            raise ValueError(utils.value_err(
                diagram, 'already has a composite'))

        if name is None:
            name = '⟨{}⟩'.format(str(diagram.name))
        if compositorname is None:
            compositorname = 'comp({})'.format(str(diagram.name))
        for x in (name, compositorname):
            if x in self.generators:
                raise ValueError(utils.value_err(
                    x, 'name already in use'))

        composite = self.add(
                name,
                diagram.input,
                diagram.output,
                **kwargs)
        compositor = self.add(
                compositorname,
                diagram,
                composite)

        self._generators[compositorname].update({
            'composite': name})
        self._compositors.update({
            compositorname: {
                'shape': diagram.shape,
                'mapping': diagram.mapping}
            })

        return composite, compositor

    def make_composite(self, generatorname, diagram,
                       compositorname=None):
        """
        Given a generator and a round diagram, it makes the first
        the weak composite of the second by adding a compositor, and
        returns the compositor as a diagram.

        The compositor can be given a custom name.
        If the diagram to be composed is named :code:`'a'`, the
        default name is :code:`'comp(a)'`.

        In the theory of diagrammatic sets, a weak composite is
        witnessed by a weakly invertible compositor.
        In the implementation, we take an "invert when necessary"
        approach, where compositors are not invertible by default, and
        should be inverted when needed.

        Notes
        -----
        A cell (a diagram whose shape is an atom) is treated as already
        having itself as a composite, witnessed by a unit cell; this
        method can only be used on non-atomic diagrams.

        Arguments
        ---------
        generatorname : :class:`hashable`
            Name of the generator that should be its composite.
        diagram : :class:`Diagram`
            The diagram to compose.
        compositorname : :class:`hashable`, optional
            Name of the compositor.

        Returns
        -------
        compositor : :class:`Diagram`
            The diagram picking the compositor.

        Raises
        ------
        :class:`ValueError`
            If the diagram is not round, or already has a composite, or
            the diagram and the generator do not have matching boundaries.
        """
        if isinstance(generatorname, Diagram):
            generatorname = generatorname.name
        utils.typecheck(diagram, {
            'type': Diagram,
            'st': lambda x: x.shape.isround,
            'why': 'composable diagrams must have round shape'})

        if diagram.hascomposite:
            raise ValueError(utils.value_err(
                diagram, 'already has a composite'))

        if compositorname is None:
            compositorname = 'comp({})'.format(str(diagram.name))
        if compositorname in self.generators:
            raise ValueError(utils.value_err(
                compositorname, 'name already in use'))

        generator = self[generatorname]
        compositor = self.add(
                compositorname,
                diagram,
                generator)

        self._generators[compositorname].update({
            'composite': generatorname})
        self._compositors.update({
            compositorname: {
                'shape': diagram.shape,
                'mapping': diagram.mapping}
            })

        return compositor

    def remove(self, generatorname):
        """
        Removes a generator, together with all other generators
        that depend on it.

        Arguments
        ---------
        generatorname : :class:`hashable`
            Name of the generator to remove.
        """
        if isinstance(generatorname, Diagram):
            generatorname = generatorname.name
        to_remove = [*self.generators[generatorname]['cofaces']]
        for x in to_remove:
            self.remove(x)

        dim = self[generatorname].dim
        self._by_dim[dim].remove(generatorname)
        if len(self._by_dim[dim]) == 0:
            self._by_dim.pop(dim, None)

        if 'inverse' in self.generators[generatorname]:
            inverse = self.generators[generatorname]['inverse']
            self.generators[inverse].pop('inverse', None)
            self.generators[inverse].pop('linvertor', None)
            self.generators[inverse].pop('rinvertor', None)

        if 'inverts' in self.generators[generatorname]:
            fst, snd = self.generators[generatorname]['inverts']
            linvertor = self.generators[fst]['linvertor']
            for x in (fst, snd):
                self.generators[x].pop('inverse', None)
                self.generators[x].pop('linvertor', None)
                self.generators[x].pop('rinvertor', None)
            self.generators[linvertor].pop('inverts', None)

        if 'composite' in self.generators[generatorname]:
            # This does not remove the composite!
            self.compositors.pop(generatorname, None)

        for x in self.generators[generatorname]['faces']:
            self.generators[x]['cofaces'].remove(generatorname)

        self.generators.pop(generatorname, None)

    def update(self, generatorname, **kwargs):
        """
        Updates the optional arguments of a generator.

        Arguments
        ---------
        generatorname : :class:`hashable`
            Name of the generator to update.

        Keyword arguments
        -----------------
        **kwargs
            Any arguments to update.

        Raises
        ------
        :class:`AttributeError`
            If the optional argument uses a private keyword.
        """
        if isinstance(generatorname, Diagram):
            generatorname = generatorname.name
        for key in self._PRIVATE:
            if key in kwargs:
                raise AttributeError(key, 'private attribute')
        self.generators[generatorname].update(kwargs)

    def copy(self):
        """
        Returns a copy of the object.

        Returns
        -------
        copy : :class:`DiagSet`
            A copy of the object.
        """
        new = DiagSet()
        new._generators = self.generators.copy()
        new._by_dim = self.by_dim.copy()
        new._compositors = self.compositors.copy()
        return new

    @staticmethod
    def yoneda(shape):
        """
        Alternative constructor creating a diagrammatic set from
        a :class:`shapes.Shape`.

        Mathematically, diagrammatic sets are certain sheaves on the
        category of shapes and maps of shapes; this constructor
        implements the Yoneda embedding of a shape.
        This has an `n`-dimensional generator for each `n`-dimensional
        element of the shape.

        Arguments
        ---------
        shape : :class:`shapes.Shape`
            A shape.

        Returns
        -------
        yoneda : :class:`DiagSet`
            The Yoneda-embedded shape.
        """
        utils.typecheck(shape, {'type': Shape})
        yoneda = DiagSet()
        for x in shape:
            atom = shape.atom_inclusion(x)
            yoneda._generators.update({
                x: {
                    'shape': atom.source,
                    'mapping': atom.mapping,
                    'faces': {y for y in shape.faces(x)},
                    'cofaces': {y for y in shape.cofaces(x)}
                    }
                })
        for n, n_size in enumerate(shape.size):
            yoneda._by_dim.update({
                n: {
                    El(n, k) for k in range(n_size)
                    }
                })
        return yoneda

    # Internal methods
    def _update_generators(self, name, shape, mapping, **kwargs):
        self._generators.update({
                name: {
                    'shape': shape,
                    'mapping': mapping,
                    'faces': set(),
                    'cofaces': set(),
                    **kwargs
                    }
                })

        if shape.dim in self._by_dim:
            self._by_dim[shape.dim].add(name)
        else:
            self._by_dim[shape.dim] = {name}

        if shape.dim > 0:
            for x in mapping[-2]:
                self.generators[x]['cofaces'].add(name)
                self.generators[name]['faces'].add(x)


class Diagram:
    """
    Class for diagrams, that is, mappings from a shape to an
    "ambient" diagrammatic set.

    To create a diagram, we start from *generators*
    of a diagrammatic set, returned by the :meth:`DiagSet.add`
    method or requested with indexer operators.

    Then we produce other diagrams in two main ways:

    - pulling back a diagram along a map of shapes
      (:meth:`pullback`), or
    - pasting together two diagrams along their boundaries
      (:meth:`paste`, :meth:`to_inputs`, :meth:`to_outputs`).

    In practice, the direct use of :meth:`pullback`, which requires
    an explicit shape map, can be avoided in common cases by using
    :meth:`unit`, :meth:`lunitor`, :meth:`runitor`, or the
    specialised :class:`SimplexDiagram.simplex_degeneracy`,
    :class:`CubeDiagram.cube_degeneracy`, and
    :class:`CubeDiagram.cube_connection` methods.

    Notes
    -----
    Initialising a :class:`Diagram` directly creates an empty
    diagram in a given diagrammatic set.

    Arguments
    ---------
    ambient : :class:`DiagSet`
        The ambient diagrammatic set.
    """
    def __init__(self, ambient):
        utils.typecheck(ambient, {'type': DiagSet})

        self._shape = Shape.empty()
        self._ambient = ambient
        self._mapping = []
        self._name = ''

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
        return isinstance(other, Diagram) and \
                self.shape == other.shape and \
                self.ambient == other.ambient and \
                self.mapping == other.mapping

    def __len__(self):
        return len(self.shape)

    def __getitem__(self, element):
        if element in self.shape:
            return self.mapping[element.dim][element.pos]
        raise ValueError(utils.value_err(
            element, 'not an element of the shape'))

    def __contains__(self, item):
        return item in self.shape

    def __iter__(self):
        return iter(self.shape)

    @property
    def name(self):
        """
        Returns the name of the diagram.

        Returns
        -------
        name : :class:`hashable`
            The name of the diagram.
        """
        return self._name

    @property
    def shape(self):
        """
        Returns the shape of the diagram.

        Returns
        -------
        shape : :class:`shapes.Shape`
            The shape of the diagram.
        """
        return self._shape

    @property
    def ambient(self):
        """
        Returns the ambient diagrammatic set.

        Returns
        -------
        ambient : :class:`DiagSet`
            The ambient diagrammatic set.
        """
        return self._ambient

    @property
    def mapping(self):
        """
        Returns the data specifying the mapping of shape elements to
        generators.

        The mapping is specified as a list of lists, similar to
        :class:`ogposets.OgMap`, in the following way:
        :code:`mapping[n][k] == s` if the diagram sends :code:`El(n, k)`
        to the generator named :code:`s`.

        Returns
        -------
        mapping : :class:`list[list[hashable]]`
            The data specifying the diagram's assignment.
        """
        return self._mapping

    @property
    def layers(self):
        """
        Returns the layering of the diagram corresponding to the current
        layering of the shape.

        Returns
        -------
        layers : :class:`list[Diagram]`
            The current layering.
        """
        if not hasattr(self.shape, '_layering'):
            return [self]
        return [
                self.pullback(f, 'layer {} of {}'.format(
                    str(n), self.name))
                for n, f in enumerate(self.shape._layering)]

    @property
    def rewrite_steps(self):
        """
        Returns the sequence of rewrite steps associated to the current
        layering of the diagram.

        The :code:`0`-th rewrite step is the input boundary of the diagram.
        For :code:`n > 0`, the :code:`n`-th rewrite step is the output
        boundary of the :code:`(n-1)`-th layer.

        Returns
        -------
        rewrite_steps : :class:`list[Diagram]`
            The current sequence of rewrite steps.
        """
        rewrite_steps = [
                *[layer.input for layer in self.layers],
                self.layers[-1].output
                ]
        for n, step in enumerate(rewrite_steps):
            step.rename('step {} of {}'.format(
                    str(n), self.name))
        return rewrite_steps

    @property
    def dim(self):
        """
        Shorthand for :code:`shape.dim`.
        """
        return self.shape.dim

    @property
    def isdegenerate(self):
        """
        Returns whether the diagram is *degenerate*, that is, its
        image has dimension strictly lower than the dimension of its shape.

        Returns
        -------
        isdegenerate : :class:`bool`
            :code:`True` if and only if the diagram is degenerate.
        """
        if self.dim <= 0:
            return False
        for x in self.mapping[-1]:
            if self.ambient[x].dim == self.dim:
                return False
        return True

    @property
    def isround(self):
        """
        Shorthand for :code:`shape.isround`.
        """
        return self.shape.isround

    @property
    def iscell(self):
        """
        Shorthand for :code:`shape.isatom` (a *cell* is a diagram
        whose shape is an atom).
        """
        return self.shape.isatom

    @property
    def isinvertiblecell(self):
        """
        Returns whether the diagram is an invertible cell.

        A cell is invertible if either

        - it is degenerate, or
        - its image is an invertible generator.

        Returns
        -------
        isinvertiblecell : :class:`bool`
            :code:`True` if and only if the diagram is an invertible cell.
        """
        if self.iscell:
            if self.isdegenerate:
                return True
            if 'inverse' in self.ambient.generators[self.mapping[-1][0]]:
                return True
        return False

    @property
    def hascomposite(self):
        """
        Returns whether the diagram has a composite.

        Returns
        -------
        hascomposite : :class:`bool`
            :code:`True` if and only if the diagram has a composite.
        """
        if self.iscell:
            return True
        if self._find_compositor() is None:
            return False
        return True

    def rename(self, name):
        """
        Renames the diagram.

        Arguments
        ---------
        name : :class:`hashable`
            The new name for the diagram.
        """
        self._name = name

    # Methods for creating new diagrams
    def paste(self, other, dim=None, **params):
        """
        Given two diagrams and :code:`k` such that the output
        :code:`k`-boundary of the first is equal to the input
        :code:`k`-boundary of the second, returns their pasting along
        the matching boundaries.

        Arguments
        ---------
        fst : :class:`Diagram`
            The first diagram.
        snd : :class:`Diagram`
            The second diagram.
        dim : :class:`int`, optional
            The dimension of the boundary along which they will be pasted
            (default is :code:`min(fst.dim, snd.dim) - 1`).

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to also return the cospan of inclusions of the two
            diagrams' shapes into the pasting (default is :code:`False`).

        Returns
        -------
        paste : :class:`Diagram`
            The pasted diagram.
        paste_cospan : :class:`ogposets.OgMapPair`, optional
            The cospan of inclusions of the two diagrams' shapes into
            their pasting.

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match.
        """
        cospan = params.get('cospan', False)

        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        if dim is None:
            dim = min(self.dim, other.dim) - 1

        paste_cospan = Shape.paste(
                self.shape, other.shape, dim, cospan=True)

        shape = paste_cospan.target
        mapping = Diagram._paste_fill_mapping(
                self, other, paste_cospan)
        name = '({}) #{} ({})'.format(
                str(self.name), str(dim), str(other.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        if cospan:
            return pasted, paste_cospan
        return pasted

    def to_outputs(self, positions, other, dim=None, **params):
        """
        Returns the pasting of another diagram along a round subshape of
        the output :code:`k`-boundary, specified by the positions of its
        :code:`k`-dimensional elements.

        Arguments
        ---------
        positions : :class:`list[int]` | :class:`int`
            The positions of the outputs along which to paste. If given
            an integer :code:`n`, interprets it as the list :code:`[n]`.
        other : :class:`Diagram`
            The other diagram to paste.
        dim : :class:`int`, optional
            The dimension of the boundary along which to paste
            (default is :code:`self.dim - 1`)

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the two diagrams'
            shapes into the pasting (default is :code:`False`).

        Returns
        -------
        to_outputs : :class:`Shape`
            The pasted diagram.
        paste_cospan : :class:`ogposets.OgMapPair`, optional
            The cospan of inclusions of the two diagrams' shapes into
            their pasting.

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match, or the pasting does not produce
            a well-formed shape.
        """
        cospan = params.get('cospan', False)
        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        if isinstance(positions, int):
            positions = [positions]
        if dim is None:
            dim = self.dim-1
        paste_cospan = self.shape.to_outputs(
                positions, other.shape, dim, cospan=True)

        shape = paste_cospan.target
        mapping = Diagram._paste_fill_mapping(
                self, other, paste_cospan)
        name = '{{{} -> {}{}}}({})'.format(
                str(other.name), str(dim),
                str(sorted(positions)), str(self.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        if cospan:
            return pasted, paste_cospan
        return pasted

    def to_inputs(self, positions, other, dim=None, **params):
        """
        Returns the pasting of another diagram along a round subshape of
        the input :code:`k`-boundary, specified by the positions of its
        :code:`k`-dimensional elements.

        Arguments
        ---------
        positions : :class:`list[int]` | :class:`int`
            The positions of the inputs along which to paste. If given
            an integer :code:`n`, interprets it as the list :code:`[n]`.
        other : :class:`Diagram`
            The other diagram to paste.
        dim : :class:`int`, optional
            The dimension of the boundary along which to paste
            (default is :code:`self.dim - 1`)

        Keyword arguments
        -----------------
        cospan : :class:`bool`
            Whether to return the cospan of inclusions of the two diagrams'
            shapes into the pasting (default is :code:`False`).

        Returns
        -------
        to_inputs : :class:`Shape`
            The pasted diagram.
        paste_cospan : :class:`ogposets.OgMapPair`, optional
            The cospan of inclusions of the two diagrams' shapes into
            their pasting.

        Raises
        ------
        :class:`ValueError`
            If the boundaries do not match, or the pasting does not produce
            a well-formed shape.
        """
        cospan = params.get('cospan', False)

        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        if isinstance(positions, int):
            positions = [positions]
        if dim is None:
            dim = self.dim-1
        paste_cospan = self.shape.to_inputs(
                positions, other.shape, dim, cospan=True)

        shape = paste_cospan.target
        mapping = Diagram._paste_fill_mapping(
                other, self, paste_cospan)
        name = '({}){{{}{} <- {}}}'.format(
                str(self.name), str(sorted(positions)),
                str(dim), str(other.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        if cospan:
            return pasted, paste_cospan
        return pasted

    def rewrite(self, positions, diagram):
        """
        Returns the diagram representing the application of a
        higher-dimensional rewrite to a subdiagram, specified
        by the positions of its top-dimensional elements.

        This is in fact an alias for :meth:`to_outputs` in the top
        dimension, reflecting the intuitions of higher-dimensional
        rewriting in this situation.

        Arguments
        ---------
        positions : :class:`list[int]` | :class:`int`
            The positions of the top-dimensional elements to rewrite.
            If given an integer :code:`n`, interprets it as the list
            :code:`[n]`.
        diagram : :class:`Diagram`
            The diagram representing the rewrite to apply.

        Returns
        -------
        rewrite : :class:`Shape`
            The diagram representing the application of the rewrite to
            the given positions.
        """
        return self.to_outputs(positions, diagram, self.dim)

    def pullback(self, shapemap, name=None):
        """
        Returns the pullback of the diagram along a shape map.

        Arguments
        ---------
        shapemap : :class:`shapes.ShapeMap`
            The map along which to pull back.
        name : :class:`hashable`, optional
            The name to give to the new diagram.

        Returns
        -------
        pullback : :class:`Diagram`
            The pulled back diagram.

        Raises
        ------
        :class:`ValueError`
            If the target of the map is not equal to the diagram shape.
        """
        utils.typecheck(shapemap, {
            'type': ShapeMap,
            'st': lambda x: x.target == self.shape,
            'why': 'target does not match diagram shape'})

        shape = shapemap.source
        mapping = [
                [
                    self.mapping[x.dim][x.pos]
                    for x in n_data
                ]
                for n_data in shapemap.mapping]
        if name is None:
            name = 'pullback of {}'.format(str(self.name))
        return Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)

    def boundary(self, sign, dim=None):
        """
        Returns the boundary of a given orientation and dimension.

        This is, by definition, the pullback of a diagram along the
        inclusion map :code:`self.shape.boundary(sign, dim)`.

        Arguments
        ---------
        sign : :class:`str`
            Orientation: :code:`'-'` for input, :code:`'+'` for output.
        dim : :class:`int`, optional
            Dimension of the boundary (default is :code:`self.dim - 1`).

        Returns
        -------
        boundary : :class:`Diagram`
            The requested boundary.
        """
        if dim is None:
            dim = self.dim - 1
        sign = utils.mksign(sign)
        name = '∂[{},{}]({})'.format(
                sign, str(dim), str(self.name))
        return self.pullback(self.shape.boundary(
            sign, dim), name)

    @property
    def input(self):
        """
        Alias for :code:`boundary('-')`.
        """
        return self.boundary('-')

    @property
    def output(self):
        """
        Alias for :code:`boundary('+')`.
        """
        return self.boundary('+')

    def unit(self):
        """
        Returns the unit on the diagram: a degenerate diagram one
        dimension higher, with input and output equal to the diagram.

        This is, by definition, the pullback of the diagram along
        :code:`self.shape.inflate()`.

        Returns
        -------
        unit : :class:`Diagram`
            The unit diagram.
        """
        name = '1({})'.format(str(self.name))
        return self.pullback(self.shape.inflate(), name)

    def lunitor(self, sign='-', positions=None):
        """
        Returns a left unitor on the diagram:
        a degenerate diagram one dimension higher, with one
        boundary equal to the diagram, and the other equal to the
        diagram with units pasted to some of its inputs.

        Arguments
        ---------
        sign : :class:`str`, optional
            The boundary on which the units are: :code:`'-'` (default)
            for input, :code:`'+'` for output.
        positions : :class:`list[int]` | :class:`int`
            The positions of the inputs to which a unit is attached
            (default is *all* of the inputs). If given
            an integer :code:`n`, interprets it as the list :code:`[n]`.

        Returns
        -------
        lunitor : :class:`Diagram`
            The left unitor diagram.

        Raises
        ------
        :class:`ValueError`
            If the positions do not correspond to inputs.
        """
        sign = utils.mksign(sign)
        dim = self.dim - 1

        all = self.shape.all()
        output = all.boundary('+')

        collapsed = output
        if positions is not None:
            if isinstance(positions, int):
                positions = [positions]
            input = all.boundary('-')
            notcollapsed = GrSubset(
                    GrSet(
                        *[El(dim, k) for k in positions]),
                    self.shape)
            if not notcollapsed.issubset(input):
                raise ValueError(utils.value_err(
                    positions, 'not in the input boundary'))

            notcollapsed = notcollapsed.closure()
            collapsed = collapsed.union(
                    input.difference(notcollapsed).closure())

        unitor_map = self.shape.inflate(
                collapsed)
        if positions is None:
            name = 'λ({})'.format(str(self.name))
        else:
            name = 'λ{}({})'.format(
                    str(sorted(positions)),
                    str(self.name))

        if sign == '-':
            unitor_map = unitor_map.dual(self.dim + 1)
        else:
            name = '({})⁻¹'.format(name)
        return self.pullback(unitor_map, name)

    def runitor(self, sign='-', positions=None):
        """
        Returns a right unitor on the diagram:
        a degenerate diagram one dimension higher, with one
        boundary equal to the diagram, and the other equal to the
        diagram with units pasted to some of its outputs.

        Arguments
        ---------
        sign : :class:`str`, optional
            The boundary on which the units are: :code:`'-'` (default)
            for input, :code:`'+'` for output.
        positions : :class:`list[int]` | :class:`int`
            The positions of the outputs to which a unit is attached
            (default is *all* of the outputs). If given
            an integer :code:`n`, interprets it as the list :code:`[n]`.

        Returns
        -------
        runitor : :class:`Diagram`
            The right unitor diagram.

        Raises
        ------
        :class:`ValueError`
            If the positions do not correspond to outputs.
        """
        sign = utils.mksign(sign)
        dim = self.dim - 1

        all = self.shape.all()
        input = all.boundary('-')

        collapsed = input
        if positions is not None:
            if isinstance(positions, int):
                positions = [positions]
            output = all.boundary('+')
            notcollapsed = GrSubset(
                    GrSet(
                        *[El(dim, k) for k in positions]),
                    self.shape)
            if not notcollapsed.issubset(output):
                raise ValueError(utils.value_err(
                    positions, 'not in the output boundary'))

            notcollapsed = notcollapsed.closure()
            collapsed = collapsed.union(
                    output.difference(notcollapsed).closure())

        unitor_map = self.shape.inflate(
                collapsed)
        if positions is None:
            name = 'ρ({})'.format(str(self.name))
        else:
            name = 'ρ{}({})'.format(
                    str(sorted(positions)),
                    str(self.name))

        if sign == '+':
            unitor_map = unitor_map.dual(self.dim + 1)
            name = '({})⁻¹'.format(name)
        return self.pullback(unitor_map, name)

    @property
    def inverse(self):
        """
        Returns the inverse of an invertible cell.

        Returns
        -------
        inverse : :class:`Diagram`
            The inverse cell.

        Raises
        ------
        :class:`ValueError`
            If the diagram is not an invertible cell.
        """
        if not self.isinvertiblecell:
            raise ValueError(utils.value_err(
                self, 'not an invertible cell'))

        if not self.isdegenerate:
            top = self.mapping[-1][0]
            top_inv = self.ambient.generators[top]['inverse']

            if self.shape == self.ambient.generators[top]['shape']:
                return self.ambient[top_inv]

        reordering = Shape.dual(self.shape, self.dim,
                                reordering=True)
        shape = reordering.source
        mapping = [
                [self[x] for x in n_data]
                for n_data in reordering.mapping
                ]
        if not self.isdegenerate:
            mapping[-1][0] = top_inv
        name = '({})⁻¹'.format(self.name)

        return Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)

    @property
    def rinvertor(self):
        """
        Returns the right invertor for an invertible cell.

        Returns
        -------
        rinvertor : :class:`Diagram`
            The right invertor.

        Raises
        ------
        :class:`ValueError`
            If the diagram is not an invertible cell.
        """
        if not self.isinvertiblecell:
            raise ValueError(utils.value_err(
                self, 'not an invertible cell'))

        top = self.mapping[-1][0]
        if not self.isdegenerate:
            top_rinvertor = self.ambient.generators[top]['rinvertor']

            if self.shape == self.ambient.generators[top]['shape']:
                return self.ambient[top_rinvertor]

        inverse = self.inverse
        rpaste, rpaste_cospan = self.paste(inverse, cospan=True)
        unit = self.input.unit()

        atom_cospan = rpaste.shape.atom(unit.shape, cospan=True)
        shape = atom_cospan.target
        mapping = Diagram._paste_fill_mapping(
                rpaste, unit, atom_cospan)

        if not self.isdegenerate:
            mapping[-1][0] = top_rinvertor
        else:
            mapping[-1][0] = top

        name = 'inv({}, {})'.format(self.name, inverse.name)

        return Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)

    @property
    def linvertor(self):
        """
        Returns the left invertor for an invertible cell.

        Returns
        -------
        linvertor : :class:`Diagram`
            The left invertor.

        Raises
        ------
        :class:`ValueError`
            If the diagram is not an invertible cell.
        """
        if not self.isinvertiblecell:
            raise ValueError(utils.value_err(
                self, 'not an invertible cell'))

        top = self.mapping[-1][0]
        if not self.isdegenerate:
            top_linvertor = self.ambient.generators[top]['linvertor']

            if self.shape == self.ambient.generators[top]['shape']:
                return self.ambient[top_linvertor]

        inverse = self.inverse
        lpaste, lpaste_cospan = inverse.paste(self, cospan=True)
        unit = self.output.unit()

        atom_cospan = lpaste.shape.atom(unit.shape, cospan=True)
        shape = atom_cospan.target
        mapping = Diagram._paste_fill_mapping(
                lpaste, unit, atom_cospan)

        if not self.isdegenerate:
            mapping[-1][0] = top_linvertor
        else:
            mapping[-1][0] = top

        name = 'inv({}, {})'.format(inverse.name, self.name)

        return Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)

    @property
    def composite(self):
        """
        Returns the composite of the diagram, if it exists.

        Returns
        -------
        composite : :class:`Diagram`
            The composite.

        Raises
        ------
        :class:`ValueError`
            If the diagram does not have a composite.
        """
        if not self.hascomposite:
            raise ValueError(utils.value_err(
                self, 'does not have a composite'))

        if self.iscell:
            return self
        compositorname = self._find_compositor()
        name = self.ambient.generators[compositorname]['composite']
        return self.ambient[name]

    @property
    def compositor(self):
        """
        Returns the compositor of the diagram, if it exists.

        Returns
        -------
        compositor : :class:`Diagram`
            The compositor.

        Raises
        ------
        :class:`ValueError`
            If the diagram does not have a composite.
        """
        if not self.hascomposite:
            raise ValueError(utils.value_err(
                self, 'does not have a compositor'))

        if self.iscell:
            return self.unit()
        compositorname = self._find_compositor()
        return self.ambient[compositorname]

    def generate_layering(self):
        """
        Assigns a layering to the diagram, iterating through all
        the layerings, and returns it.

        Returns
        -------
        layers : :class:`list[Diagram]`
            The generated layering.
        """
        self.shape.generate_layering()
        return self.layers

    def hasse(self, **params):
        """
        Bound version of :meth:`hasse.draw`.

        Calling :code:`x.hasse(**params)` is equivalent to calling
        :code:`hasse.draw(x, **params)`.
        """
        from rewalt.hasse import draw
        return draw(self, **params)

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

    # Alternative constructors
    @staticmethod
    def yoneda(shapemap, name=None):
        """
        Alternative constructor creating a diagram from
        a :class:`shapes.ShapeMap`.

        Mathematically, diagrammatic sets are certain sheaves on the
        category of shapes and maps of shapes; this constructor
        implements the Yoneda embedding of a map of shapes.

        Arguments
        ---------
        shapemap : :class:`shapes.Shape`
            A map of shapes.
        name : :class:`hashable`, optional
            The name of the generated diagram.

        Returns
        -------
        yoneda : :class:`Diagram`
            The Yoneda-embedded map.
        """
        utils.typecheck(shapemap, {
            'type': ShapeMap})
        return Diagram._new(
                shapemap.source,
                DiagSet.yoneda(shapemap.target),
                shapemap.mapping,
                name)

    @staticmethod
    def with_layers(fst, *layers):
        """
        Given a non-zero number of diagrams that can be pasted sequentially
        in the top dimension, returns their pasting.

        Arguments
        ---------
        fst : :class:`Diagram`
            The first diagram.
        *layers : :class:`Diagram`
            Any number of additional diagrams.

        Returns
        -------
        with_layers : :class:`Diagram`
            The pasting of all the diagrams in the top dimension.

        Raises
        ------
        :class:`ValueError`
            If the diagrams are not pastable.
        """
        utils.typecheck(fst, {'type': Diagram})
        dim = fst.dim

        diagram = fst
        for x in layers:
            utils.typecheck(x, {
                'type': Diagram,
                'st': lambda x: x.dim == dim,
                'why': 'expecting diagram of dimension {}'.format(
                    str(dim))})
            diagram, cospan = diagram.paste(x, cospan=True)

        return Diagram._new(
                diagram.shape,
                diagram.ambient,
                diagram.mapping,
                diagram.name)

    # Internal methods
    @staticmethod
    def _new(shape, ambient, mapping, name=None):
        def diagramclass():
            if isinstance(shape, shapes.Point):
                return PointDiagram
            if isinstance(shape, shapes.Arrow):
                return ArrowDiagram
            if isinstance(shape, shapes.Cube):
                return CubeDiagram
            if isinstance(shape, shapes.Simplex):
                return SimplexDiagram
            return Diagram

        new = Diagram.__new__(diagramclass())
        new._shape = shape
        new._ambient = ambient
        new._mapping = mapping
        new._name = name
        return new

    def _find_compositor(self):
        description = {
                'shape': self.shape,
                'mapping': self.mapping}
        for x in self.ambient.compositors:
            if self.ambient.compositors[x] == description:
                return x
        return None

    @staticmethod
    def _paste_fill_mapping(fst, snd, paste_cospan):
        shape = paste_cospan.target
        mapping = [
            [None for _ in n_data]
            for n_data in shape.face_data
            ]

        for x in fst.shape:
            y = paste_cospan.fst[x]
            mapping[y.dim][y.pos] = fst.mapping[x.dim][x.pos]
        for x in snd.shape:
            y = paste_cospan.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = snd.mapping[x.dim][x.pos]
            else:
                if mapping[y.dim][y.pos] != snd.mapping[x.dim][x.pos]:
                    raise ValueError(utils.value_err(
                            snd, 'boundary does not match '
                            'boundary of {}'.format(repr(fst))))
        return mapping


class SimplexDiagram(Diagram):
    """
    Subclass of :class:`Diagram` for diagrams whose shape is an
    oriented simplex.

    The methods of this class provide an implementation of the
    structural maps of a simplicial set.
    """
    def simplex_face(self, k):
        """
        Returns one of the faces of the simplex.

        This is, by definition, the pullback of the diagram along
        :code:`self.shape.simplex_face(k)`.

        Arguments
        ---------
        k : :class:`int`
            The index of the face, ranging from :code:`0` to
            :code:`self.dim`.

        Returns
        -------
        simplex_face : :class:`Diagram`
            The face.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
        face_map = self.shape.simplex_face(k)
        name = 'd[{}]({})'.format(
                str(k), str(self.name))
        return self.pullback(face_map, name)

    def simplex_degeneracy(self, k):
        """
        Returns one of the degeneracies of the simplex.

        This is, by definition, the pullback of the diagram along
        :code:`self.shape.simplex_degeneracy(k)`.

        Arguments
        ---------
        k : :class:`int`
            The index of the degeneracy, ranging from :code:`0` to
            :code:`self.dim`.

        Returns
        -------
        simplex_degeneracy : :class:`Diagram`
            The degeneracy.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
        degeneracy_map = self.shape.simplex_degeneracy(k)
        name = 's[{}]({})'.format(
                str(k), str(self.name))
        return self.pullback(degeneracy_map, name)


class CubeDiagram(Diagram):
    """
    Subclass of :class:`Diagram` for diagrams whose shape is an
    oriented cube.

    The methods of this class provide an implementation of the
    structural maps of a cubical set with connections.
    """
    def cube_face(self, k, sign):
        """
        Returns one of the faces of the cube.

        This is, by definition, the pullback of the diagram along
        :code:`self.shape.cube_face(k, sign)`.

        Arguments
        ---------
        k : :class:`int`
            Index of the face, ranging from :code:`0` to
            :code:`self.dim - 1`.
        sign : :class:`str`
            Side: :code:`'-'` or :code:`'+'`.

        Returns
        -------
        cube_face : :class:`Diagram`
            The face.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
        sign = utils.mksign(sign)
        face_map = self.shape.cube_face(k, sign)
        name = 'δ[{},{}]({})'.format(
                str(k), sign, str(self.name))
        return self.pullback(face_map, name)

    def cube_degeneracy(self, k):
        """
        Returns one of the degeneracies of the cube.

        This is, by definition, the pullback of the diagram along
        :code:`self.shape.cube_degeneracy(k)`.

        Arguments
        ---------
        k : :class:`int`
            The index of the degeneracy, ranging from :code:`0` to
            :code:`self.dim`.

        Returns
        -------
        cube_degeneracy : :class:`Diagram`
            The degeneracy.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
        degeneracy_map = self.shape.cube_degeneracy(k)
        name = 'σ[{}]({})'.format(
                str(k), str(self.name))
        return self.pullback(degeneracy_map, name)

    def cube_connection(self, k, sign):
        """
        Returns one of the connections of the cube.

        This is, by definition, the pullback of the diagram along
        :code:`self.shape.cube_connection(k, sign)`.

        Arguments
        ---------
        k : :class:`int`
            Index of the connection, ranging from :code:`0` to
            :code:`self.dim - 1`.
        sign : :class:`str`
            Side: :code:`'-'` or :code:`'+'`.

        Returns
        -------
        cube_face : :class:`Diagram`
            The connection.

        Raises
        ------
        :class:`ValueError`
            If the index is out of range.
        """
        sign = utils.mksign(sign)
        connection_map = self.shape.cube_connection(k, sign)
        name = 'γ[{},{}]({})'.format(
                str(k), sign, str(self.name))
        return self.pullback(connection_map, name)


class ArrowDiagram(SimplexDiagram, CubeDiagram):
    pass


class PointDiagram(SimplexDiagram, CubeDiagram):
    """
    Subclass of :class:`Diagram` for diagrams whose shape is a point.
    """
    def degeneracy(self, shape):
        """
        Given a shape, returns the unique degenerate diagram
        of that shape over the point.

        This is, by definition, the pullback of the point diagram along
        :code:`self.shape.terminal()`.

        Arguments
        ---------
        shape : :class:`shapes.Shape`
            The shape of the degenerate diagram.

        Returns
        -------
        degeneracy : :class:`Diagram`
            The degenerate diagram.
        """
        utils.typecheck(shape, {'type': Shape})
        return self.pullback(
                shape.terminal(), self.name)
