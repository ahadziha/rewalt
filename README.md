# rewalt

1. *(archaic)* to overturn, throw down
2. a library for **rew**riting, **al**gebra, and **t**opology, developed in Tallinn (aka **Reval**)

![CoolDiagram](https://raw.githubusercontent.com/ahadziha/rewalt/main/docs/_static/img/readme_1.png)

## About

`rewalt` is a toolkit for **higher-dimensional diagram rewriting**, with applications in

- **higher** and **monoidal category theory**,
- **homotopical algebra**,
- **combinatorial topology**,

and more. Thanks to its visualisation features, it can also be used as a structure-aware **string diagram** editor, supporting [TikZ](//tikz.net/) output so the string diagrams can be directly embedded in your LaTeX files.

![BlackWhiteDiagram](https://raw.githubusercontent.com/ahadziha/rewalt/main/docs/_static/img/readme_2.png)

It implements [diagrammatic sets](//arxiv.org/abs/2007.14505), which, by the "higher-dimensional rewriting" paradigm, double as a model of

- *higher-dimensional rewrite systems*, and of
- *directed cell complexes*.

This model is "topologically sound": a diagrammatic set built in `rewalt` presents a finite CW complex, and a diagram constructed in the diagrammatic set presents a valid homotopy in this CW complex.

A diagrammatic set can be seen as a generalisation of a *simplicial set* or of a *cubical set* with many more "cell shapes". As a result, `rewalt` also contains a *full implementation* of finitely presented **simplicial sets** and **cubical sets with connections**.

## Getting started

`rewalt` is available for Python 3.7 and higher. You can install it with the command

```shell
pip install rewalt
```

Then you should take a look at the [documentation](//rewalt.readthedocs.io/), which includes several worked examples from category theory, algebra, and homotopy theory.

## Usage

The [docs/notebooks/](//github.com/ahadziha/rewalt/tree/main/docs/notebooks/) folder contains several worked examples in the form of Jupyter notebooks.

For example, this is how you create a single-sorted algebraic signature with one binary operation $m$ and one constant $u$, then represent the term $m(u, -)$ as a string diagram oriented bottom-to-top.

```python
from rewalt import DiagSet

X = DiagSet()
pt = X.add('pt', draw_label=False)
a = X.add('a', pt, pt, draw_label=False)  # the sort

m = X.add('m', a.paste(a), a)  # binary operation
u = X.add('u', pt.unit(), a)  # constant

m.to_inputs(0, u).draw()
```

![ExampleAlgebra](https://raw.githubusercontent.com/ahadziha/rewalt/main/docs/_static/img/readme_3.png)

This is how you construct a 3-dimensional diagram shape as an "oriented cylinder" whose bases are 2-simplices, then output its *oriented face poset* in the form of a Hasse diagram with *magenta* edges for *input* faces, and *blue* edges for *output* faces.

```python
from rewalt import Shape

twosimplex = Shape.simplex(2)
arrow = Shape.arrow()

cylinder = arrow * twosimplex  # Gray product of arrow and 2-simplex
cylinder.hasse(labels=False)
```

![ExampleTopology](https://raw.githubusercontent.com/ahadziha/rewalt/main/docs/_static/img/readme_4.png)

## Testing

You can run all tests with the command

```shell
pytest
```

## Documentation

The latest documentation is hosted on [Read the Docs](//rewalt.readthedocs.io/).

If you want to build a local copy of the documentation, first install the required dependencies:

```shell
pip install -r docs/requirements.txt
```

Then run

```shell
cd docs/
make clean
make html
```

You will then find the documentation under `docs/_build/`.

## Further reading

For a first introduction to the ideas of higher-dimensional rewriting, diagrammatic sets, and "topological soundness", you may want to watch these presentations at the [CIRM meeting on *Higher Structures*](//cirmbox.cirm-math.fr/s/8a8DXyFA4bzaSNF) and at the [GETCO 2022 conference](//youtu.be/UlVZPiJ87kw).

A nice overview of the general landscape of higher-dimensional rewriting is Yves Guiraud's [mémoire d'habilitation](//webusers.imj-prg.fr/~yves.guiraud/articles/hdr.pdf).

So far there are two papers on the theory of diagrammatic sets: [the first one](//arxiv.org/abs/2007.14505) containing the foundations, [the second one](//arxiv.org/abs/2101.10361) containing some developments applied to categorical universal algebra.

A description and complexity analysis of some of the data structures and algorithms behind `rewalt` will be published in the [proceedings of ACT 2022](https://msp.cis.strath.ac.uk/act2022/).

## License

`rewalt` is distributed under the BSD 3-clause license; see [`LICENSE`](//github.com/ahadziha/rewalt/tree/main/LICENSE).

## Contributing

Currently, the only active developer of `rewalt` is [Amar Hadzihasanovic](//ioc.ee/~amar).

Contributions are welcome. Please reach out either by sending me an email, or by [opening an issue](https://github.com/ahadziha/rewalt/issues/new).
