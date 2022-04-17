import networkx as nx
import discopy

import rewal
from rewal import utils
from rewal.strdiags import StrDiag


def to_discopy(diagram):
    utils.typecheck(diagram, {
        'type': rewal.diagrams.Diagram,
        'st': lambda x: 0 <= x.dim <= 2,
        'why': 'expecting a diagram of dimension 2 or less'})

    if diagram.dim == 0:
        return discopy.monoidal.Id(
                discopy.monoidal.Ty())

    # Get input and output
    strinput = StrDiag(diagram.boundary('-', 1))
    input_sort = [x for x in nx.topological_sort(strinput.graph)
                  if x in strinput.nodes]
    inputs = discopy.monoidal.Ty()
    for x in input_sort:
        if not strinput.nodes[x]['isdegenerate']:
            inputs = inputs @ discopy.monoidal.Ty(
                    strinput.nodes[x]['label'])

    if diagram.dim == 1:
        return discopy.monoidal.Id(inputs)

    stroutput = StrDiag(diagram.boundary('+', 1))
    output_sort = [x for x in nx.topological_sort(stroutput.graph)
                   if x in stroutput.nodes]
    outputs = discopy.monoidal.Ty()
    for x in output_sort:
        if not stroutput.nodes[x]['isdegenerate']:
            outputs = outputs @ discopy.monoidal.Ty(
                    stroutput.nodes[x]['label'])

    # Produce a good layering of the diagram
    prelayers = diagram.layers
    layers = []
    for layer in prelayers:
        if layer.dim < 2:
            continue
        if layer.shape.size[2] > 1:
            layer.generate_layering()
            for x in layer.layers:
                layers.append(x)
        else:
            layers.append(layer)

    boxes = []
    offsets = []
    for layer in layers:
        strlayer = StrDiag(layer)
        for x in strlayer.nodes:
            node = x
        if strlayer.nodes[node]['isdegenerate']:
            continue

        total_graph = strlayer.widthgraph
        total_graph.add_edges_from(strlayer.graph.edges)
        sort = list(nx.topological_sort(total_graph))
        in_faces = layer.shape.faces(node, '-')
        out_faces = layer.shape.faces(node, '+')

        offset, idx = 0, 0
        while sort[idx] not in in_faces:
            if not strlayer.wires[sort[idx]]['isdegenerate']:
                offset += 1
            idx += 1
        offsets.append(offset)

        input_range = range(idx, idx+len(in_faces))
        output_range = range(
                idx+len(in_faces)+1,
                idx+len(in_faces)+1+len(out_faces))

        def extend_ty(ty, wire):
            if not strlayer.wires[wire]['isdegenerate']:
                return ty @ discopy.monoidal.Ty(
                        strlayer.wires[wire]['label'])
            return ty

        box_in = discopy.monoidal.Ty()
        for k in input_range:
            box_in = extend_ty(box_in, sort[k])
        box_out = discopy.monoidal.Ty()
        for k in output_range:
            box_out = extend_ty(box_out, sort[k])

        box = discopy.monoidal.Box(
                discopy.monoidal.Ty(strlayer.nodes[node]['label']),
                box_in,
                box_out)
        boxes.append(box)

    return discopy.monoidal.Diagram(
            inputs, outputs, boxes, offsets)
