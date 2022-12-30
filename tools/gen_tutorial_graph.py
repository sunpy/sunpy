from pydot import Dot, Edge, Node

graph = Dot("my_graph", graph_type="graph", bgcolor="white")

core_kwargs = {'style': 'filled', 'fillcolor': '#9bdd99'}
optional_kwargs = {'style': 'filled', 'fillcolor': '#9fc5e8'}
external_kwargs = {'style': 'filled', 'fillcolor': '#efa7a7'}

# Core curriculum
graph.add_node(Node('Installing', **core_kwargs))
graph.add_node(Node('Units', **core_kwargs))
graph.add_edge(Edge('Installing', 'Units', dir='forward'))

graph.add_node(Node('Times', **core_kwargs))
graph.add_edge(Edge('Units', 'Times', dir='forward'))

graph.add_node(Node('Coordinates', **core_kwargs))
graph.add_edge(Edge('Times', 'Coordinates', dir='forward'))

graph.add_node(Node('Acquiring data', **core_kwargs))
graph.add_edge(Edge('Coordinates', 'Acquiring data', dir='forward'))

graph.add_node(Node('Maps', **core_kwargs))
graph.add_edge(Edge('Acquiring data', 'Maps', dir='forward'))

graph.add_node(Node('Timeseries', **core_kwargs))
graph.add_edge(Edge('Acquiring data', 'Timeseries', dir='forward'))

# Specific acquiring data tutorials
graph.add_node(Node('Acquiring JSOC data', **optional_kwargs))
graph.add_edge(Edge('Acquiring data', 'Acquiring JSOC data', dir='forward'))

graph.add_node(Node('Searching HEK', **optional_kwargs))
graph.add_edge(Edge('Acquiring data', 'Searching HEK', dir='forward'))

graph.add_node(Node('Working with GOES XRS', **optional_kwargs))
graph.add_edge(Edge('Timeseries', 'Working with GOES XRS', dir='forward'))

# External tutorials
graph.add_node(Node('AIAPy tutorial', **external_kwargs))
graph.add_edge(Edge('Maps', 'AIAPy tutorial', dir='forward'))

graph.add_node(Node('Querying Helioviewer.org', **external_kwargs))
graph.add_edge(Edge('Acquiring data', 'Querying Helioviewer.org', dir='forward'))

graph.write_png("tutorial_structure.png")
