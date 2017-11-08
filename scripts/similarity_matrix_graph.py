import numpy as np
import json
import os
import re
from bokeh.plotting import figure, save, output_file
from bokeh.models import HoverTool, ColumnDataSource, LinearColorMapper


def generate_similarity_graph(sim_input_json_file, sim_out_html_file, title="similarity matrix"):
    data = json.load(open(sim_input_json_file))
    records = data['records']
    sim_array = data['similarities']
    return _generate_similarity_graph(records, sim_array, sim_out_html_file, title=title)


def generate_similarity_graph_go(sim_input_json_file, sim_out_html_file, title="similarity matrix"):
    data = json.load(open(sim_input_json_file))

    # get feature names and sort them
    names = set(map(lambda d: d['product1']['name'], data))
    names.update(set(map(lambda d: d['product2']['name'], data)))
    names = sorted(names)

    # generate records and similarities arrays
    records = map(lambda n: {'name': n}, names)
    sim_array = map(lambda d: {'val': d['similarity'], 'x': names.index(d['product1']['name']),
                               'y': names.index(d['product2']['name'])}, data)

    return _generate_similarity_graph(records, sim_array, sim_out_html_file, title=title)


def _generate_similarity_graph(records, sim_array, sim_out_html_file, title="similarity matrix"):
    N = len(records)
    similarities = np.zeros((N, N))
    max_sim = 0
    for val in sim_array:
        sim = val['val'] if val['x'] != val['y'] else 0
        if sim > max_sim:
            max_sim = sim
        similarities[val['x']][val['y']] = similarities[val['y']][val['x']] = sim

    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    mapper = LinearColorMapper(palette=colors, low=0.0, high=max_sim)

    names = [str(i) for i in range(0, N)]

    xname = []
    yname = []
    xlabel = []
    ylabel = []
    for i in range(0, N):
        for j in range(0, N):
            xname.append(str(i))
            xlabel.append(records[i]['name'])
            yname.append(str(j))
            ylabel.append(records[j]['name'])

    source = ColumnDataSource(data=dict(
        xname=xname,
        yname=yname,
        xlabel=xlabel,
        ylabel=ylabel,
        similarity=similarities.flatten(),
    ))

    p = figure(title=title,
               x_axis_location="above", tools="hover,save",
               x_range=list(reversed(names)), y_range=names)

    p.plot_width = 800
    p.plot_height = 800
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = np.pi / 3

    p.rect('xname', 'yname', 0.9, 0.9, source=source, line_color=None, hover_line_color='black',
           fill_color={'field': 'similarity', 'transform': mapper})

    p.select_one(HoverTool).tooltips = [
        ('x', '@xlabel'),
        ('y', '@ylabel'),
        ('sim', '@similarity'),
    ]

    output_file(sim_out_html_file, title=title)
    save(p)  # save the plot


if __name__ == "__main__":
    base_path = '/home/victor/Escritorio/Genotipado_Alternativo/colocalizacion/go_distances'

    file_pattern = 'go_similarities_(.+)\.json'
    prog = re.compile(file_pattern)

    sim_files = [f for f in os.listdir(base_path) if os.path.isfile(os.path.join(base_path, f)) and prog.match(f)]

    for f in sim_files:
        sim_input_file = os.path.join(base_path, f)
        suffix = prog.match(f).groups()[0]
        sim_out_file = os.path.join(base_path, "sim_matrix_%s.html" % suffix)
        generate_similarity_graph_go(sim_input_file, sim_out_file, title="Similarity matrix %s" % suffix)

    # generate_similarity_graph(sim_input_file, sim_out_file)
