# Author: Aaron Graubert  https://github.com/agraubert
import nbformat as nbf
import base64
import io
import json
import sys

def trim(docstring):
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return '\n'.join(trimmed)

def encode_figure(figure, **kwargs):
    img = io.BytesIO()
    figure.savefig(img, **kwargs)
    img.seek(0,0)
    return nbf.v4.new_output(
        'display_data',
        {
            'text/plain': [repr(figure)],
            'image/png': base64.b64encode(img.read()).decode()
        }
    )

def encode_dataframe(df, n, **kwargs):
    return nbf.v4.new_output(
        'execute_result',
        {
            'text/plain': [df.to_string()],
            'text/html': [df.to_html(**kwargs)]
        },
        execution_count=n
    )

def encode_output(obj, n):

    return nbf.v4.new_output(
        'execute_result',
        {'text/plain': [repr(obj)]},
        execution_count=n
    ) if obj is not None else None

class Notebook(object):
    """
    Wrapper to nbformat Notebook
    """
    def __init__(self, header=None):
        self.nb = nbf.v4.new_notebook()
        if header is not None:
            self.add_markdown_cell(header, '---', 'Created by the nb_encode api')
        self.exec_count = 1

    def add_markdown_cell(self, *lines):
        lines = [line.rstrip()+'\n' for line in lines]
        lines[-1] = lines[-1][:-1]
        self.nb['cells'].append(nbf.v4.new_markdown_cell(lines))

    def add_code_cell(self, source, *outputs, **kwargs):
        if isinstance(source, list):
            source = '\n'.join(line.rstrip() for line in source)
        self.nb['cells'].append(nbf.v4.new_code_cell(
            source,
            execution_count=self.exec_count,
            outputs=[
                encode_output(output, self.exec_count)
                if not isinstance(output, nbf.notebooknode.NotebookNode)
                else output
                for output in outputs
                if output is not None
            ],
            **kwargs
        ))
        self.exec_count += 1

    def write(self, dest):
        if isinstance(dest, str):
            with open(dest, 'w') as w:
                nbf.write(self.nb, w)
        else:
            nbf.write(self.nb, dest)

def encode_plot_cell(cell, source, result, figure):
    img = io.BytesIO()
    figure.savefig(img)
    img.seek(0,0)
    img = base64.b64encode(img.read())
    output_cell = nbf.v4.new_code_cell(
        source,
        outputs=[
            nbf.v4.new_output(
                'execute_result',
                {
                    'text/plain': [result]
                },
                execution_count=cell
            ),
            nbf.v4.new_output(
                'display_data',
                {
                    'text/plain': [repr(figure)],
                    'image/png': img.decode()
                }
            )
        ]
    )
    return output_cell

def encode_standard_cell(cell):
    source = eval('_i%d'%cell)
    try:
        result = repr(eval('_%d'%cell))
    except:
        result = None
    output_cell = nbf.v4.new_code_cell(
        source,
        outputs=([
            nbf.v4.new_output(
                'execute_result',
                {'text/plain': [result]},
                execution_count=cell
            )
        ] if result is not None else [])
    )
    return output_cell
