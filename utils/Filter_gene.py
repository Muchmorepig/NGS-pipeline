import os
import click
from SimpleProcessBar import *


@click.command()
@click.argument('gene_txt', type=click.File('r'))
@click.argument('input_matrix', type=click.File('r'))
@click.argument('output_matrix', type=click.File('w'))

def Filter_gene(gene_txt, input_matrix, output_matrix):
    gl = [line.strip() for line in gene_txt.readlines()]
    il = [line for line in input_matrix.readlines()]
    ol = []
    x = SimpleProcessBar(len(gl))
    for gene in gl:
        x.incr()
        ol.append(
            list(
                filter(None,
                       map(lambda s: ol.append(s)
                           if gene in s else None, il))))
    final = list(filter(None, ol))
    [output_matrix.write(line) for line in final]
    click.echo('\nResults In: %s' % (os.getcwd() + '/' + output_matrix.name))


if __name__ == '__main__':
    Filter_gene()
