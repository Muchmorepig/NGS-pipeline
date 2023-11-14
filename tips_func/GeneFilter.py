import os
import click
from SimpleProcessBar import *


@click.command()
@click.argument('gene_txt', type=click.File('r'))
@click.argument('input_matrix', type=click.File('r'))
@click.argument('output_matrix', type=click.File('w'))
def gene_filter(gene_txt, input_matrix, output_matrix):
    print('Read matrix')
    gl = [line.strip() for line in gene_txt.readlines()]
    # ol = []
    print('Begin to filter')
    header = os.popen('head -n 1 {0}'.format(input_matrix.name)).read()
    # ol.append(header)
    output_matrix.write(header)
    x = SimpleProcessBar(
        int(
            os.popen('wc -l {0}'.format(
                input_matrix.name)).read().split()[0]))

    for line in input_matrix.readlines():
        x.incr()
        l = line.strip().split()[0]
        if l in gl:
            output_matrix.write(line)

    # [output_matrix.write(line) for line in ol]
    click.echo('\nResults In: %s' % (os.getcwd() + '/' + output_matrix.name))


if __name__ == '__main__':
    gene_filter()
