import click
import os
import h5tobinnedPq
import nanoDocAnalysis
import makeIndex
import nanoDocAnalysis
import pathlib

@click.group()
def cmd():
    pass



@cmd.command()
@click.option('-i', '--input')
@click.option('-o', '--output')
@click.option('-r', '--ref')
@click.option('-t', '--thread', default=4)
def formatFile(input,output,ref,thread):

    print(input,output)
    click.echo('make index')
    #make our put dir
    if not os.path.exists(output):
        os.mkdir(output)
    indexf = output +"/" + "index.txt"

    if not os.path.exists(indexf):
        makeIndex.mkIdx(input,indexf,thread)
    #
    click.echo('make parquet')

    h5tobinnedPq.makeParquet(indexf,output,ref,thread)


@cmd.command()
@click.option('-w', '--wight')
@click.option('-p', '--param')
@click.option('-r', '--ref')
@click.option('-rraw', '--refraw')
@click.option('-traw', '--tgraw')
@click.option('-o', '--output')
@click.option('-chrom', '--chrom',default="")
@click.option('-s', '--start',default=1)
@click.option('-e', '--end',default=-1)
@click.option('-st', '--strand',default="+")
@click.option('-minreadlen', '--minreadlen',default=200)
def analysis(wight,param,ref,refraw,tgraw,output,chrom,start,end,strand,minreadlen):

    click.echo('modification call')
    p_sub = pathlib.Path(output)
    if not os.path.exists(p_sub.parent):
        os.mkdir(p_sub.parent)
    print(refraw)
    print(tgraw)
    nanoDocAnalysis.modCall(wight,param, ref, refraw,tgraw, output, chrom, chrom, start, end, strand, minreadlen)


def main():
    cmd()


if __name__ == '__main__':
    main()