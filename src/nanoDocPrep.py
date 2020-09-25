import click
import SamplingPlan
import to5merpq
import os
import initialtrainingCNN
import make5merwight
import itertools
import os

@click.group()
def cmd():
    pass


@cmd.command()
@click.option('-rraw', '--refraw')
@click.option('-r', '--ref')
@click.option('-ssize', '--samplesize',default=12000)
@click.option('-o', '--out')
def make5mer(refraw,ref,samplesize,out):

    click.echo('make5mer')
    #samplingplan
    indexf = refraw+"/index.txt"
    path_w = out + "/sampleplingplan.txt"
    if not os.path.exists(path_w):
        SamplingPlan.makeSamplePlan(ref,indexf,samplesize,path_w)
    #make 5mer
    to5merpq.mkpq(path_w,ref,refraw,out)

@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-o', '--outwight')
@click.option('-ssize', '--samplesize',default=1200)
@click.option('-epochs', '--epochs',default=500)
def traincnn(in5mmer,outwight,samplesize,epochs):

    click.echo('trainCNN')
    initialtrainingCNN.main(in5mmer,outwight,samplesize,epochs)

@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-o', '--outwight')
@click.option('-wight', '--bestwight')
@click.option('-ssize', '--samplesize',default=12000)
@click.option('-epochs', '--epochs',default=3)
def traindoc(in5mmer,outwight,bestwight,samplesize,epochs):

    click.echo('trainDoc')
    nucs = ('A','T','C','G')
    for n1,n2,n3,n4,n5 in itertools.product(nucs, nucs, nucs,nucs, nucs):

        nuc = n1+n2+n3+n4+n5
        print('training doc',nuc)
        make5merwight.train(in5mmer,outwight,nuc,bestwight,samplesize,epochs)

def main():
    cmd()


if __name__ == '__main__':
    main()