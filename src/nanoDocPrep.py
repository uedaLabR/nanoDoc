import click
import SamplingPlan
import to5merpq
import os
import initialtrainingCNN
import make5merwight
import itertools

@click.group()
def cmd():
    pass


@cmd.command()
@click.option('-rraw', '--refRaw')
@click.option('-ref', '--refFastaFile')
@click.option('-ssize', '--samplesize')
@click.option('-out', '--out')
def make5mer(refRaw,refFastaFile,samplesize,out):

    click.echo('make5mer')
    #samplingplan
    indexf = refRaw+"/index.txt"
    path_w = out + "/sampleplingplan.txt"
    os.mkdir(path_w)
    SamplingPlan.makeSamplePlan(refFastaFile,indexf,samplesize,path_w)
    #make 5mer
    to5merpq.mkpq(path_w,refFastaFile,refRaw,out)

@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-out', '--outwight')
@click.option('-epochs', '--epochs',defult=500)
def trainCNN(in5mmer,outwight,epochs):

    click.echo('trainCNN')
    initialtrainingCNN.main(in5mmer,outwight,epochs)

@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-out', '--outwight')
@click.option('-wight', '--bestwight')
@click.option('-epochs', '--epochs',defult=3)
def trainDoc(in5mmer,outwight,bestwight,epochs):

    click.echo('trainDoc')
    nucs = ('A','T','C','G')
    for n1,n2,n3,n4,n5 in itertools.product(nucs, nucs, nucs,nucs, nucs):

        nuc = n1+n2+n3+n4+n5
        make5merwight.train(in5mmer,outwight,nuc,bestwight,epochs)

def main():
    cmd()


if __name__ == '__main__':
    main()