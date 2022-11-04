import argparse

from lsl.misc.wisdom import show, make

parser = argparse.ArgumentParser(
    description='build LSL-specific FFTW wisdom files', 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
parser.add_argument('-s', '--show', action='store_true', 
                    help='show information about the avaliable LSL wisdom')
args = parser.parse_args()

if args.show:
    show()
else:
    make()
