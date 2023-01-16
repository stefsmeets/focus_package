#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import input
from builtins import range
from builtins import object
import sys
import os
import subprocess as sp
import time
import argparse
from multiprocessing import cpu_count

from .__init__ import __version__


class ProgressBar(object):
    """docstring for ProgressBar"""
    def __init__(self):
        
        self.reverse = False
        self.a = ' '
        self.b = '='
        self.width = 10 # length of cycle, default is 2 seconds


    def loop(self):
        for i in range(self.width):
            if self.reverse:
                i = self.width-i
            time.sleep(0.1)
            printer('['+self.b*(i)+self.a*(self.width-i)+']')
        self.reverse = not self.reverse
        for i in range(self.width):
            if self.reverse:
                i = self.width-i
            time.sleep(0.1)
            printer('['+self.a*(self.width-i)+self.b*i+']')
        self.a,self.b = self.b,self.a
        
    def clear(self):
        printer('')


def resize_terminal(height=40,width=80):
    print("\x1b[8;{};{};t".format(height,width))


def printer(data):
    """Print things to stdout on one line dynamically"""
    sys.stdout.write("\r\x1b[K"+data.__str__())
    sys.stdout.flush()


def check_file(fn):
    print('Checking', fn, '...', end=' ')
    try:
        f = open(fn)
    except IOError as e:
        print(e)
        sys.exit()
    else:
        print('OK')
        f.close()


def run(options, filenames):
    
    n_procs = options.procs
    n_cores = cpu_count()

    # work-around for argparse not being able to parse '-' in options
    focus_args = ' '.join(options.focus_args).replace('+','-') 
    
    assert options.runs >= 0, 'Expected a positive number of runs'
    n_runs = str(options.runs)

    assert n_procs >= 0, 'Expected a positive number of processors'
    if n_procs*len(filenames) > n_cores:
        print('{} cores detected, are you sure you want to run {}*{} processes? [y/n]'.format(n_cores,len(filenames),n_procs))
        confirm = input(' >> [y] ')
        if 'n' in confirm:
            sys.exit()

    processes = []

    for fn in filenames:
        check_file(fn)

    print('Starting processes...')
    for i in range(n_procs):
        for fn in filenames:
            #fnout = fn.replace('.inp','_{}.out'.format(i))
            fnout = os.path.splitext(fn)[0] + '_{}.out'.format(i)
                
            cmd = ['focus', fn, n_runs]
            print('     >> focus', fn, n_runs, focus_args, '...', end=' ')

            if focus_args:
                cmd.append(focus_args)

            if not options.dry_run:
                p = sp.Popen(cmd,shell=False,stdout=open(fnout,'wb'))
                processes.append(p)
                print('started')
            else:
                print('not started')
            
        if n_procs != 1:
            time.sleep(1) # prevent equal random seeds
        
    if options.resize: # resize window to fit output
        resize_terminal(height=5+(i+1)*len(filenames))
    
    pb = ProgressBar()
    while any(p.poll() == None for p in processes):
        pb.loop()
    pb.clear()


def main():

    usage = """"""


    description = """
Program for running several focus instances in parallel. It will start Np 
processes for every file specified.

""" 
    
    epilog = 'Version: {}'.format(__version__)
    
    parser = argparse.ArgumentParser(#usage=usage,
                                    description=description,
                                    epilog=epilog, 
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    version=__version__)
    
    
    parser.add_argument("args", 
                        type=str, metavar="FILE",nargs='*',
                        help="Paths to input files.")


    parser.add_argument("-p", "--procs", metavar='Np',
                        action="store", type=int, dest="procs",
                        help="Number of processes to use (default = 1)")

    parser.add_argument("-r", "--runs", metavar='Nr',
                        action="store", type=int, dest="runs",
                        help="Number of runs per process (default = 200)")

    parser.add_argument("-a", "--args", metavar='arg',
                        action="store", type=str, nargs='+', dest="focus_args",
                        help="Arguments to pass directly to focus. Due to limitations of the option parser, please use '+' instead of '-' to prefix every command.")
    
    parser.add_argument("-d", "--dry",
                        action="store_true", dest="dry_run",
                        help="Runs the program, but doesn't start any processes.")

    
    parser.set_defaults(procs=1,
                        runs=200,
                        focus_args="",
                        dry_run=False,
                        resize=False,
                        )
    
    options = parser.parse_args()
    args = options.args

    if not args:
        parser.print_help()
        sys.exit()

    run(options,args)


if __name__ == '__main__':
    main()