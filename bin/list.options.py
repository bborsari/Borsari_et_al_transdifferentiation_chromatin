#!/usr/bin/env python


import sys, ast, getopt, types


def main(argv):
    arg_dict={}
    switches=['mark', 'sample', 'db', 'Loci', 'Method', 'log', 'window', 'type', 'rep1', 'Rep2', 'outFolder']
    singles=''.join([x[0]+':' for x in switches])
    long_form=[x+'=' for x in switches]
    d={x[0]+':':'--'+x for x in switches}

    try:            
        opts, args = getopt.getopt(argv, singles, long_form)
    except getopt.GetoptError:          
        print "Error: bad arg"                       
        sys.exit(2)

    for opt, arg in opts:
        if opt[1]+':' in d: o=d[opt[1]+':'][2:]
        elif opt in d.values(): o=opt[2:]
        else: o =''
        if o and arg:
            arg_dict[o]=ast.literal_eval(arg)
        if not o:    
            print opt, arg, "Error: bad arg"
            sys.exit(2)
                 
    if 'type' not in arg_dict:
	arg_dict['type']='pileupSignal'

    for m in arg_dict['mark']:
	for s in arg_dict['sample']:
		if 'rep1' not in arg_dict:
			rep1="".join((s,m,"X1"))
		else:
			rep1=arg_dict['rep1'][(s, m)]

                if 'Rep2' not in arg_dict:
                        Rep2="".join((s,m,"X2"))
		else:
			Rep2=arg_dict['Rep2'][(s, m)]
			
		for d in arg_dict['db']:
			for j in range(0, len(arg_dict['window'])):
				for M in arg_dict['Method']:
					for l in arg_dict['log']:
						out="".join((arg_dict['outFolder'][0], "/", s, m, ".", d.split("/")[-1], ".", arg_dict['window'][j], ".", M, ".", l, ".correlation.tsv"))
						options="\t".join((m, s, d, arg_dict['Loci'][j], M, l, arg_dict['window'][j], arg_dict['type'], rep1, Rep2, out))
						print options

if __name__ == '__main__':
    main(sys.argv[1:])        
