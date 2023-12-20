import sys, getopt, os

'''
###############################################################################
#   General Information
###############################################################################


'''

### Global variables


### Functions


def cr2res_do():

    return 

'''
###############################################################################
#   Main
###############################################################################
'''
def main(argv):
    # Parse Inputs
    inputfile = ''
    opts, args = getopt.getopt(argv,"hi:",["ifile="])
    for opt, arg in opts:
        if opt == '-h':
            print('cr2res_gen_eop_table.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg

    # Check if input file is provided and exists
    if not os.path.isfile(inputfile):
        print("Please provide an existing input FITS file")
        sys.exit()

    # Opening FITS file
    f = open(inputfile)






    # Closing file
    f.close()


if __name__ == "__main__":
   main(sys.argv[1:])



