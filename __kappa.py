
#!/usr/bin/python

from sequence import Sequence

from seqfileparser import SequenceFileParser, SequenceFileParserException


## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================

if __name__=="__main__":
    import argparse 
 
    parser = argparse.ArgumentParser()

    parser.add_argument("--sequence-file","-s", help="File containing sequence protein sequence")
    parser.add_argument("--kappa-only","-k", help="Only print the sequence's kappa value",action='store_true')

    
    args = parser.parse_args()

    
    if args.sequence_file:

        parserMachine = SequenceFileParser()
        SeqFileParser = SequenceFileParser() # create a sequence file parsing object

        try:
            seq=parserMachine.parseSeqFile(args.sequence_file)
        except SequenceFileParserException:
            print "ok.."
            exit
            

        SequenceObject = Sequence(seq)

        if args.kappa_only:

            print str(len(seq)) +", "+ str(SequenceObject.Fplus()) +", " + str(SequenceObject.Fminus()) + ", "+str(SequenceObject.kappa()) + ", " + str(SequenceObject.deltaMax()) + ", " + str(SequenceObject.delta())
            
                
        else:
            print SequenceObject.toFileString()
            


