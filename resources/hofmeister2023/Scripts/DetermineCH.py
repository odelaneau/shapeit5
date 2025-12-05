#!/usr/bin/python3
#29-Aug-2022 Diogo Ribeiro @ UNIL

import gzip
import argparse

#===============================================================================
DESC_COMMENT = "Script to determine compound heterozygous in phased UKB exome data."
SCRIPT_NAME = "DetermineCH.py"
#===============================================================================

"""
#===============================================================================
@author: Diogo Ribeiro
@date: 29 August 2022
@copyright: Copyright 2022, University of Lausanne
"Script to determine compound heterozygous in phased UKB exome data."
#===============================================================================
"""

class DetermineCH(object):

    def __init__(self, svcfFile, mappingFile, outputFile, probCutoff):
#    def __init__(self, svcfFile, mappingFile, outputFile, nRands, probCutoff):

        self.svcfFile = svcfFile
        self.mappingFile = mappingFile
        self.outputFile = outputFile
#         self.nRands = nRands
        self.probCutoff = probCutoff


    def read_variant_gene_mapping(self):
        """
        Read file with mapping between variants and genes. Uses variant ID.
        
        Example format (TSV):
            chr1_169854517_G_GAT    ENSG00000000457 SCYL3   pLoF
            chr1_169854521_CTT_C    ENSG00000000457 SCYL3   pLoF
            chr1_169854550_TC_T     ENSG00000000457 SCYL3   pLoF        

        Only the first two columns are used.
        """
    
        print("read_variant_gene_mapping: reading %s file" % self.mappingFile)
        mappingFile = gzip.open(self.mappingFile,"rb")
        
        geneVar = {} # key -> gene, val -> set of variants
        varGene = {} # key -> var, val -> gene
        geneMapping = {} # key -> gene ENSG, val -> gene Name
        count = 0
        for line in mappingFile:
            count+=1
            if count % 10000 == 0:
                print("read_variant_gene_mapping: read %s entries.." % count)
            
            line = line.strip().decode('utf-8')
            spl = str(line).split("\t")
            
            var = spl[0]
            ensg = spl[1]
            geneName = spl[2]
            
            if ensg not in geneVar:
                geneVar[ensg] = set()
            geneVar[ensg].add(var)
            if var not in varGene:
                varGene[var] = set()
            varGene[var].add(ensg)

            geneMapping[ensg] = geneName
            
        
        print("read_variant_gene_mapping: read %i entries (%i genes, %i variants)" % (count, len(geneVar), len(varGene)) )
        self.geneVar = geneVar
        self.varGene = varGene
        self.geneMapping = geneMapping


    def read_variant_individual_mapping(self):
        """
        Read file with mapping variants and individuals in sparse VCF format.
        
        Example format (TSV):
            chr20     19475       chr20_19475_G_A   G       A       .       .       AC=9;AN=759352  1|0     3794;4435;1441;1910;518
            chr20     20032       chr_20032_G_A   G       A       .       .       AC=9;AN=759352  0|1     38;3012;1326
            chr20     23157       chr_23157_C_T   C       T       .       .       AC=10;AN=759352 1|0     2057;528;3070;5114
        
        Last column is the list of individuals harboring the variant. The column before the last is the phase of the variant.
        """

        print("read_variant_individual_mapping: reading %s file" % self.svcfFile)
        svcfFile = gzip.open(self.svcfFile,"rb")
        
        varInd = {} # key -> variant, val -> haplotype, then individual
        varProb = {} # key -> variant, val -> haplotype, then prob
        indVar = {} # key -> individual, val -> (variant,prob)
        varAF = {} # key -> variant, val -> allele frequency
        count = 0
        for line in svcfFile:
            count+=1
            if count % 10000 == 0:
                print("read_variant_individual_mapping: read %s entries" % count)
            
            line = line.strip().decode('utf-8')
            spl = line.split("\t")
            var = spl[2]
            hap = spl[9]
            inds = spl[10].split(";")
            probs = spl[11].split(";")
            info = spl[7]
            infoSpl = info.split(";")
            ac = int(infoSpl[0].split("=")[1])
            an = int(infoSpl[1].split("=")[1])
            af = ac / float(an)
            if var not in varAF:
                varAF[var] = af
        
            if var not in varInd:
                varInd[var] = {}
                varProb[var] = {}
            if hap not in varInd[var]:
                varInd[var][hap] = []
                varProb[var][hap] = []

            assert(len(inds) == len(probs)), line
          
            for i in range(len(inds)):
                ind = inds[i]
                varInd[var][hap].append(ind)
                varProb[var][hap].append(probs[i])
                if ind not in indVar:
                    indVar[ind] = {}
                if hap not in indVar[ind]:
                    indVar[ind][hap] = set()
                indVar[ind][hap].add(var+":"+probs[i])
                  
        print("read_variant_individual_mapping: read %s entries, %s variants and %s individuals" % (count, len(varInd),len(indVar)))
        self.varInd = varInd
        self.varProb = varProb
        self.indVar = indVar
        self.varAF = varAF

        # outFileInds = open(self.outputFile + ".inds", "w")
        # outFileInds.write("individual\tvariant\thaplotype\tprob\tgene\n")
        # for ind in indVar:
        #     for hap in indVar[ind]:
        #         for v in indVar[ind][hap]:
        #             var,prob = v.split(":")
        #             if var in self.varGene:
        #                 for gene in self.varGene[var]:
        #                     outFileInds.write("%s\t%s\t%s\t%s\t%s\n" % (ind, var, hap, prob, gene))
        # outFileInds.close()


    def determine_ch(self):
        """
        Determine individuals with CH for each gene.
        
        Writes output files.        
        """

        outFileGene = open(self.outputFile, "w")
        outFileGene.write("gene\tgeneName\tnumCH\tnum2mut\ttotalVars\texpectedCH\n")
        
        outFileDetail = open(self.outputFile + ".detail", "w")
        outFileDetail.write("gene\tgeneName\tindividual\tvars01\tvars10\n")

        outFileVars = open(self.outputFile + ".vars", "w")
        outFileVars.write("gene\tvar\taf\tinds10\tinds01\tprobs10\tprobs01\n")

        ## Analysis per gene
        found=0
        notFound=0
        chGenes = 0
        chEvents = 0
        count = 0
        #seenVars = set()
        for gene in sorted(set(self.geneVar)):
            count+=1
            if count % 100 == 0:
                print("determine_ch: processed %s genes.." % count)
            
            variants = self.geneVar[gene]
            inds1List = []
            inds2List = []
            for var in sorted(variants):
                #if var in seenVars: print("determine_ch: warning: %s variant seen in a previous gene" % var)
                #seenVars.add(var)
                if var in self.varInd:
                    found+=1
                    if "1|0" in self.varInd[var]: 
                        samples1 = list(self.varInd[var]["1|0"])
                        probs1 = list(self.varProb[var]["1|0"])

                        if self.probCutoff:
                            filt1 = [i for i in range(len(probs1)) if float(probs1[i]) >= self.probCutoff ]
                            if len(filt1) != len(samples1): # Apply filter
                                samples1 = [samples1[i] for i in range(len(samples1)) if i in filt1]
                                probs1 = [probs1[i] for i in range(len(probs1)) if i in filt1]

                        indText1 = ",".join(samples1)
                        probText1 = ",".join(probs1)
                        inds1List = inds1List + list(samples1)
                    else: 
                        indText1 = ""
                        probText1 = ""
                        
                    if "0|1" in self.varInd[var]: 
                        samples2 = self.varInd[var]["0|1"]
                        probs2 = list(self.varProb[var]["0|1"])

                        if self.probCutoff:
                            filt2 = [i for i in range(len(probs2)) if float(probs2[i]) >= self.probCutoff ]
                            if len(filt2) != len(samples2): # Apply filter
                                samples2 = [samples2[i] for i in range(len(samples2)) if i in filt2]
                                probs2 = [probs2[i] for i in range(len(probs2)) if i in filt2]
                                                
                        indText2 = ",".join(samples2)
                        probText2 = ",".join(probs2)
                        inds2List = inds2List + list(samples2)
                    else: 
                        indText2 = ""
                        probText2 = ""

                    outFileVars.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, var, self.varAF[var], indText1, indText2, probText1, probText2))
        #           else:
        #               print("%s only present in one haplotype" % var)
                else:
                    notFound+=1

            union = inds1List + inds2List
            
            
            expected = 0
            for ind in set(union):
                o = union.count(ind)
                if o > 1:
                    p = 1 - 0.5 ** (o-1)
                    expected+=p

            # Observed            
            overlap = set([x for x in inds1List if x in inds2List])
            
            numCH = len(overlap)
            if numCH > 0: chGenes+=1

            # Count individuals with 2 or more mutations (regardless of strand)
            inds2mut = 0
            for ind in set(union):
                v1=[]
                v2=[]
                if "1|0" in self.indVar[ind]: 
                    v1 = [x for x in self.indVar[ind]["1|0"] if x.split(":")[0] in variants]
                    if self.probCutoff: v1 = [v for v in v1 if float(v.split(":")[1]) >= self.probCutoff]
                if "0|1" in self.indVar[ind]: 
                    v2 = [x for x in self.indVar[ind]["0|1"] if x.split(":")[0] in variants]
                    if self.probCutoff: v2 = [v for v in v2 if float(v.split(":")[1]) >= self.probCutoff]
                indV = v1 + v2
                if len(indV) > 1:    inds2mut+=1

            outFileGene.write("%s\t%s\t%i\t%s\t%s\t%s\n" % (gene, self.geneMapping[gene], numCH, inds2mut, len(variants), expected))
        
            for ind in overlap:
                chEvents+=1
                v1 = []
                v2 = []
                if "1|0" in self.indVar[ind]: v1 = [x for x in self.indVar[ind]["1|0"] if x.split(":")[0] in variants]
                if "0|1" in self.indVar[ind]: v2 = [x for x in self.indVar[ind]["0|1"] if x.split(":")[0] in variants]
                outFileDetail.write("%s\t%s\t%s\t%s\t%s\n" % (gene, self.geneMapping[gene], ind, ";".join(v1), ";".join(v2))) 

        print("determine_ch: %i gene-variants were considered, %i gene-variants in the mapping were not found in the svcf file" % (found,notFound))
        print("determine_ch: processed %s genes" % count)
        print("determine_ch: %i CH events found, %i genes" % (chEvents, chGenes))
        
        outFileGene.close()
        outFileDetail.close()
        outFileVars.close()


    def run(self):
        """
        Run functions in order
        """

        self.read_variant_gene_mapping()
        
        self.read_variant_individual_mapping()

        self.determine_ch()


if __name__ == "__main__":


    print ("STARTING " + SCRIPT_NAME)

    #===============================================================================
    # Get input arguments
    #===============================================================================
    parser = argparse.ArgumentParser(description= DESC_COMMENT) 

    # positional args
    parser.add_argument('svcfFile', metavar='svcfFile', type=str,
                         help='Sparse VCF file from "vcf_to_sparse.py".')
    parser.add_argument('mappingFile', metavar='mappingFile', type=str,
                         help='Variant-gene mapping file.')
    parser.add_argument('outputFile', metavar='outputFile', type=str,
                         help='File where output will be written.')
    parser.add_argument('--probCutoff', metavar='probCutoff', type=float, default = 0,
                         help='Minimum number of haplotype posterior probability to be considered. (Default = 0, i.e. no filter)')
            
    #gets the arguments
    args = parser.parse_args( )

    #===============================================================================
    # Run analysis / processing
    #===============================================================================

    # Initialise class    
    run = DetermineCH( args.svcfFile, args.mappingFile, args.outputFile, args.probCutoff)

    run.run()

    # Stop the chrono      
    print( "FINISHED " + SCRIPT_NAME )


