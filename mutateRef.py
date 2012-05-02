#!/usr/bin/python

import os.path
import sys
import optparse
import random

import referenceClass
from referenceClass import *

import variantsClass
from variantsClass import *

import contigs
from contigs import *

__author__ = "alistair ward"
__version__ = "version 0.02"
__date__ = "april 2012"

def main():

  # Parse the command line options
  usage = "Usage: mutateRef.py [options] -r <ref>.fa"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-r", "--reference",
                    action="store", type="string",
                    dest="refFile", help="input reference fasta")
  parser.add_option("-o", "--output",
                    action="store", type="string",
                    dest="output", help="output file name stub")
  parser.add_option("-d", "--define-event",
                    action="append", type="string",
                    dest="definedVariants", help="specified variant in the form <type>_<length>_<ref sequence>_<position>")
  parser.add_option("-l", "--length",
                    action="store", type="int",
                    dest="sequenceLength", help="length of lines in output reference file (default: 60)")
  parser.add_option("-f", "--flank-length",
                    action="store", type="int",
                    dest="flankLength", help="length of flanking sequence in alternate contigs (default: 100)")
  parser.add_option("-s", "--snps",
                    action="store", type="int",
                    dest="numberSnps", help="number of SNPs to generate (default: 0)")
  parser.add_option("-m", "--mnps",
                    action="store", type="int",
                    dest="numberMnps", help="number of MNPs to generate (default: 0)")
  parser.add_option("-e", "--maximum-mnp-length",
                    action="store", type="int",
                    dest="maximumMnpLength", help="maximum length of inserted MNPs (default: 4)")
  parser.add_option("-i", "--indels",
                    action="store", type="int",
                    dest="numberIndels", help="number of indels to generate (default: 0)")
  parser.add_option("-n", "--maximum-indel-length",
                    action="store", type="int",
                    dest="maximumIndelLength", help="maximum length of inserted indel (default: 5)")

  (options, args) = parser.parse_args()

  # Check the command line arguments and ensure that necessary options are
  # supplied and set the defaults for options not selected.
  setDefaultOptions(parser, options)

  # Define a reference and a variants object.
  ref = reference(options)
  var = variants(options)

  # Parse and evaluate any user defined variants to be included.
  parseDefinedVariants(parser, options, var)

  # Open the input fasta file.
  ref.openInputFasta(options)

  # Get information on the input reference file (e.g. the number
  # of reference sequence and the length of the sequences).
  for line in ref.filehandle:
    ref.line = line.rstrip("\n")
    ref.getSequenceInformation()

  # The number of sequences started at 0 to ensure that the lists were
  # populated in 0-based coordinates.  Add one to the number so that the
  # correct number of sequences is represented.
  ref.numberOfSequences = ref.numberOfSequences + 1

  # Close the input fastq file.
  ref.closeInputFasta()

  # If there are specific variants defined, check that they fall within
  # the reference sequences available and define the whole fasta
  # coordinates (the coordinates that run from 0 to the total length
  # of the reference fasta file).
  if len(var.definedVariantTypes) > 0: var.checkDefinedVariants(ref)

  # Randomly pick positions to mutate based on the input parameters.
  if var.numberSnps    != 0: var.createSnps(ref)
  if var.numberMnps    != 0: var.createMnps(ref)
  if var.numberIndels  != 0: var.createIndels(ref)

  # MNPs and deletions were created with each mutated base appearing in
  # the var.variantType hash in order to easily check that different
  # variants did not overlap.  Now that all variants have been created,
  # there is only a need to keep entries at the anchor base position.
  var.removeExtraEntries()

  # Mutate the reference sequence.
  # Reopen the fasta file
  ref.openInputFasta(options)

  # Open the output files.  If there are no variants, only an output fasta
  # file is required.  If there are variants, then the output file stub
  # in the options needs to be defined and an output reference including
  # all of the mutations will be produced along with a file containing all
  # the alternate contigs.
  ref.openOutputFiles(options)

  # Read in line by line, check if any mutations are required and write
  # out the new reference as well as a file containing mutation positions
  # and types.  Work through each reference sequence in turn.
  if len(var.variantType) != 0: print >> sys.stdout, "Mutating the reference."
  currentRefID = -1

  # Define a contigs object to hold all of the information about alternate
  # contigs.
  altContigs = contigs()

  # Read through the input reference file creating contigs as necessary and
  # output the final reference.
  for line in ref.filehandle:
    line = line.rstrip("\n")
    if line.startswith(">"):

      # Write out any of the sequence from the previous reference sequence before
      # building the output sequence for the new one.
      ref.writeOutputSequence(True)

      print >> ref.outputFilehandle, line
      currentRefID += 1

      # Define the contig information for the new reference sequence.
      hasVariants = altContigs.defineAltContigInformation(var, currentRefID);
      sequenceEnd = -1
    else:
      sequenceStart = sequenceEnd + 1
      sequenceEnd   = sequenceStart + len(line) - 1

      # If there are no variants in this reference sequence, just send the
      # new line to the output sequence.
      if hasVariants == False:
        ref.outputSequence += line
      else:

        # Get the position of the next variant.  If the last variant has been
        # completed, the array altContigs.sortedVariantPositions will be empty.
        if len(altContigs.sortedVariantPositions) == 0:
          ref.outputSequence += line
          hasVariants = False
        else:

          # Loop through each variant in turn and determine the course of action.
          # If the variant is beyond this segment of reference sequence, there
          # will be no variants in this segment, so it can be written to the
          # output.
          processVariants = True
          writtenToOutput = False
          variantID       = 0
          while (processVariants == True) and (variantID < len(altContigs.sortedVariantPositions)):
            if len(altContigs.sortedVariantPositions) == 0: break
            variantPosition = altContigs.sortedVariantPositions[variantID]

            # Check if alleles are under construction (i.e. the alleles span multiple
            # lines in the input fasta file).  If so, work on constructing the alleles.
            if altContigs.constructingAlleles[variantPosition] == True:
              altContigs.continueConstructingAlleles(var, ref, currentRefID, variantPosition, line, sequenceStart)
              variantID += 1

              # If this is the last variant, write the rest of this line of reference
              # sequence to the output sequence.
              if variantID >= len(altContigs.sortedVariantPositions):
                altContigs.lastVariant(var, ref, currentRefID, variantPosition, line, sequenceStart)
            else:
  
              # Check if the next variant allele starts in this line of sequence.  Also
              # check if the start of the leading flank for the next variant allele is
              # contained in this line of sequence.
              alleleInThisSequence     = True if (variantPosition <= sequenceEnd) else False
              alleleInPreviousSequence = True if (variantPosition < sequenceStart) else False
              flankInSequence          = True if ((variantPosition - var.flankLength) <= sequenceEnd) else False

              # If the next allele and the start of the leading flank for the next allele
              # are not contained in this line of the reference sequence, add the line to
              # the output sequence, but none of the contigs need to be updated.
              if (alleleInThisSequence == False) and (flankInSequence == False):
                processVariants = False

                # Write to the output sequence.
                altContigs.writeToOutputSequence(var, ref, currentRefID, variantPosition, variantID, line, sequenceStart)
  
              # If the next allele position is not in the line of sequence, but the start
              # of the leading flank is, add the required sequence to the contig.
              elif (alleleInThisSequence == False) and (flankInSequence == True):
                altContigs.sequence[variantPosition] += line[(max(0, (variantPosition - var.flankLength - sequenceStart))):len(line)]

                # Write to the output sequence.
                altContigs.writeToOutputSequence(var, ref, currentRefID, variantPosition, variantID, line, sequenceStart)
                variantID += 1

                # If this is the last variant, write the rest of this line of reference
                # sequence to the output sequence.
                if variantID >= len(altContigs.sortedVariantPositions):
                  altContigs.lastVariant(var, ref, currentRefID, variantPosition, line, sequenceStart)

              # If the allele was in the previous sequence, just add trailing flank and
              # check if the contig is complete (if the allele construction was incomplete
              # this would have been caught already).
              elif alleleInPreviousSequence == True:
                altContigs.addTrailingFlank(ref, var, currentRefID, variantPosition, line)

                # If this is the last variant, write the rest of this line of reference
                # sequence to the output sequence.
                variantID += 1
                if variantID >= len(altContigs.sortedVariantPositions):
                  altContigs.lastVariant(var, ref, currentRefID, variantPosition, line, sequenceStart)
  
              # If the next allele appears in this line of reference sequence, add the sequence
              # up to the variant position to the contig and the output sequence and then
              # determine the alleles.
              elif (alleleInThisSequence == True):
                start = max(0, (variantPosition - var.flankLength - sequenceStart))
                end   = min(len(line), (variantPosition - sequenceStart))
                altContigs.sequence[variantPosition] += line[start:end]

                # Write to the output sequence.
                altContigs.writeToOutputSequence(var, ref, currentRefID, variantPosition, variantID, line, sequenceStart)

                # Set the length of the output sequence.  It is possible that due to the position
                # of the variant, there is insufficient sequence to complete the flanking sequences
                # to the requested lengths.  Determination of a completed contig depends on
                # comparing the length of the contig to its calculated length.  In the absence of
                # sufficient leading flank, this calculated length will be smaller than expected.
                # Set the value to the length of the leading flank and the allele length etc. will
                # be added later.
                altContigs.contigLength[variantPosition] = len(altContigs.sequence[variantPosition])
                altContigs.constructAlleles(ref, var, currentRefID, variantPosition, line, sequenceStart, sequenceEnd)
                variantID += 1

                # If this was the last variant in this reference sequence, the loop will be
                # terminated and so the remainder of the sequence contained in line needs to be
                # written to the output sequence.
                if variantID >= len(altContigs.sortedVariantPositions):
                  altContigs.lastVariant(var, ref, currentRefID, variantPosition, line, sequenceStart)

          # Each of the variants has been evaluated for this line of the input reference sequence.  If
          # any of the contigs were completed, purge them from the data structures.
          i = 0;
          while i < len(altContigs.sortedVariantPositions):
            if altContigs.completed[altContigs.sortedVariantPositions[i]] == True:
              altContigs.clearContigs(var, currentRefID, altContigs.sortedVariantPositions[i])
            else:
              i += 1

          # If all of the variants have been removed, set hasVariants to False.
          if len(altContigs.sortedVariantPositions) == 0: hasVariants = False

      # Write the output reference sequence to file.
      ref.writeOutputSequence(False)

  # Write out any remaining sequence.
  ref.writeOutputSequence(True)

  # Close the files.
  ref.closeInputFasta()
  ref.closeOutputFiles()

  # End the program
  return 0

def setDefaultOptions(parser, options):

  # Check that a single reference fastq file is given.
  if options.refFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput reference fasta file (--reference, -r) is required."
    exit(1)

  # Set the number of inserted variants to zero for variant types not specified in the
  # command line.
  options.numberSnps    = 0 if options.numberSnps == None else options.numberSnps
  options.numberMnps    = 0 if options.numberMnps == None else options.numberMnps
  options.numberIndels  = 0 if options.numberIndels == None else options.numberIndels

  # If a sequence length for each line in the fastq file is not
  # defined, default to 60.
  options.sequenceLength = 60 if options.sequenceLength == None else options.sequenceLength

  # If a flank length is not defined, use the default of 100.
  options.flankLength = 100 if options.flankLength == None else options.flankLength

  # If MNPs are to be included and the maximum MNP length isn't set, default
  # to a maximum length of 4.
  options.maximumMnpLength = 4 if options.maximumMnpLength == None else options.maximumMnpLength

  # If indels are to be included and the maximum indel length isn't set, default
  # to a maximum length of 5.
  options.maximumIndelLength = 5 if options.maximumIndelLength == None else options.maximumIndelLength

def parseDefinedVariants(parser, options, var):

  # Check if any variants have been defined.
  if options.definedVariants != None:
    for line in options.definedVariants:

      # There should be three underscores in the string.  If not,
      # throw an exception.
      if line.count("_") != 3:
        print >> sys.stderr, "ERROR: Unexpected format for specified variant."
        print >> sys.stderr, "  Expected: <snp/mnp/ins/del>_<length>_<ref sequence ID>_<position>"
        print >> sys.stderr, "  Got:", line
        exit(1)

      # Break the string apart and determine the variant, type, length,
      # reference sequence ID and position.
      #
      # First, deteremine the variant type.
      variantType = line[0:line.index("_")]
      if (variantType != "snp") & (variantType != "mnp") & (variantType != "ins") & (variantType != "del"):
        print >> sys.stderr, "ERROR: unknown variant type"
        print >> sys.stderr, "  Allowed types: snp, mnp, ins or del"
        print >> sys.stderr, "  Got:", variantType
        exit(1)
      var.definedVariantTypes.append(variantType)

      # Now determine the variant length.
      line = line[line.index("_") + 1:len(line)]
      variantLength = int(line[0:line.index("_")])
      var.definedVariantLengths.append(variantLength)

      # If the variant type is a SNP, check that the length is 1.  If it is
      # an MNP, check that the length is greater than 1.
      if (variantType == "snp") & (variantLength != 1):
        print >> sys.stderr, "ERROR: A SNP is defined with length not equal to 1 ( length =", variantLength, ")"
        exit(1)
      if (variantType == "mnp") & (variantLength < 2):
        print >> sys.stderr, "ERROR: An MNP is defined with length less than 2."
        exit(1)

      # Now determine the reference sequence ID.
      line = line[line.index("_") + 1:len(line)]
      refID = int(line[0:line.index("_")]) - 1
      var.definedVariantRefIDs.append(refID)

      # Finally, determine the position.
      line = line[line.index("_") + 1:len(line)]
      position = int(line[0:len(line)]) - 1
      var.definedVariantPositions.append(position)

if __name__ == "__main__":
  main()
