#!/usr/bin/python

from __future__ import print_function
import os.path
import sys
import re

class contigs:

# Constructor.
  def __init__(self):
    self.name                  = {}
    self.sequence              = {}
    self.referenceAllele       = {}
    self.alternateAllele       = {}
    self.referenceAlleleLength = {}
    self.contigLength          = {}
    self.completed             = {}
    self.constructingAlleles   = {}

# Loop over all of the variants for this reference sequence
# and set up all the information required.
  def defineAltContigInformation(self, var, refID):
    if var.variantType.has_key(refID):
      for variantPosition in var.variantType[refID].keys():
        self.name[variantPosition]                  = ""
        self.sequence[variantPosition]              = ""
        self.referenceAllele[variantPosition]       = ""
        self.alternateAllele[variantPosition]       = ""
        self.referenceAlleleLength[variantPosition] = ""
        self.contigLength[variantPosition]          = 0
        self.completed[variantPosition]             = False
        self.constructingAlleles[variantPosition]   = False

      self.sortedVariantPositions = sorted(var.variantType[refID].keys())

    return var.variantType.has_key(refID)

# Build the reference and alternate alleles.
  def constructAlleles(self, ref, var, refID, position, sequence, sequenceStart, sequenceEnd):

    # Determine the length of the alternate contig.  This will be used
    # for checking if the contig is complete.  Also determine the length
    # of the reference allele.
    if var.variantType[refID][position][0] == "del": alternateAlleleLength = 1
    else: alternateAlleleLength = var.variantType[refID][position][1]
    self.contigLength[position] += alternateAlleleLength + var.flankLength

    if var.variantType[refID][position][0] == "ins": self.referenceAlleleLength[position] = 1
    else: self.referenceAlleleLength[position] = var.variantType[refID][position][1]

    # Check if the entirety of the reference allele falls within this
    # line of reference sequence (a deletion, for instance, may have
    # deleted bases present in the next line of reference sequence.
    if (position + self.referenceAlleleLength[position] - 1) <= sequenceEnd:
      self.constructingAlleles[position] = False
      self.referenceAllele[position]     = sequence[(position - sequenceStart):(position + self.referenceAlleleLength[position] - sequenceStart)]
      self.alternateAllele[position]     = var.generateAltAllele(refID, position, self.referenceAllele[position])
      self.allelesComplete(var, ref, refID, position, sequence, sequenceStart)
    else:
      self.constructingAlleles[position] = True
      self.referenceAllele[position]     = sequence[(position - sequenceStart):len(sequence)]

# If the alleles straddled lines in the reference sequence. continue constructing
# the alleles.
  def continueConstructingAlleles(self, var, ref, refID, position, sequence, sequenceStart):

    # Update the reference allele.
    self.referenceAllele[position] += sequence[0:min((self.referenceAlleleLength[position] - len(self.referenceAllele[position])),len(sequence))]

    # If the allele is complete, calculate the alternate allele and upate the contig.
    if self.referenceAlleleLength[position] == len(self.referenceAllele[position]):
      self.constructingAlleles[position] = False
      self.alternateAllele[position]     = var.generateAltAllele(refID, position, self.referenceAllele[position])
      self.allelesComplete(var, ref, refID, position, sequence, sequenceStart)

# When the alleles are constructed, update the contig, output sequence and write
# the variants to screen.
  def allelesComplete(self, var, ref, refID, position, sequence, sequenceStart):

    # Add the alternate allele to the output reference sequence and
    # to the alternate contig and print the alleles to screen.
    ref.outputSequence += self.alternateAllele[position]
    self.sequence[position] += self.alternateAllele[position]
    print("Generated variant at: ", ref.sequenceName[refID], ":", position + 1, " (", sep="", end="", file=sys.stdout)
    print(self.referenceAllele[position], "->", self.alternateAllele[position], ")", sep="", file=sys.stdout)

    # Now add trailing sequence to the contig.  If the whole flank has
    # the contig is complete and can be written out.
    lastElement = min(len(sequence), (position + self.referenceAlleleLength[position] + var.flankLength - sequenceStart))
    self.sequence[position] += sequence[(position + self.referenceAlleleLength[position] - sequenceStart):lastElement]
    if len(self.sequence[position]) == self.contigLength[position]: self.completeContig(var, ref, refID, position)

# When a contig has been complete, write the contig to the otuput contig file
# and delete the contig from the sortedVariantPositions list as well as the
# var.variantType object.
  def completeContig(self, var, ref, refID, position):

    # Define the contig name.
    self.name[position] = '>alternate_contig:' + ref.sequenceName[refID] + "_" + str(position + 1) + "_" + self.referenceAllele[position]
    self.name[position] += "_" + self.alternateAllele[position] + ' dna:alternate_contig'
    print(self.name[position], file=ref.contigFilehandle)
    print(self.sequence[position], file=ref.contigFilehandle)
    self.completed[position] = True

# If the last variant was reached, output the remainder of the line of input
# reference sequence to the output.
  def lastVariant(self, var, ref, refID, position, sequence, sequenceStart):
    if var.variantType[refID][position][0] == "ins": variantLength = 1
    else: variantLength = var.variantType[refID][position][1]
    start = position - sequenceStart + variantLength      
    ref.outputSequence += sequence[max(0, start):len(sequence)]

# For contigs requiring more trailing flank, add the flank and check for
# completion.
  def addTrailingFlank(self, ref, var, refID, position, sequence):

    # Determine the amount of flank still required.
    numberMissingBases = self.contigLength[position] - len(self.sequence[position])
    self.sequence[position] += sequence[0:min((numberMissingBases), len(sequence))]
    if len(self.sequence[position]) == self.contigLength[position]: self.completeContig(var, ref, refID, position)

# Clear the variants and contig information for a particular.
# reference sequence.
  def clearContigs(self, var, refID, position):
    del(var.variantType[refID][position])
    del(self.name[position])
    del(self.sequence[position])
    del(self.referenceAllele[position])
    del(self.alternateAllele[position])
    self.sortedVariantPositions.pop(0)

# Add sequence to the output reference.  Check if there was a previous variant
# in this reference sequence.  If so, the sequence up to and including that variant
# is already included in the output reference.
  def writeToOutputSequence(self, var, ref, refID, position, variantID, sequence, sequenceStart):
    start = 0
    if variantID > 0:
      oldPosition = self.sortedVariantPositions[variantID - 1]
      if var.variantType[refID][oldPosition][0] == "ins": variantLength = 1
      else: variantLength = var.variantType[refID][oldPosition][1]
      start = self.sortedVariantPositions[variantID - 1] - sequenceStart + variantLength
    ref.outputSequence += sequence[max(0, start):min((position - sequenceStart), len(sequence))]
