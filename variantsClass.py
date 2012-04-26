#!/usr/bin/python

import os.path
import sys
import random
import re

class variants:

# Constructor.
  def __init__(self, options):

    # Define a list that contains user specified variants.
    self.numberDefinedVariants   = 0
    self.definedVariantTypes     = []
    self.definedVariantLengths   = []
    self.definedVariantRefIDs    = []
    self.definedVariantPositions = []

    # Define the hash table for storing the created variants
    self.variantType             = {}

    # Set the number of the different variant types to the values
    # specified in the command line.
    self.numberSnps              = options.numberSnps
    self.numberMnps              = options.numberMnps
    self.numberIndels            = options.numberIndels

    # Set some other parameters associated with the variants and
    # specified on the command line.
    self.maximumMnpLength        = options.maximumMnpLength
    self.maximumIndelLength      = options.maximumIndelLength
    self.flankLength             = options.flankLength

    # Define the string that stores sequence to be written to the
    # output file and stores sequence for constructing flanking
    # sequence in the alternate contigs.
    self.storedSequence          = ""
    self.storedFlank             = ""

    # Define variables that keep track of the position within the
    # reference sequences.
    self.refStart                = -1
    self.refEnd                  = -1
    self.seqStart                = -1
    self.seqEnd                  = -1

    # Generate an array that assigns a number to each nucleotide.
    self.nucleotideArray         = []
    self.nucleotideArray.append("A")
    self.nucleotideArray.append("C")
    self.nucleotideArray.append("G")
    self.nucleotideArray.append("T")

# Check the defined variants.
  def checkDefinedVariants(self, ref):

    # Loop over the variants.
    for i, variantType in enumerate(self.definedVariantTypes):
      refID    = self.definedVariantRefIDs[i]
      position = self.definedVariantPositions[i]
      length   = self.definedVariantLengths[i]

      # Check that the specified reference ID exists.
      if refID >= ref.numberOfSequences:
        print >> sys.stderr, "ERROR: Specified reference ID does not exist for variant: ", variantType, "at ", refID + 1, ":", position + 1
        print >> sys.stderr, "  Max reference ID (1-based): ", ref.numberOfSequences
        exit(1)

      # Check that the position specified exists within the reference
      # sequence.  This depends on the variant type.  Also, check if
      # another variant overlaps with the defined variant.  If so, throw
      # an exception.  If not, add the variant to the list of mutations
      # to apply to the reference.
      error   = False
      overlap = False

      # SNPs just need to have the specified position exist.
      if variantType == "snp":
        if (position < 0) | (position > ref.sequenceEnd[refID]): error = True
        if self.variantType.has_key(refID):
          if self.variantType[refID].has_key(position): overlap = True

        if overlap == False:
          try: self.variantType[refID][position] = ("snp", 1)
          except:
            self.variantType[refID] = {}
            self.variantType[refID][position] = ("snp", 1)

      # MNPs need to have the specified position plus length - 1 exist.  Only
      # in this way are there enough bases left in the reference sequence to
      # mutate length bases.
      elif variantType == "mnp":
        if (position < 0) | ( (position + length - 1) > ref.sequenceEnd[refID]): error = True
        for i in range(position, position + length):
          if self.variantType.has_key(refID):
            if self.variantType[refID].has_key(position): overlap = True

        if overlap == False:
          try: self.variantType[refID][position] = ("mnp", length)
          except:
            self.variantType[refID] = {}
            self.variantType[refID][position] = ("mnp", length)
          for i in range(position + 1, position + length):
            self.variantType[refID][i] = ("mnp", "remove")

      # Insertions just require the specified position to exist.  For insertions and
      # deletions, increase the length by one as the alleles contain an anchor base
      # as well as the inserted/deleted bases.
      elif variantType == "ins":
        length += 1
        if (position < 0) | (position > ref.sequenceEnd[refID]): error = True
        if self.variantType.has_key(refID):
          if self.variantType[refID].has_key(position): overlap = True

        if overlap == False:
          try: self.variantType[refID][position] = ("ins", length + 1)
          except:
            self.variantType[refID] = {}
            self.variantType[refID][position] = ("ins", length + 1)

      # Deletions require the specified base (the anchor base) and the next
      # length bases to exist so that they can be deleted.
      elif variantType == "del":
        length += 1
        if (position < 0) | ( (position + length) > ref.sequenceEnd[refID]): error = True
        for i in range(position, position + length + 1):
          if self.variantType.has_key(refID):
            if self.variantType[refID].has_key(i): overlap = True
 
        if overlap == False:
          try: self.variantType[refID][position] = ("del", length)
          except:
            self.variantType[refID] = {}
            self.variantType[refID][position] = ("del", length)
          for i in range(position + 1, position + length):
            self.variantType[refID][i] = ("del", "remove")
      
      # If the position is invalid, terminate the program,
      if error == True:
        print >> sys.stderr, "ERROR: Specified position does not exist in reference ID: ", variantType, "at ", refID + 1, ":", position + 1
        exit(1)

      # If the variant overlaps with an already defined variant, throw an
      # exception.
      if overlap == True:
        print >> sys.stderr, "ERROR: Overlap found with another specified variant: ", variantType, "at ", refID + 1, ":", position + 1
        exit(1)

# Create SNPs.
  def createSnps(self, ref):
    print >> sys.stdout, "Creating SNPs."
    for i in range(self.numberSnps):

      # Pick a reference sequence at random and then a position at random within the
      # sequence.
      snpRefID    = random.randrange(0, ref.numberOfSequences - 1)
      snpPosition = random.randrange(0, ref.sequenceEnd[snpRefID], 1)

      # If a SNP already exists at this position, pick a new position.
      if self.variantType.has_key(snpRefID):
        while self.variantType[snpRefID].has_key(snpPosition):
          snpPosition = random.randrange(0, ref.sequenceEnd[snpRefID] - 1, 1)
      else: self.variantType[snpRefID] = {}

      self.variantType[snpRefID][snpPosition] = ("snp", 1)

# Create MNPs.
  def createMnps(self, ref):
    print >> sys.stdout, "Creating MNPs."
    for i in range(self.numberMnps):
      mnpRefID    = random.randint(0, ref.numberOfSequences - 1)
      mnpPosition = random.randrange(0, ref.sequenceEnd[mnpRefID] - 1, 1)
      mnpLength   = random.randint(2, self.maximumMnpLength)

      # If the MNP overlaps an already defined variant, pick a new position.
      while True:
        overlap = False
        if self.variantType.has_key(mnpRefID):
          for i in range(mnpPosition, mnpPosition + mnpLength):
            if self.variantType[mnpRefID].has_key(i): overlap = True
        else: self.variantType[mnpRefID] = {}

        # If an overlap occurs, pick a new position.
        if overlap == True:
          mnpPosition = random.randrange(0, ref.sequenceEnd[mnpRefID] - 1, 1)
        else: break

      # Treat the MNP as a set of consecutive SNPs so that the routines for
      # generating SNPs can be used to generate MNPs.
      self.variantType[mnpRefID][mnpPosition] = ("mnp", mnpLength)
      for i in range(mnpPosition + 1, mnpPosition + mnpLength):
        self.variantType[mnpRefID][i] = ("mnp", "remove")

# Create indels.
  def createIndels(self, ref):
    print >> sys.stdout, "Creating indels."
    for i in range(self.numberIndels):
      indelRefID    = random.randint(0, ref.numberOfSequences - 1)
      indelPosition = random.randrange(0, ref.sequenceEnd[indelRefID] - 1, 1)
      indelLength   = random.randint(1, self.maximumIndelLength)
      indelType     = random.randint(0, 1)

      # Check for mutations already defined at this position.
      if self.variantType.has_key(indelRefID):
        while self.variantType.has_key(indelPosition):
          indelPosition = random.randrange(0, ref.sequenceEnd[indelRefID] - 1, 1)
      else: self.variantType[indelRefID] = {}

      # Insertions.
      if indelType == 0: self.variantType[indelRefID][indelPosition] = ("ins", indelLength)

      # Deletions
      elif indelType == 1:

        # For deletions, check that the deleted bases do not span another
        # polymorphism..
        while True:
          overlap = False
          for i in range(indelPosition, indelPosition + indelLength + 1):
            if self.variantType[indelRefID].has_key(i): overlap = True

          # If an overlap occurs, pick a new indel location.
          if overlap == True:
            indelPosition = random.randrange(0, ref.sequenceEnd[indelRefID] - 1, 1)
          else: break

        self.variantType[indelRefID][indelPosition] = ("del", indelLength)
        for i in range(indelPosition + 1, indelPosition + indelLength + 1):
          self.variantType[indelRefID][i] = ("del", "remove")

# Remove all entries except for the anchor base position of each variant.
  def removeExtraEntries(self):
    for variantRefID in self.variantType.keys():
      for variant in self.variantType[variantRefID].keys():
        if self.variantType[variantRefID][variant][1] == "remove": del self.variantType[variantRefID][variant]

# Find the coordinates of the current reference sequence with respect
# to the coordinates of the reference treated as a single sequence (used
# for randomly assigning variant positions in the whole reference).
  def findCoordinates(self, line):
    self.refStart = self.refEnd + 1
    self.refEnd   = self.refStart + len(line) - 1
    self.seqStart = self.seqEnd + 1
    self.seqEnd   = self.seqStart + len(line) - 1

# Generate the alternate allele from the reference.
  def generateAltAllele(self, refID, position, refAllele):
    altAllele = ""

    # First determine the variant type.
    variantClass  = self.variantType[refID][position][0]
    variantLength = self.variantType[refID][position][1]

    # Generate the alternate allele based on the class.
    if variantClass == "snp":       altAllele = self.generateSnp(refAllele)
    elif variantClass == "mnp":     altAllele = self.generateMnp(refAllele, variantLength)
    elif variantClass == "ins":     altAllele = self.generateInsertion(refAllele, variantLength - 2)
    elif variantClass == "del":     altAllele = self.generateDeletion(refAllele, variantLength)

    # Return the alternate allele.
    return altAllele

# Generate a SNP.
  def generateSnp(self, refAllele):
    altAllele = refAllele
 
    # Ensure that the new nucleotide is difference from the reference.
    while refAllele == altAllele: altAllele = self.nucleotideArray[random.randint(0, 3)]

    return altAllele
    
# Generate a MNP.
  def generateMnp(self, refAllele, length):
    altAllele = ""
  
    for i in range(0, length):
      allele = self.nucleotideArray[random.randint(0, 3)]

      # Ensure that the new nucleotide is difference from the reference.
      while allele == refAllele[i]: allele = self.nucleotideArray[random.randint(0, 3)]
  
      altAllele = altAllele + allele

    return altAllele
    
# Generate a insertion.
  def generateInsertion(self, refAllele, length):
    altAllele = refAllele

    # Generate the sequence to insert.
    for i in range(0, length):
      insertedNucleotide = self.nucleotideArray[random.randint(0, 3)]
      altAllele          = altAllele + insertedNucleotide

    return altAllele
    
# Generate a deletion.
  def generateDeletion(self, refAllele, length):
    altAllele = refAllele[0]

    return altAllele
