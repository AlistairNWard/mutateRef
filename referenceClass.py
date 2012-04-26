#!/usr/bin/python

import os.path
import sys
import re

class reference:

# Constructor.
  def __init__(self, options):
    self.numberOfSequences  = -1
    self.sequenceName       = {}
    self.sequenceEnd        = {}
    self.sequenceLength     = options.sequenceLength
    self.outputSequence     = ""

# Open the input fasta file.
  def openInputFasta(self, options):
    if options.refFile == "stdin":
      self.filehandle = sys.stdin
      self.filename   = "stdin"
    else:
      try: self.filehandle = open(options.refFile, "r")
      except IOError:
        print >> sys.stderr, "Failed to find file: ", options.refFile
        exit(1)
      self.filename = os.path.abspath(options.refFile)

# Open the output reference files.
  def openOutputFiles(self, options):
    output = "output" if options.output == None else options.output
    refFile    = output + ".fa"
    contigFile = output + ".contigs.fa"
  
    self.outputFilehandle = open(refFile, 'w')
    self.contigFilehandle = open(contigFile, 'w')
  
    self.outputFilename = os.path.abspath(refFile)
    self.contigFilename = os.path.abspath(contigFile)

# Close the input fasta file.
  def closeInputFasta(self):
    self.filehandle.close()

# Close the output reference files.
  def closeOutputFiles(self):
    self.outputFilehandle.close()
    self.contigFilehandle.close()

# Collate information on the input fasta file.
  def getSequenceInformation(self):

  # Check if this is sequence or header.
    if self.line.startswith(">"):
      self.numberOfSequences                    = self.numberOfSequences + 1

      # Find the first space in the line and determine the reference name.
      firstSpace = self.line.index(" ")
      self.sequenceName[self.numberOfSequences] = self.line[1:firstSpace]

      self.sequenceEnd[self.numberOfSequences]  = 0
    else:
      self.sequenceEnd[self.numberOfSequences]  = self.sequenceEnd[self.numberOfSequences] + len(self.line) - 1

# Write out any of the output reference sequence that can be written.
  def writeOutputSequence(self, clear):
    while len(self.outputSequence) > self.sequenceLength:
      print >> self.outputFilehandle, self.outputSequence[0:self.sequenceLength]
      self.outputSequence = self.outputSequence[self.sequenceLength:len(self.outputSequence)]

    # If this is the end of the reference sequence and there is still some
    # sequence left to be output, output it now.
    if (clear == True) & (len(self.outputSequence) > 0):
      print >> self.outputFilehandle, self.outputSequence
      self.outputSequence = ""
