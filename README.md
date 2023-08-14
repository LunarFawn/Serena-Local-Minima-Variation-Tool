# Serena-Local-Minima-Variation-Tool

Installation Instructions:
Serna is a pip installable package

To install package:
pip install .

To uninstall package pass:
pip uninstall serena-local-minima-variation-tool

User Instructions:

Serena consists of my take and thoughts on how to design the framework for transporting RNA bioinformatics to and from other applications, as well as provide access to novel RNA analysis algorithms and methodologies.

The Sara2 Secondary Structure is the backbone of the framework. It is an object that contains everything needed to represent a RNA secondary structure in dot bracket notation with the primary structure sequence and both free engery and stack energy.

The sara2 secondary structure list is the primary vehicle for performing analysis of the structural ensemble. It is populated with each secondary structure found during for the ensemble. This feeds the SingleEnsembleGroup and then the MultipleEnsemble groups as methods for transporting all the information about the secondary structures of an ensemble around in a object oriented framework.

See unit tests for examples of implementation.

This version 1.0.1 and it is still in beta release
