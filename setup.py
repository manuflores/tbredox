import setuptools


setuptools.setup(
	name = "tbredox", 
	version = "0.0.1", 
	author = "Emanuel Flores, Adrian Jinich", 
	author_email = "manuflores {at} caltech {dot} edu",
	description = "Tools for extracting knowledge from RNAseq and TnSeq.", 
	packages = setuptools.find_packages(), 
	classifiers = [
		"Programming Language :: Python :: 3", 
		"License :: OSI Approved :: MIT License", 
		"Operating System :: OS Independent"
	]
)
