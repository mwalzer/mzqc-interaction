#!/usr/bin/env python
__author__ = 'walzer'
import os
import sys
import logging
from os import listdir as ld
from os.path import isfile, isdir, join, basename, splitext
from lxml import etree
from datetime import datetime, timedelta
from collections import defaultdict
from typing import List, Dict, Union, Any
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
import hashlib
from mzqc import MZQCFile as qc
import pronto
import click
import logging

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
INFO = '''
A simple qcML to mzQC converter in python. 
'''

def print_help():
	"""
	Print the help of the tool
	:return:
	"""
	ctx = click.get_current_context()
	click.echo(ctx.get_help())
	ctx.exit()

def safe_get_xpath(tree: etree.ElementTree, path: str, attrib: str):
	p = tree.xpath(path)
	if p:
		e = next(iter(p), None)
		if e: 
			return e.attrib.get(attrib) 

def sha256fromfile(abs_file_path: str) -> str:
	"""
	sha256fromfile will create a sha256 digest from the file at given path.

	To preserve memory and speed up the digest,
	the file is digested with the help of a memoryview and hashlib.sha256().update.

	Parameters
	----------
	abs_file_path : str
			The absolute path to the file to digest

	Returns
	-------
	str
			The cast or unchanged argument

	Raises
	------
	FileNotFoundError
			If abs_file_path is not a file  
	"""
	sha = hashlib.sha256()
	b = bytearray(128 * 1024)
	mv = memoryview(b)

	with open(abs_file_path, 'rb', buffering=0) as f:
		for n in iter(lambda: f.readinto(mv), 0):
			sha.update(mv[:n])
	return sha.hexdigest()

def pull_mzml_data(mzml_path):
	tree = etree.parse(mzml_path)

	tics_s = [(ms1.getparent().getparent().getparent().attrib.get("id",None),
	    		ms1.attrib.get("value",None),
				ms1.getparent().getparent().getparent().attrib.get("defaultArrayLength",None),
				next(iter(ms1.getparent().getparent().getparent().xpath(".//*[@accession='MS:1000285']"))).attrib.get("value",None))\
				for ms1 in tree.xpath(".//*[@accession='MS:1000016']") if len(ms1.getparent().getparent().getparent().xpath(".//*[@accession='MS:1000579']"))==1]
	# .xpath(".//*[local-name() means relative to calling node
	tics_s = pd.DataFrame(tics_s, columns=["MS:1000767", "MS:1000894", "MS:1003059", "MS:1000285"])

	nid_rt_lookup = [(ms1.getparent().getparent().getparent().attrib.get("id",None),
	    		ms1.attrib.get("value",None)) for ms1 in tree.xpath(".//*[@accession='MS:1000016']")]
	nid_rt_lookup = pd.DataFrame(nid_rt_lookup, columns=["MS:1000767", "MS:1000894"])
	nid_rt_lookup["MS:1000894"] = nid_rt_lookup["MS:1000894"].astype(float) * 60

	# Instrument Type
	rpg = [cvp for cvp in tree.xpath("//*[local-name() = 'referenceableParamGroup']")[0] if cvp.tag.endswith("cvParam")]
	psi_ms_url = "https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v4.1.130/psi-ms.obo"
	ms = pronto.Ontology(psi_ms_url, import_depth=0)
	cv_instruments  = {x.id for x in ms['MS:1000031'].subclasses().to_set()}
	mzml_instrument = {tag.attrib.get("accession",None) for tag in rpg if tag.attrib.get("accession",None) in cv_instruments}
	if len(mzml_instrument) > 1:
		logging.warn("Provided mzML has more than one instrument registered, ignoring all but first.")
	itype = ms.get(next(iter(mzml_instrument)))

	return itype, nid_rt_lookup

def convert_qcml2mzqc(qcml_input, pull=False, rawid=False):
	tree = etree.parse(qcml_input)

	completiontime = datetime.strptime(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='MS:1000747']")[0].attrib.get('value'), "%Y-%m-%d")

	# pull metric values
	ms1 = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000006']")[0].attrib.get('value'))
	ms2 = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000007']")[0].attrib.get('value'))
	nc = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000008']")[0].attrib.get('value'))
	mzrange = list(map(float, tree.xpath("//*[local-name() = 'attachment'][@accession='QC:0000009']/*[local-name() = 'table']/*[local-name() = 'tableRowValues']")[0].text.split()))
	rtrange = list(map(float, tree.xpath("//*[local-name() = 'attachment'][@accession='QC:0000012']/*[local-name() = 'table']/*[local-name() = 'tableRowValues']")[0].text.split()))
	fc = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000046']")[0].attrib.get('value'))
	idfc = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000058']")[0].attrib.get('value'))
	jump = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000059']")[0].attrib.get('value'))
	slump = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000060']")[0].attrib.get('value'))

	psms = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000029']")[0].attrib.get('value'))  # as recorded in table
	id_peptides = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000030']")[0].attrib.get('value'))  # peptidoforms
	id_proteins = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000032']")[0].attrib.get('value'))  # proteoforms
	id_uniq_proteins = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000033']")[0].attrib.get('value'))  # accessions
	# id_uniq_peptides = int(tree.xpath("//*[local-name() = 'qualityParameter'][@accession='QC:0000031']")[0].attrib.get('value'))  # no fitting term?

	proto_tics = [''.join(row.itertext()).split() for row in tree.xpath("//*[local-name() = 'attachment'][@accession='QC:0000056']")[0][0]]
	proto_tics.reverse()
	col_names = proto_tics.pop()
	proto_tics.reverse()
	tics = pd.DataFrame(proto_tics, columns = col_names)
	tics.rename({"MS:1000894_[sec]": "MS:1000894", "peak_count": "MS:1003059"}, axis='columns', inplace=True)
	tics.drop(columns=["S/N"], inplace=True)
	# tics = pd.to_numeric(tics[tics.columns].stack(), errors='coerce').unstack()  # to indiscriminant
	tics["MS:1000894"] = tics["MS:1000894"].astype(float)
	tics["MS:1000285"] = tics["MS:1000285"].astype(int)
	tics["MS:1003059"] = tics["MS:1003059"].astype(int)
	# missing native id come in trough table join later if mzML pull is provided

	proto_iip = [''.join(row.itertext()).split() for row in tree.xpath("//*[local-name() = 'attachment'][@accession='QC:0000018']")[0][0]]
	proto_iip.reverse()
	col_names = proto_iip.pop()
	proto_iip.reverse()
	iip = pd.DataFrame(proto_iip, columns = col_names)
	iip.rename({"MS:1000894_[sec]": "MS:1000894"}, axis='columns', inplace=True)
	iip["MS:1000894"] = iip["MS:1000894"].astype(float)
	iip["MS:1000927"] = iip["MS:1000927"].astype(float)
	# missing native id come in trough table join later if mzML pull is provided

	proto_ppm = [''.join(row.itertext()).split() for row in tree.xpath("//*[local-name() = 'attachment'][@accession='QC:0000038']")[0][0]]
	proto_ppm.reverse()
	col_names = proto_ppm.pop()
	proto_ppm.reverse()
	# proto_ppm = [a + [None] * (len(col_names) - len(a)) for a in proto_ppm]
	ppm = pd.DataFrame(proto_ppm, columns = col_names)
	ppm.rename({"RT": "MS:1000894", "delta_ppm": "MS:4000072", 
				"PeptideSequence": "MS:1003169",
				"Charge": "MS:1000041"}, axis='columns', inplace=True)
	ppm.drop(columns=["Score","MZ", "TheoreticalWeight", "Oxidation_(M)", "Acetyl_(N-term)"], inplace=True)   
	ppm["MS:1000894"] = ppm["MS:1000894"].astype(float)
	ppm["MS:4000072"] = ppm["MS:4000072"].astype(float)
	ppm["MS:1000041"] = ppm["MS:1000041"].astype(int)

	if not rawid:
		expected_peptides = ['YAEAVTR','STLTDSLVC(Carbamidomethyl)K','SLADELALVDVLEDK',
							'NPDDITNEEYGEFYK','LAVDEEENADNNTK','FEELNMDLFR',
							'EAALSTALSEK','DDVAQTDLLQIDPNFGSK','RFPGYDSESK',
							'EATTEFSVDAR','EQFLDGDGWTSR','TPAQFDADELR','LGDLYEEEMR',
							'EVSTYIK','FAFQAEVNR']  # QC2 sample peptides as defined in MS:4000078
		ppm = ppm[ppm["MS:1003169"].isin(expected_peptides)]
			
	mzml_path = tree.xpath("//*[local-name() = 'qualityParameter'][@accession='MS:1000577']")[0].attrib.get('value')
	if not mzml_path.lower().endswith('.mzML', ):
		mzml_path = '.'.join([mzml_path,'mzML'])
	mzid_path = '.'.join([os.path.splitext(mzml_path)[0], "mzid.gz"])
	mzqc_path = '.'.join([os.path.splitext(mzml_path)[0], "mzqc"])

	psi_ms_url = "https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v4.1.130/psi-ms.obo"
	ms = pronto.Ontology(psi_ms_url, import_depth=0)

	if pull:
		chksm = sha256fromfile(mzml_path)
		itype, nid_rt_lookup = pull_mzml_data(mzml_path)

		tics['merge'] = np.round(tics["MS:1000894"], decimals=6)
		iip['merge'] = np.round(iip["MS:1000894"], decimals=6)
		nid_rt_lookup['merge'] = np.round(nid_rt_lookup["MS:1000894"], decimals=6)
		
		tics = pd.merge(tics, nid_rt_lookup, how="inner", on="merge", sort=True,)
		tics.drop(columns=['merge','MS:1000894_y'],inplace=True)
		tics.rename({"MS:1000894_x": "MS:1000894"}, axis='columns', inplace=True)
		iip = pd.merge(iip, nid_rt_lookup, how="inner", on="merge", sort=True,)
		iip.drop(columns=['merge','MS:1000894_y'],inplace=True)
		iip.rename({"MS:1000894_x": "MS:1000894"}, axis='columns', inplace=True)
	else:
		chksm = None
		instruments  = ms['MS:1000031'].subclasses().to_set()
		instrumentname = tree.xpath("//*[local-name() = 'qualityParameter'][@accession='MS:1000031']")[0].attrib.get('value')
		itype = {t.name: t for t in instruments}.get(instrumentname, ms['MS:1000031'])
		logging.warning("Conversion without mzML access may result in some incomplete metrics (e.g. for instrument native spectrum identifiers).")

	cv_map = {
		'MS:4000071': nc,
		'MS:4000059': ms1,
		'MS:4000060': ms2,
		'MS:4000102': fc,
		'MS:4000103': idfc,
		'MS:4000097': jump,
		'MS:4000098': slump,
		'MS:4000104': tics.where((pd.notnull(tics)), None).to_dict(orient='list'), 
		'MS:4000105': iip.where((pd.notnull(iip)), None).to_dict(orient='list'), 
		'MS:4000078': ppm.where((pd.notnull(ppm)), None).to_dict(orient='list'),
		'MS:4000069': mzrange,
		'MS:4000070': rtrange,
		'MS:1003251': psms,
		'MS:1003250': id_peptides,
		'MS:1003328': id_proteins,
		#'MS:???': id_uniq_peptides,
		'MS:1002404': id_uniq_proteins,
	}
	qcl = [qc.QualityMetric(accession=acc,name=ms.get(acc,ms["MS:4000001"]).name,description=str(ms.get(acc,ms["MS:4000001"]).definition),value=val,) for acc, val in cv_map.items()] 

	# assemble mzQC
	infi1 = qc.InputFile(name=mzml_path, location=mzml_path, fileFormat=qc.CvParameter(accession="MS:1000584", name="mzML format"))
	if chksm:
		infi1.fileProperties.append(qc.CvParameter(accession="MS:1003151", name="SHA-256", value=chksm))
	infi1.fileProperties.append(qc.CvParameter(accession=itype.id, name=itype.name))
	infi1.fileProperties.append(qc.CvParameter(accession="MS:1000747", name="completion time", value=completiontime))
	infi2 = qc.InputFile(name=mzid_path, location=mzid_path, 
		      fileFormat=qc.CvParameter(accession="MS:1002073", name="mzIdentML format"))
	anso1 = qc.AnalysisSoftware(accession="MS:1002251", name="qcloud", version="not latest version", uri="https://github.com/proteomicsunitcrg/qcloud2-pipeline")
	anso2 = qc.AnalysisSoftware(accession="MS:1003357", name="qcml2mzqc_converter", version="0", uri="https://github.com/MS-Quality-Hub/mzqclib-manuscript")
	meta = qc.MetaDataParameters(inputFiles=[infi1, infi2],analysisSoftware=[anso1, anso2], label="qcloud conversion")
	rq = qc.RunQuality(metadata=meta, qualityMetrics=qcl)
	cv = qc.ControlledVocabulary(name="PSI-MS", uri="https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v4.1.130/psi-ms.obo", version="v4.1.130")
	mzqc = qc.MzQcFile(version="1.0.0", description="mzQC converted from qcloud qcML", contactName="mwalzer", 
			contactAddress="https://github.com/MS-Quality-Hub/mzqclib-manuscript", runQualities=[rq], controlledVocabularies=[cv]) 
	return mzqc_path, mzqc

@click.command
@click.help_option('-h', '--help')
@click.argument('qcml_input', type=click.Path(exists=True,readable=True) )  # help="The file with the spectra to analyse"
@click.option('--log', type=click.Choice(['debug', 'info', 'warn'], case_sensitive=False),
	default='warn', show_default=True,
	required=False, help="Log detail level. (verbosity: debug>info>warn)")
@click.option("--pull", is_flag=True, show_default=True, default=False, help="Try to pull extra metadata from mzML (must be present in same folder and named in the qcml file).")
def qcml2mzqc_converter(qcml_input, log, pull):
	"""
	Will convert an OpenMS::QCCalculator qcML for qcloud into a mzQC v1.0.0 file.
	"""
	# set loglevel - switch to match-case for py3.10+
	lev = {'debug': logging.DEBUG,
		'info': logging.INFO,
		'warn': logging.WARN }
	logging.basicConfig(format='%(levelname)s:%(message)s', level=lev[log])

	if not qcml_input:
		print_help()

	mzqc_path, qcobj = convert_qcml2mzqc(qcml_input=qcml_input, pull=pull)

	try:
		with open(mzqc_path, "w") as file:
			file.write(qc.JsonSerialisable.ToJson(qcobj, readability=1))
	except Exception as e:
		click.echo(e)
		click.echo("Failed to write converted file next to input. Adjust permissions.")
		print_help()

if __name__ == '__main__':
	qcml2mzqc_converter()