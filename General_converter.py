#!/usr/bin/env python
# coding: utf-8

import os
import re
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from splash import Spectrum, SpectrumType, Splash
import time
import json

atom = {"H":1.008, "C": 12.011, "O":15.999, "N": 14.007, "P":30.976, "S": 32.065}
inputdir = "/Users/donghanli/Documents/GitHub/convert_MassBank_format/KI_GIAR/"
inputfile = "KI-GIAR_zic-HILIC_Pos_v0.90.msp"
outputdir ="/Users/donghanli/Documents/GitHub/convert_MassBank_format/KI_GIAR/KI_Gumma/"
with open(inputdir+inputfile, "r") as f:
    readdata = f.read()
filelist = readdata.split("\n\n")
filelist = filelist[:-2]
with open(inputdir+"classlist.json","r") as f:
    classlist =json.loads(f.read())
data = []
faillist = []
successlist=[]
finallist =[]
name, formula, smiles, mol, inchi, instrmt_type, msn="","","","","","",""

def dateformat(date):
    date = date.replace("/",".")
    date = date.replace("-",".")
    y,m,d = date.split(".")
    MBdate = y+"."+format(int(m),"02d")+"."+format(int(d),"02d")
    return(MBdate)

def calculate_weight(chem):     # <string> formula.
	s = re.findall('([A-Z][a-z]?)([0-9]*)', chem)
	compoundweight = 0
	for element, count in s:
		count = int(count or '1')
		compoundweight += atom[element] * count
	return compoundweight      # <float> formula weight

def searchsingle(tokeyword,fromkeyword1, fromkeyword2="\n", search = "ON"):
    temp = fromkeyword1
    tempfield = ""
    if search == "ON":
        temp = data.split(fromkeyword1)[1].split(fromkeyword2)[0]
        if temp.upper() == "POSITIVE" or temp.upper() =="NEGATIVE":
            temp = temp.upper()
        elif temp == "MS1":
            temp = "MS"
        tempfield = tokeyword+temp
    elif search =="OFF":
        tempfield = tokeyword+fromkeyword1
    finallist.append(tempfield)
    return temp

def searchmulti(tokeyword, fromkeyword, ignore = ""):
    searchlen = len(fromkeyword)
    ignorelen = len(ignore)
    temp = [k for k in data.split("\n") if k[:searchlen] == fromkeyword]
    temp = list(set(temp))
    if ignore != "":
        for j in reversed(range(len(temp))):
            if temp[j][-ignorelen:]==ignore:
                del temp[j]
    tempfield = "\n".join(temp)
    tempfield = tempfield.replace(fromkeyword, tokeyword)
    finallist.append(tempfield)
    return temp

def peakprocess(beginline, endline="END", annotation=""):
    peakX=[]
    peakY=[]
    anno=[]
    if endline =="END":
        peakinfo = data.split(beginline)[1].split("\n")[1:]
    else:
        peakinfo1 = data.split(beginline)[1].split("\n")[1:]
        peakinfo2 = "\n".join(peakinfo1).split(endline)[0]
        peakinfo = peakinfo2.split("\n")[:-1]
    if annotation != "":
        tempinfo = [peakinfo[0]]
        for i in range(1,len(peakinfo)):
            if (peakinfo[i-1].split(annotation)[0] != peakinfo[i].split(annotation)[0]):
                tempinfo.append(peakinfo[i])
        peakinfo = tempinfo
        for i in range(len(peakinfo)):
            calpeakinfo = peakinfo[i].split(annotation)[0]
            calpeakinfo = calpeakinfo.strip()
            calpeakinfo = calpeakinfo.replace(" ","\t")
            peakX.append(float(calpeakinfo.split("\t")[0]))
            peakY.append(float(calpeakinfo.split("\t")[1]))
        pairpeak = list(sorted(zip(peakX, peakY)))
        spectrum = Spectrum(pairpeak, SpectrumType.MS)
        SPNO = "PK$SPLASH: "+ Splash().splash(spectrum)
        pknum= "PK$NUM_PEAK: " + str(len(pairpeak))
        maxY=max(peakY)
        pklist="PK$PEAK: m/z int. rel.int.\n"
        for i in range(len(pairpeak)):
            pklist+=("  "+str(pairpeak[i][0])+" "+str(pairpeak[i][1])+" "+str(round(pairpeak[i][1]/maxY*999))+"\n")
        pklist+="//"
        pkanno ="PK$ANNOTATION: m/z annotation\n"
        for i in range(len(peakinfo)):
            anno= peakinfo[i].split("\"")[1]
            anno = anno.replace(";", "\\")
            anno = anno.replace(" ","_")
            anno = anno.replace("\\_","\\")
            pkanno +=("  "+str(pairpeak[i][0])+" "+anno+"\n")
        pkanno = pkanno[:-2]
        finallist.extend([SPNO,pkanno,pknum,pklist])
        return None
    elif annotation == "":
        for i in range(len(peakinfo)):
            calpeakinfo = peakinfo[i].strip()
            calpeakinfo = calpeakinfo.replace(" ","\t")
            peakX.append(float(calpeakinfo.split("\t")[0]))
            peakY.append(float(calpeakinfo.split("\t")[1]))
        pairpeak = list(sorted(zip(peakX, peakY)))
        spectrum = Spectrum(pairpeak, SpectrumType.MS)
        SPNO = "PK$SPLASH: "+ Splash().splash(spectrum)
        pknum= "PK$NUM_PEAK: " + str(len(pairpeak))
        maxY=max(peakY)
        pklist="PK$PEAK: m/z int. rel.int.\n"
        for i in range(len(pairpeak)):
            pklist+=("  "+str(pairpeak[i][0])+" "+str(pairpeak[i][1])+" "+str(round(pairpeak[i][1]/maxY*999))+"\n")
        pklist+="//"
        finallist.extend([SPNO,pknum,pklist])
        return None

def inserttitle(n,i,m):
    title="RECORD_TITLE: "+n+"; "+i+"; "+m
    finallist.insert(1,title)
    return None

def standardASCii(finallist):
    newstring = "\n".join(finallist)
    newstring = newstring.replace("á","a")
    newstring = newstring.replace("‐","-")
    newstring = newstring.replace("™","")
    newstring = newstring.replace(" "," ")
    newstring = newstring.replace("′","'")
    newstring = newstring.replace(" 　","")
    newstring = newstring.replace("　","")
    newstring = newstring.replace("フジッコ,","")
    newstring = newstring.replace("α","alpha")
    newstring = newstring.replace("µm","um")
    return newstring

def adjustor2(base = "SMILES"):
    if base == "SMILES":
        s = smiles
        m = Chem.MolFromSmiles(s)
        i = Chem.MolToInchi(m)
        ik= Chem.inchi.MolToInchiKey(m)
        f = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if f[-1]=="+" or f[-1]=="-":
            f = "["+f[:-1]+"]"+f[-1]
        w = str(round(calculate_weight(f),4))
    elif base == "INCHI":
        i = inchi
        m = Chem.MolFromInchi(i)
        s = Chem.MolToSmiles(m)
        ik= Chem.inchi.MolToInchiKey(m)
        f = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if f[-1]=="+" or f[-1]=="-":
            f = "["+f[:-1]+"]"+f[-1]
        w = str(round(calculate_weight(f),4))
    templist =["CH$SMILES: "+s if item.find("CH$SMILES: ")==0 else item for item in finallist]
    templist =["CH$IUPAC: "+i if item.find("CH$IUPAC: ")==0 else item for item in templist]
    templist =["CH$LINK: INCHIKEY "+ik if item.find("CH$LINK: INCHIKEY ")==0 else item for item in templist]
    templist =["CH$FORMULA: "+f if item.find("CH$FORMULA:")==0 else item for item in templist]
    templist =["CH$EXACT_MASS: "+w if item.find("CH$EXACT_MASS: ")==0 else item for item in templist]
    return templist

def adjustor():
    tempf = formula
    mol = Chem.MolFromSmiles(smiles)
    f2 = Chem.rdMolDescriptors.CalcMolFormula(mol)
    if f2[:-1] == formula:
        tempf= "["+f2[:-1]+"]"+f2[-1]
    elif f2 != formula and formula != f2[:-1]:
        faillist.append(name + " Bad SMILES/formula: ")
        print("Fail: Bad SMILES/Formula.")
        return False
    templist =["CH$FORMULA: "+tempf if item.find("CH$FORMULA:")==0 else item for item in finallist]
    templist =["CH$EXACT_MASS: "+str(round(calculate_weight(tempf),4)) if item.find("CH$EXACT_MASS: ")==0 else item for item in templist]
    return templist


k1=time.time()
smilelist=[]
for i, data in enumerate(filelist):
    finallist =[]
    name, formula, smiles, mol, inchi, instrmt_type, msn="","","","","","",""
    #data = filelist[i]
    print("i=",i)
    prefix = "TMP"
    accession = searchsingle("ACCESSION: "+prefix, format(i+1, "0"+str(8-len(prefix))+"d"), search="OFF")
    searchsingle("DATE: ","2019.06.10", search = "OFF")
    searchsingle("AUTHORS: ","AUTHORS: ")
    searchsingle("LICENSE: ","LICENSE: ")
    searchsingle("COMMENT: ","COMMENT: ")
    name = searchsingle("CH$NAME: ","NAME: ",";")
    inchikey = searchsingle("inchikey:", "INCHIKEY: ")
    inchikey = "InChIKey="+inchikey
    finallist.pop()
    try:
        class1 = classlist[inchikey]
    except KeyError:
        print("No such class")
        continue
    else:
        searchsingle("CH$COMPOUND_CLASS: ",class1, search = "OFF")
    #class1 = classlist[inchikey]
    #searchsingle("CH$COMPOUND_CLASS: ",class1, search = "OFF")
    formula = searchsingle("CH$FORMULA: ","FORMULA: ")
    searchsingle("CH$EXACT_MASS: ",str(round(calculate_weight(formula),4)),search = "OFF")
    smiles = searchsingle("CH$SMILES: ","SMILES: ")
    print(smiles)
    mol = Chem.MolFromSmiles(smiles)
    inchi = searchsingle("CH$IUPAC: ", "INCHI: ")
    searchsingle("CH$LINK: INCHIKEY ","INCHIKEY: ")
    searchsingle("AC$INSTRUMENT: ", "INSTRUMENT: ")
    instrmt_type = searchsingle("AC$INSTRUMENT_TYPE: ", "INSTRUMENTTYPE: ")
    msn = searchsingle("AC$MASS_SPECTROMETRY: MS_TYPE ", "MSLEVEL: ")
    searchsingle("AC$MASS_SPECTROMETRY: ION_MODE ","IONMODE: ")
    searchsingle("AC$MASS_SPECTROMETRY: IONIZATION ", "IONISATION: ")
    searchsingle("AC$MASS_SPECTROMETRY: COLLISION_ENERGY ", "COLLISIONENERGY: ")
    searchsingle("AC$CHROMATOGRAPHY: COLUMN_NAME ","COLUMN: ")
    searchsingle("AC$CHROMATOGRAPHY: FLOW_GRADIENT ", "FlowGradient: ")
    searchsingle("AC$CHROMATOGRAPHY: FLOW_RATE ","FlowRate: ")
    searchsingle("AC$CHROMATOGRAPHY: RETENTION_TIME ", "RETENTIONTIME: ")
    searchsingle("AC$CHROMATOGRAPHY: SOLVENT A ", "SolventA: ")
    searchsingle("AC$CHROMATOGRAPHY: SOLVENT B ", "SolventB: ")
    searchsingle("MS$FOCUSED_ION: PRECURSOR_TYPE ", "PRECURSORTYPE: ")
    searchsingle("MS$FOCUSED_ION: PRECURSOR_M/Z ", "PRECURSORMZ: ")
    peakprocess("Num Peaks", annotation="\"")
    inserttitle(name,instrmt_type,msn)
    finallist = adjustor()
    #finallist = adjustor2("INCHI")
    if finallist != False:
        finalstr = standardASCii(finallist)
        print(prefix+accession+" success.")
        smilelist.append(smiles)
        with open(outputdir+prefix+accession+".txt", "w") as f:
            f.write(finalstr)
k2=time.time()
print(k2-k1)
with open(outputdir+"faillist.txt", "w") as f:
    f.writelines(faillist)
