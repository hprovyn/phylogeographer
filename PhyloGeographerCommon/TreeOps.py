# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 08:04:08 2018

@author: hunte
"""

import xml
import xml.etree.ElementTree as ET

def parseTreeFileReturnRoot(fil):
    if isXml(fil):
        return parseTreeXML(fil)
    else:
        if isJson(fil):
            return parseTreeJSON(fil)
        else:
            hierarchy = parseTreeParentChildren(fil)
            return (getRoot(hierarchy), hierarchy, {})
    
def isXml(fil):
    fr = open(fil,'r')
    if fr.readline()[0] == "<":
        return True
    else:
        return False
    
def isJson(fil):
    fr = open(fil,'r')
    if fr.readline()[0] == "{":
        return True
    else:
        return False
    
def parseTreeParentChildren(fil):
    hierarchy = {}
    fr = open(fil,'r')
    for line in fr:
        splits = line.split(",")
        parent = splits[0]
        for i in range(len(splits) - 1):
            child = splits[i + 1].replace("\n","")
            hierarchy[child] = parent
    return hierarchy

def getRoot(hierarchy):
    keys = [key for key in hierarchy]
    node = keys[0]
    while node in hierarchy:
        node = hierarchy[node]
    return node

def parseTreeXML(fil):
    hierarchy = {}
    ybp = {}
    root = xml.etree.ElementTree.parse(fil).getroot()
    ybp[root.attrib["name"]] = float(root.attrib["formed"])
    recurseTreeXML(root, hierarchy, ybp)
    return (root.attrib["name"], hierarchy, ybp)

def parseTreeXMLTMRCA(fil):
    hierarchy = {}
    ybp = {}
    root = xml.etree.ElementTree.parse(fil).getroot()
    ybp[root.attrib["name"]] = float(root.attrib["formed"])
    tmrca = {}
    tmrca[root.attrib["name"]] = float(root.attrib["tmrca"])
    recurseTreeXMLTMRCA(root, hierarchy, ybp, tmrca)
    return (root.attrib["name"], hierarchy, ybp, tmrca)

def recurseTreeXML(node, hierarchy, ybp):
    for child in node:    
        if "*" not in child.attrib["name"]:
            hierarchy[child.attrib["name"]] = node.attrib["name"]
            if "formed" in child.attrib and child.attrib["formed"] != "":
                print(child.attrib["formed"])
                ybp[child.attrib["name"]] = float(child.attrib["formed"])
            else:
                ybp[child.attrib["name"]] = node.attrib["formed"]
            recurseTreeXML(child, hierarchy, ybp)

def recurseTreeXMLTMRCA(node, hierarchy, ybp, tmrca):
    for child in node:    
        if "*" not in child.attrib["name"]:
            hierarchy[child.attrib["name"]] = node.attrib["name"]
            if "formed" in child.attrib and child.attrib["formed"] != "":
                print(child.attrib["formed"])
                ybp[child.attrib["name"]] = float(child.attrib["formed"])
            else:
                ybp[child.attrib["name"]] = node.attrib["formed"]
            if "tmrca" in child.attrib and child.attrib["tmrca"] != "":
                print(child.attrib["tmrca"])
                tmrca[child.attrib["name"]] = float(child.attrib["tmrca"])
            else:
                tmrca[child.attrib["name"]] = ybp[child.attrib["name"]]
            recurseTreeXMLTMRCA(child, hierarchy, ybp, tmrca)
import json

def parseTreeJSON(fil):
    hierarchy = {}
    ybp = {}
    root = json.load(open(fil))
    ybp[root["name"]] = float(root["formed"])
    recurseTreeJson(root, hierarchy, ybp)
    return (root["name"], hierarchy, ybp)

def recurseTreeJson(node, hierarchy, ybp):
    if "children" in node:
        for child in node["children"]:
            if "*" not in child["name"]:
                hierarchy[child["name"]] = node["name"]
                theybp = ybp[node["name"]]
                if "formed" in child:
                    theybp = float(child["formed"])
                ybp[child["name"]] = theybp
                recurseTreeJson(child, hierarchy, ybp)

def getChildren(clade, childParents):
    children = []
    for child in childParents:
        if childParents[child] == clade:
            children.append(child)

    return children

def createXML(theroot, childParents, ybp):
    tree = ET.ElementTree()
    root = ET.Element('branch')
    root.set("name", theroot)
    root.set("formed", str(ybp[theroot]))
    addChildrenRecursively(root, childParents, ybp)
    ET.dump(root)
    tree._setroot(root)
    return tree

def createXMLTMRCA(theroot, childParents, ybp, tmrca):
    tree = ET.ElementTree()
    root = ET.Element('branch')
    root.set("name", theroot)
    root.set("formed", str(ybp[theroot]))
    root.set("tmrca", str(tmrca[theroot]))
    print(tmrca)
    addChildrenRecursivelyTMRCA(root, childParents, ybp, tmrca)
    ET.dump(root)
    tree._setroot(root)
    return tree

def createXMLTMRCA_latlon(theroot, childParents, ybp, tmrca, latLons):
    tree = ET.ElementTree()
    root = ET.Element('branch')
    root.set("name", theroot)
    root.set("formed", str(ybp[theroot]))
    root.set("tmrca", str(tmrca[theroot]))
    if theroot in latLons:
        root.set("lat", str(round(latLons[theroot][0],2)))
        root.set("lon", str(round(latLons[theroot][1],2)))
    print(tmrca)
    addChildrenRecursivelyTMRCA_latlon(root, childParents, ybp, tmrca, latLons)
    ET.dump(root)
    tree._setroot(root)
    return tree

def addChildrenRecursively(node,childParents, ybp):
    for child in getChildren(node.attrib['name'], childParents):
        childNode = ET.SubElement(node, 'branch')
        childNode.set("name", child)
        childNode.set("formed",str(ybp[child]))
        if len(getChildren(child, childParents)) > 0:
            addChildrenRecursively(childNode, childParents, ybp)

def addChildrenRecursivelyTMRCA(node,childParents, ybp, tmrca):
    for child in getChildren(node.attrib['name'], childParents):
        childNode = ET.SubElement(node, 'branch')
        childNode.set("name", child)
        childNode.set("formed",str(ybp[child]))
        childNode.set("tmrca", str(tmrca[child]))
        if len(getChildren(child, childParents)) > 0:
            addChildrenRecursivelyTMRCA(childNode, childParents, ybp, tmrca)

def addChildrenRecursivelyTMRCA_latlon(node,childParents, ybp, tmrca, latLons):
    for child in getChildren(node.attrib['name'], childParents):
        childNode = ET.SubElement(node, 'branch')
        childNode.set("name", child)
        childNode.set("formed",str(ybp[child]))
        childNode.set("tmrca", str(tmrca[child]))
        if child in latLons:
            childNode.set("lat", str(round(latLons[child][0],2)))
            childNode.set("lon", str(round(latLons[child][1],2)))
        if len(getChildren(child, childParents)) > 0:
            addChildrenRecursivelyTMRCA_latlon(childNode, childParents, ybp, tmrca, latLons)
            
def convertToXML(projectName):
    #projectName = "R-U106"
    inputFile = "C:\\PhyloGeographer\\preprocessing\\data\\projects\\" + projectName + "\\" + projectName + "-tree.txt"
    outputFile = "C:\\PhyloGeographer\\preprocessing\\data\\projects\\" + projectName + "\\" + projectName + "-tree.xml"
    root,childParents,ybp = parseTreeFileReturnRoot(inputFile)
    writeXML(root, childParents, ybp, outputFile)
    
def writeXML(root, childParents, ybp, outputFile):
    tree = createXML(root, childParents, ybp)
    tree.write(outputFile)
    
def writeXMLTMRCA(root, childParents, ybp, tmrca, outputFile):
    tree = createXMLTMRCA(root, childParents, ybp, tmrca)
    tree.write(outputFile)

def writeXMLTMRCA_latlon(root, childParents, ybp, tmrca, latLons, outputFile):
    tree = createXMLTMRCA_latlon(root, childParents, ybp, tmrca, latLons)
    tree.write(outputFile)
    
def createHTML(theroot, childParents):
    html = "<html><body>"
    for child in childParents.keys():
        html = html + "<a href=" + 'https:\\phylogeographer.com\mygrations?hg=' + theroot + '&clade=' + child + '>' + child + '</a><br>'
    html = html + '</body></html>'
    return html

def createHTMLRecursive(theroot, childParents):
    html = "<html><body>"
    html = html + '<style>body {font-family: "Times New Roman", Times, monospace;} </style>'
    html = html + "<div>"
    html = html + "<a target='_blank' href=" + 'https:\\phylogeographer.com\mygrations?hg=' + theroot + '&clade=' + theroot + '>' + theroot + '</a><br>'
    html = createHTMLRecurse(theroot, childParents, theroot, html, 0)
    html = html + '</div>'
    html = html + '</body></html>'
    return html

def createHTMLRecurse(theroot, childParents, clade, html, level):
    thisLevel = level + 1
    children = getChildren(clade, childParents)
    if len(children) > 0:
        spacer = ""
        for i in range(thisLevel):
            spacer = spacer + "&nbsp;&nbsp;"
        for child in children:
            html = html + "<div>" + spacer
            html = html + "<a target='_blank' href=" + 'https:\\phylogeographer.com\mygrations?hg=' + theroot + '&clade=' + child + '>' + child + '</a><br>'
            html = createHTMLRecurse(theroot, childParents, child, html, thisLevel)
            html = html + "</div>"
    return html
    
def writeHTML(root, childParents, outputFile):
    html = createHTMLRecursive(root, childParents)
    w = open(outputFile,"w+")
    w.write(html)
    w.close()
    
def getTreeFilename(projectName):
    return "C:\\PhyloGeographer\\YFull\\" + projectName + "_YFullxml.txt"

import numpy as np

def upstream(clade, hier):
    if clade in hier:
        parent = hier[clade]
        return np.append([parent], upstream(parent, hier))
    return []

def getAllowableUpstreamFromHaplogroupName(haplogroup, subcladeHaplogroups):
    haploroot, haplohier, haploybp = parseTreeXML(getTreeFilename(haplogroup))
    return getAllowableUpstream(subcladeHaplogroups, haplohier, haploroot)
    
def getAllowableUpstream(subcladeHaplogroups, hierarchy, root):
    allowableUpstream = set([])
    for clade in subcladeHaplogroups:
        theupstream = upstream(clade, hierarchy)
        for up in theupstream:
            allowableUpstream.add(up)
    allowableUpstream.remove(root)
    return allowableUpstream.union(subcladeHaplogroups)
      
def addAllowableUpstreamToHierarchy(hierarchy, ybp, haplogroup, subcladeHaplogroups, prefixToAppend):
    haploroot, haplohier, haploybp = parseTreeXML(getTreeFilename(haplogroup))
    
    for clade in getAllowableUpstream(subcladeHaplogroups, haplohier, haploroot):
        theparent = haplohier[clade]
        if theparent != haploroot:
            theparent = prefixToAppend + theparent
        hierarchy[prefixToAppend + clade] = theparent
        ybp[prefixToAppend + clade] = haploybp[clade]
    