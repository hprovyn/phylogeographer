# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 19:03:28 2019

@author: hunte
"""

from PhyloGeographerCommon import TreeOps
import pandas as pd
import numpy as np

from shapely.geometry import Polygon, Point, LinearRing

worldBounds = {'minLat':-56, 'maxLat':75.0, 'minLon':-180.0, 'maxLon':180.0}
ignoreAussieBounds = {'minLat':-56, 'maxLat':0.0, 'minLon':100.0, 'maxLon':180.0}
ignoreAmericasBounds = {'minLat':-56, 'maxLat':75.0, 'minLon':-180.0, 'maxLon':-52.0}
ignoreSAmericaBounds = {'minLat': -56, 'maxLat': 15.0, 'minLon': -140.0, 'maxLon': -7.0}
includeAussie = ["A0-T","C","D","H","M","O","Q","R2","S","F-Y27277","K-Y28299"]
includeAmericas = {"Q": ["Q-Z780", "Q-M3", "Q-F746"],
                   "C": ["C-1373"]}

from PhyloGeographerCommon import DownloadGoogleDrive

def inBounds(lat, lon, bounds):
    return lat > bounds["minLat"] and lat < bounds["maxLat"] and lon > bounds["minLon"] and lon < bounds["maxLon"]

def inBoundsForProject(lat, lon, project, clade, hierarchy):
    if inBounds(lat, lon, worldBounds):
        if inBounds(lat, lon, ignoreAussieBounds):
            if project not in includeAussie:
                return False
            else:
                return True
        else:
            if inBounds(lat, lon, ignoreAmericasBounds) or inBounds(lat, lon, ignoreSAmericaBounds):
                if project not in includeAmericas or all(cladeIsEqualOrDownstream(clade, c, hierarchy) == False for c in includeAmericas[project]):
                    return False
        return True
    else:
        return False
    
def getYFKitsFileName(projectName):
    return "C:\\PhyloGeographer\\YFull\\" + projectName + "-YF-kits.txt"

def isBasal(negatives, clade, hier):
    return set(TreeOps.getChildren(clade, hier))==set(negatives)


def getAncientKits(fil):
    ancientKits = {}
    k = pd.read_csv(fil)
    ids = k["id"]
    clade = k["clade"]
    lats = k["latitude"]
    lons = k["longitude"]
    years = k["ybp"]
    negs = k["negatives"]
    hgs = k["hg"]
    
    excludes = []
    for i in range(len(k)):
        if str(clade[i]) != "nan":
            ancientKit = {}
            
            theid = ids[i]
            if theid not in excludes:
                ancientKit["latitude"] = lats[i]
                ancientKit["longitude"] = lons[i]
                ancientKit["ybp"] = 2019 - years[i]
                ancientKit["clade"] = clade[i]
                ancientKit["negatives"] = negs[i]
                ancientKit["id"] = theid
                hg = hgs[i]
                if hg not in ancientKits:
                    ancientKits[hg] = {}
                ancientKits[hg][theid] = ancientKit
                
                print(theid, "added")
            
    return ancientKits

 
def downloadAndGetAncientKits():
    ancientDownloadFile = "C:\\PhyloGeographer\\2019 Refactor\\AncientDNA.csv"
    ancientFileId = "1uSY4GeZVmzIiJ444d0ndXSeZBvlFkYKFysSK3jDJyds"
    DownloadGoogleDrive.downloadAncient("AncientDNA", ancientFileId, ancientDownloadFile)
    return getAncientKits(ancientDownloadFile)

def cladeIsEqualOrDownstream(clade, eqOrBelow, hierarchy):
    if clade == eqOrBelow:
        return True
    else:
        if clade in hierarchy:
            return cladeIsEqualOrDownstream(hierarchy[clade], eqOrBelow, hierarchy)
        else:
            return False

allAddedToProject = {}
ancientAddedToProject = {}

def parseKitsCSV(fil, hierarchy, project, kitsMap, allAncientKits):
    allAddedToProject[project] = []
    ancientAddedToProject[project] = []
    def addKitIfUseable(kit, clade, project, hierarchy, n):
        if clade in hierarchy:
            if inBoundsForProject(kit["latitude"], kit["longitude"], project, clade, hierarchy):
                if len(TreeOps.getChildren(clade, hierarchy)) > 0:
                    if kit["id"] == "ANCIENT_1374" or (str(n) != 'nan' and str(n) != "" and isBasal(str(n).split(":"), clade, hierarchy)):
                        return kit
                else:
                    return kit
        return None
                    
    def addToKitsMap(kit, clade, isAncient):
        if clade in kitsMap:
            kitsMap[clade].append(kit)
        else:
            kitsMap[clade] = [kit]
        allAddedToProject[project].append(kit["id"])
        if isAncient:
            ancientAddedToProject[project].append(kit["id"])
            
    k = pd.read_csv(fil)
    ids = k["id"]
    lat = k["latitude"]
    lon = k["longitude"]
    clade = k["clade"]
    neg = k["negatives"]
    ybp = k["ybp"]
    
    ancientKits = {}
    addedAncientYFullIds = set([])
    if project in allAncientKits:
        ancientKits = allAncientKits[project]
    if project == "J2":
        ancientKits = {**allAncientKits["J2"], **allAncientKits["J-M241"], **allAncientKits["J2b"], **allAncientKits["J-M205"], **allAncientKits["J-M410"]}
    if project == "A0-T":
        ancientKits = allAncientKits["A0-T"]
        allA0T = ["A1", "A1b", "BT", "CT", "DE", "CF","F","GHIJK","HIJK","IJK","IJ","LT","K2","K2b1","P","P1","NO","R","I"]
        for a in allA0T:
            if a in allAncientKits:
                ancientKits = {**ancientKits, **allAncientKits[a]}
        #ancientKits = {**allAncientKits["A0-T"], **allAncientKits["R"], **allAncientKits["R1"], **all}
    for i in range(len(ids)):
        theid = ids[i]
        thiskit = {"id": theid, "latitude": lat[i], "longitude": lon[i]}
        y = ybp[i]
        if y!= "" and str(y) != 'nan':
            thiskit["ybp"] = y
        else:
            thiskit["ybp"] = 65
        theclade= clade[i]
        
        isAncient = False
        if theid in ancientKits.keys():
            if ancientKits[theid]["clade"] != theclade:
                print("warning:", theid, "in YFull tree", theclade, "in doc", ancientKits[theid]["clade"])
            thiskit["latitude"] = ancientKits[theid]["latitude"]
            thiskit["longitude"] = ancientKits[theid]["longitude"]
            thiskit["ybp"] = ancientKits[theid]["ybp"]
            addedAncientYFullIds.add(theid)
            isAncient = True
        
        addedKit = addKitIfUseable(thiskit, theclade, project, hierarchy, neg[i])
        if addedKit is not None:
            addToKitsMap(addedKit, theclade, isAncient)

    for i in ancientKits:
        anc = ancientKits[i]
        if anc["id"] not in addedAncientYFullIds:
            clade = anc["clade"]
            addedKit = addKitIfUseable(anc, clade, project, hierarchy, anc["negatives"])
            if addedKit is not None:
                addToKitsMap(addedKit, clade, True)
                    
def getYFulltreeFileName(hg):
    return "C:\\PhyloGeographer\\YFull\\" + hg + "_YFullxml.txt"

def getSubProjectTreeFileName(hg):
    return "C:\\PhyloGeographer\\2019 Refactor\\" + hg + "_lessoutliers_treexml.txt"
def getSubProjectSamplesFileName(hg):
    return "C:\\PhyloGeographer\\2019 Refactor\\" + hg + "_samples.txt"

def loadYFSamples(hgs, kitsMap, ancientKits):
    for hg in hgs:        
        (root, hierarchy, ybp, tmrca) = TreeOps.parseTreeXMLTMRCA(getYFulltreeFileName(hg))
        parseKitsCSV(getYFKitsFileName(hg), hierarchy, hg, kitsMap, ancientKits)
        
projects = ["A00","A0","A0-T","A1a","A1b1","B","C","D","E","F-Y27277","G","H","I1","I2","J1","J2","K-Y28299", "L","M","N","O","Q","R1a","R1b","R2","S","T"]
projectsToSeparateIntoHGs = ["A00","A0","A1a","A1b1","B","C","D","E","F-Y27277","G","H","I1","I2","J1","J-M241","J-M410","J-M205","J-Z2453","K-Y28299","L","M","N","O","Q","R1a","R1b","R2","S","T"]
kitsMap = {}

def areLociSuitableForAveraging(loci):
    lons = []
    for loc in loci:
        lons.append(loc[1])
    return isSuitableForAveraging(lons)

def isSuitableForAveraging(lons):
    if len(lons) == 0:
        return True
    return max(lons) - min(lons) < 180

def convertToAcceptableLongitude(l):
    a = l % 360
    if a > 180:
        return a - 360
    else:
        return a
    
unableToResolveIDLthroughOffset = []
ableToResolveIDLthroughOffset = []

def addOffset(p, offset):
    for i in range(len(p)):
        p[i][1] = convertToAcceptableLongitude(p[i][1] + offset)

import copy

def getIDLSafeCoordsAndOffset(positionArray, clade = None):
    offset = 0
    posArrayOffset = copy.deepcopy(positionArray)
    while offset < 360 and not areLociSuitableForAveraging(posArrayOffset):
        addOffset(posArrayOffset, 1)
        offset += 1
    if offset == 360:        
        unableToResolveIDLthroughOffset.append((positionArray, clade))
        return None
    else:
        if offset != 0:
            ableToResolveIDLthroughOffset.append((positionArray, offset, posArrayOffset))
        return (posArrayOffset, offset)
    
def getAverageIDLSafe(positionArray, weights = None):
    if len(positionArray) == 0:
        return None
    else:
        (IDLSafeCoords, offset) = getIDLSafeCoordsAndOffset(positionArray)
        avg = getAverage(IDLSafeCoords, weights)
        avg[1] = convertToAcceptableLongitude(avg[1] + offset * -1)
        return avg
            
def getAverage(positionArray, weights = None):
    latsum = 0
    lonsum = 0
    for i in range(len(positionArray)):
        pos = positionArray[i]
        if weights is None:
            latsum += pos[0]
            lonsum += pos[1]
        else:
            latsum += pos[0] * weights[i]
            lonsum += pos[1] * weights[i]
    if len(positionArray) == 0:
        return None
    else:
        if weights is None:
            return [latsum / len(positionArray), lonsum / len(positionArray)]
        else:
            weightSum = sum(weights)
            return [latsum / weightSum, lonsum / weightSum]

exclusionDistRatioThreshold = 3
outlierMinAvgDistKm = 1000

allOutliers = {}

def getAverageLessOutliers(clade, positionArray, weights = None):
    if len(positionArray) > 2:
        outliers = getOutliers(positionArray, exclusionDistRatioThreshold, outlierMinAvgDistKm)
        outliers.reverse()
        posArrayCopy = positionArray.copy()
        weightsCopy = weights.copy()
        if len(outliers) > 0:
            allOutliers[clade] = []
        for outlier in outliers:
            posArrayCopy.pop(outlier)
            allOutliers[clade].append(positionArray[outlier])
            if weightsCopy is not None:
                weightsCopy.pop(outlier)
        return getAverageIDLSafe(posArrayCopy, weightsCopy)
    else:
        return getAverageIDLSafe(positionArray, weights)

def getAveragesFromCalcedAndBasal(children, clade, kitsMap, averages, tmrca):
    positionArray = []
    elapsedTimes = []
    for child in children:
        if child in averages:
            positionArray.append([averages[child][0], averages[child][1]])
            elapsedTimes.append(1 / max(tmrca[clade] - tmrca[child],500))
    if clade in kitsMap:
        for child in kitsMap[clade]:
            positionArray.append([child["latitude"], child["longitude"]])
            elapsedTimes.append(1 / max(tmrca[clade] - child["ybp"],500))
    return getAverageLessOutliers(clade, positionArray, elapsedTimes)

def calculateInitialAverages(clade, hierarchy, kitsMap, averages, tmrca):
    children = TreeOps.getChildren(clade, hierarchy)
    theAverage = None
    if len(children) > 0:
        for child in children:
            calculateInitialAverages(child, hierarchy, kitsMap, averages, tmrca)
        theAverage = getAveragesFromCalcedAndBasal(children, clade, kitsMap, averages, tmrca)
    else:
        theAverage = getAveragesFromCalcedAndBasal([], clade, kitsMap, averages, tmrca)
    if theAverage != None:
        averages[clade] = theAverage

def getClosestPointPoly(poly, point):
    pol_ext = LinearRing(poly.exterior.coords)
    d = pol_ext.project(point)
    p = pol_ext.interpolate(d)
    return list(p.coords)[0]

def refinePolyIDLSafe(positions, parentInitialAverage):
    (IDLSafeCoords, offset) = getIDLSafeCoordsAndOffset(positions + [parentInitialAverage])
    IDLSafeParentInitAvg = IDLSafeCoords[-1]
    IDLSafeCoords = IDLSafeCoords[0:-1]
    updated = refinePoly(IDLSafeCoords, IDLSafeParentInitAvg)
    updatedLon = convertToAcceptableLongitude(updated[1] + -1 * offset)
    return [updated[0],updatedLon]

def refinePoly(positions, parentInitialAverage):    
    gparentPoint = Point(parentInitialAverage)
    childPoly = Polygon(positions)
    isInside = childPoly.covers(gparentPoint)
    updated = []
    if isInside:
        updated = np.average([parentInitialAverage,getAverage(positions)],0)
    else:
        updated = getClosestPointPoly(childPoly, gparentPoint)
    return updated

import math


def getShortestGlobeSqDist(p1, p2):
    pdiff = abs(p2[1] - p1[1])
    pdiff = min(pdiff, 360 - pdiff)
    return getSqDist([p1[0],0], [p2[0],pdiff])

def getShortestGlobeDist(p1, p2):
    pdiff = abs(p2[1] - p1[1])
    pdiff = min(pdiff, 360 - pdiff)

    return getDist([p1[0],0], [p2[0],pdiff])
    
def getDist(p1, p2):
    thecos = math.cos((p1[0]-p2[0])/2 * math.pi / 180)
    return math.sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) * thecos)

def getPointDists(points):
    pointTotals = np.zeros(len(points))
    for i in range(len(points)):
        for j in range(len(points) - i - 1):
            secondIdx = j + i + 1
            dist = getShortestGlobeDist(points[i], points[secondIdx])
            pointTotals[i] += dist
            pointTotals[secondIdx] += dist
    return pointTotals

def getExclusionDistAvgs(pointTotals):
    exclusionDistAvgs = np.zeros(len(pointTotals))
    thesum = np.sum(pointTotals) / 2
    for i in range(len(pointTotals)):
        exclusionDistAvgs[i] = (thesum - pointTotals[i]) / (len(pointTotals) * len(pointTotals - 1) / 2 - (len(pointTotals) - 1))
    return exclusionDistAvgs

def getOutliers(points, relativeThreshold, minDistThresholdKm):
    minDistThreshold = minDistThresholdKm / 111 * (len(points) - 1)
    pointTotals = getPointDists(points)
    exD = getExclusionDistAvgs(pointTotals)
    outliers = []
    if len(points) == 3:
        relativeThreshold += 1
    for i in range(len(points)):
        if pointTotals[i]/(len(points)-1) / exD[i] > relativeThreshold and pointTotals[i] > minDistThreshold:
            outliers.append(i)
    return outliers
    
def getSqDist(p1, p2):
    thecos = math.cos((p1[0]-p2[0])/2 * math.pi / 180)
    thecos2 = thecos * thecos
    return (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) * thecos2

def refineLine(clade, p1, p2, parentInitialAverage):
    p1Dist = getShortestGlobeSqDist(p1, parentInitialAverage)
    p2Dist = getShortestGlobeSqDist(p2, parentInitialAverage)
    midPoint = getAverageIDLSafe([p1,p2])
    midPointDist = getShortestGlobeSqDist(midPoint, parentInitialAverage)
    if midPointDist < p1Dist and midPointDist < p2Dist:
        return midPoint
    else:
        closest = p1
        closestDist = p1Dist
        furthestDist = p2Dist
        if p1Dist > p2Dist:
            closest = p2
            closestDist = p2Dist
            furthestDist = p1Dist
        p1p2Dist = getShortestGlobeSqDist(p1,p2)
        if closestDist < midPointDist * 2:
            return closest
        if closestDist + p1p2Dist + .707 * math.sqrt(closestDist * p1p2Dist) < furthestDist:
            return closest
        else:
            return midPoint
                
def calculateSecondPassRefinement(clade, positions, parentInitialAverage):
    if len(positions) > 2:
        return refinePolyIDLSafe(positions, parentInitialAverage)
    else:
        if len(positions) == 2:
            return refineLine(clade, positions[0], positions[1], parentInitialAverage)
        else:
            return positions[0]

def filterOutliers(positions, exclusionDistRatioThreshold, outlierMinAvgDistKm):
    if len(positions) > 2:
        outliers = getOutliers(positions, exclusionDistRatioThreshold, outlierMinAvgDistKm)
        remaining = []
        for i in range(len(positions)):
            if i not in outliers:
                remaining.append(positions[i])
        return remaining
    else:
        return positions

def secondPass(clade, hierarchy, kitsMap, averages, refined):
    children = TreeOps.getChildren(clade, hierarchy)
    theRefined = None
    positionArray = []
    if len(children) > 0:
        for child in children:
            if child in averages:
                secondPass(child, hierarchy, kitsMap, averages, refined)
                positionArray.append(refined[child])
    if clade in kitsMap:
        for sample in kitsMap[clade]:
            positionArray.append([sample["latitude"], sample["longitude"]])
    if clade in hierarchy:
        parent = hierarchy[clade]
        theRefined = calculateSecondPassRefinement(clade, filterOutliers(positionArray,exclusionDistRatioThreshold,outlierMinAvgDistKm), averages[parent])
    else:
        theRefined = getAverageIDLSafe(positionArray)
    
    refined[clade] = theRefined

averages = {}

addlBasalAncientAdded = 0

def process():
    (root, hierarchy, ybp, tmrca) = TreeOps.parseTreeXMLTMRCA(getYFulltreeFileName("Adam"))

    ancientKits = downloadAndGetAncientKits()
    loadYFSamples(projects, kitsMap, ancientKits)
#    hypoAlaskan = {"id":"HypotheticalAlaskan", "latitude": 61.199838 , "longitude": -149.904093, "ybp": 10000}
#    kitsMap["Q-M1107"] = [hypoAlaskan]
#    hypoYakutianM930= {"id":"HypotheticalYakutian-M930", "latitude":  62.05776 , "longitude": 129.459328 , "ybp": 10000}
#    kitsMap["Q-M930"] = [hypoYakutianM930]
#    hypoYakutianM3 = {"id":"HypotheticalYakutian-M3", "latitude":  62.05776 , "longitude": 129.459328 , "ybp": 10000}
#    kitsMap["Q-M3"] = [hypoYakutianM3]

    calculateInitialAverages(root, hierarchy, kitsMap, averages, tmrca)
    refined = {}
    secondPass(root, hierarchy, kitsMap, averages, refined)
#    #secondRefined = {}
#    #secondPass(root, hierarchy, kitsMap, refined, secondRefined)
    separateProjectsIntoHgs(hierarchy, ybp, tmrca, kitsMap, projectsToSeparateIntoHGs, refined)   
    writeOutManifest(yfullVersion)     

def pruneTreeToSubtree(hierarchy, ybp, tmrca, subtreeRoot, refined, kitsMap, projectTips):
    prunedKitsMap = {}
    latLons = {}
    newHier = {}
    newYbp = {}
    newTmrca = {}
    newYbp[subtreeRoot] = ybp[subtreeRoot]
    newTmrca[subtreeRoot] = tmrca[subtreeRoot]
    projectTips[subtreeRoot] = {"ybp": ybp[subtreeRoot], "tmrca": tmrca[subtreeRoot]}
    if subtreeRoot in kitsMap:
        prunedKitsMap[subtreeRoot] = kitsMap[subtreeRoot]
        projectTips[subtreeRoot]["samples"] = kitsMap[subtreeRoot]
    if subtreeRoot in refined:
        latLons[subtreeRoot] = refined[subtreeRoot]
        projectTips[subtreeRoot]["latLon"] = refined[subtreeRoot]
    recursivelyPruneTreeToSubtree(hierarchy, ybp, tmrca, subtreeRoot, newHier, newYbp, newTmrca, latLons, refined, prunedKitsMap, kitsMap)
    TreeOps.writeXMLTMRCA_latlon(subtreeRoot, newHier, newYbp, newTmrca, latLons, getSubProjectTreeFileName(subtreeRoot))
    writePrunedSamplesMap(subtreeRoot, prunedKitsMap)

def kitsToStr(kits, innerDelim, outerDelim):
    kitStrings = []
    for kit in kits:
        theKit = [kit["id"], str(round(kit["latitude"],2)), str(round(kit["longitude"],2))]
        
        kitStrings.append(innerDelim.join(theKit))
    return outerDelim.join(kitStrings)

def writePrunedSamplesMap(project, kitsMap):
    with open(getSubProjectSamplesFileName(project), "w") as w:
        for clade in kitsMap:
            w.write(clade + "," + kitsToStr(kitsMap[clade], " ", ":") + "\n")
        w.close
        
def recursivelyPruneTreeToSubtree(hierarchy, ybp, tmrca, clade, newHier, newYbp, newTmrca, latLons, refined, prunedKitsMap, kitsMap):
    children = TreeOps.getChildren(clade, hierarchy)
    for child in children:
        recursivelyPruneTreeToSubtree(hierarchy, ybp, tmrca, child, newHier, newYbp, newTmrca, latLons, refined, prunedKitsMap, kitsMap)
        newHier[child] = clade
        hierarchy.pop(child)
        newYbp[child] = ybp[child]
        newTmrca[child] = tmrca[child]
        if child in kitsMap:
            prunedKitsMap[child] = kitsMap[child]
            kitsMap.pop(child)
        if child in refined:
            latLons[child] = refined[child]
            refined.pop(child)
        ybp.pop(child)
        tmrca.pop(child)

def createAdamProject(hierarchy, ybp, tmrca, kitsMap, refined, projectTips):
    prunedKitsMap = {}
    latLons = {}
    newHier = {}
    newYbp = {}
    newTmrca = {}
    subtreeRoot = "Adam"
    newYbp[subtreeRoot] = ybp[subtreeRoot]
    newTmrca[subtreeRoot] = tmrca[subtreeRoot]
    if subtreeRoot in kitsMap:
        prunedKitsMap[subtreeRoot] = kitsMap[subtreeRoot]
    if subtreeRoot in refined:
        latLons[subtreeRoot] = refined[subtreeRoot]
    recursivelyPruneTreeToSubtree(hierarchy, ybp, tmrca, subtreeRoot, newHier, newYbp, newTmrca, latLons, refined, prunedKitsMap, kitsMap)
    for k in projectTips:
        thetip = projectTips[k]
        newYbp[k] = thetip["ybp"]
        newTmrca[k] = thetip["tmrca"]
        if "samples" in thetip:
            prunedKitsMap[k] = thetip["samples"]
        if "latLon" in thetip:
            latLons[k] = thetip["latLon"]
    TreeOps.writeXMLTMRCA_latlon(subtreeRoot, newHier, newYbp, newTmrca, latLons, getSubProjectTreeFileName(subtreeRoot))
    writePrunedSamplesMap(subtreeRoot, prunedKitsMap)
    
def separateProjectsIntoHgs(hierarchy, ybp, tmrca, kitsMap, projectsToSeparateIntoHGs, refined):
    projectTips = {}
    for project in projectsToSeparateIntoHGs:
        pruneTreeToSubtree(hierarchy, ybp, tmrca, project, refined, kitsMap, projectTips)
    print(projectTips)
    createAdamProject(hierarchy, ybp, tmrca, kitsMap, refined, projectTips)
    
yfullVersion = "v7.05.00"


import datetime as dt
def writeOutManifest(yfullVersion):
    manifestFile = "C:\\PhyloGeographer\\2019 Refactor\\manifest"
    allAncient = 0
    allAdded = 0
    ancientIds = []
    for k in ancientAddedToProject:
        allAncient += len(ancientAddedToProject[k])
        for i in ancientAddedToProject[k]:
            ancientIds.append(i)
    for k in allAddedToProject:
        allAdded += len(allAddedToProject[k])
    
    with open(manifestFile, "w") as w:
        w.write("samples," + str(allAdded) + "\n")
        w.write("ancient," + str(allAncient) + "\n")
        w.write("YFull-version," + yfullVersion + "\n")
        w.write("allAncient," + ":".join(ancientIds) + "\n")
        w.write("computedDate," + str(dt.date.today().strftime("%b %d %Y")))
    w.close()
    
process()
