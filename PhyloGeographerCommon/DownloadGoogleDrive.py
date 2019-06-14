# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 09:16:47 2018

@author: hunte
"""
from googleapiclient.discovery import build
from httplib2 import Http
from oauth2client import file, client, tools
from googleapiclient.http import MediaIoBaseDownload

SCOPES = 'https://www.googleapis.com/auth/drive.readonly'

store = file.Storage('token.json')
creds = store.get()
if not creds or creds.invalid:
    flow = client.flow_from_clientsecrets('credentials.json',SCOPES)
    creds = tools.run_flow(flow, store)

print (creds)

service = build('sheets', 'v4', http=creds.authorize(Http()))

# Call the Sheets API
SPREADSHEET_ID = '1OttVBg8z-EIVKTG-PTo-XvzSEdqxKx8Uf2p6YjSe5jo'

def writeSupplementalKitsFile(kits, fil):
        headers = ['id','geocode','latitude','longitude','clade','negatives','ybp']
        with open(fil, 'w') as f:
            f.write(",".join(headers))
            f.write("\n")
            for kit in kits:
                print(kit)
                kitArr = [kit["id"],kit["geocode"],str(kit["latitude"]),str(kit["longitude"]),kit["clade"],kit["negatives"],str(kit["ybp"])]
                f.write(",".join(kitArr))
                f.write("\n")
        f.close()

def writeAncientSamplesFile(kits,fil):
    headers = ["id","latitude","longitude","clade","negatives","ybp","hg"]
    with open(fil, 'w') as f:
        f.write(",".join(headers))
        f.write("\n")
        for kit in kits:
            theid = "ANCIENT_" + kit["ID"]
            if kit["YFull-ID"] != "":
                theid = kit["YFull-ID"]
            negString = ""
            if str(kit["Negative-SNPs"]) != "":
                addedHgs = str(kit["Negative-SNPs"]).split(":")
                for i in range(len(addedHgs)):
                    addedHgs[i] = kit["Y-Haplogroup"][0] + "-" + addedHgs[i]
                negString = ":".join(addedHgs)
            kitArr = [theid,str(kit["Lat"]),str(kit["Long"]),kit["YFull-Branch"], negString,str(kit["Year"]),kit["Y-Haplogroup"]]
            f.write(",".join(kitArr))
            f.write("\n")
    f.close()
        
def downloadAncient(haplogroup, spreadsheetId, outfile, lines = 10000):
    SAMPLE_RANGE_NAME = "A1:P" + str(lines)

    result = service.spreadsheets().values().get(spreadsheetId=spreadsheetId, range=SAMPLE_RANGE_NAME).execute()
    values = result.get('values', [])
    
    def parseAlphabet(string):
        parsed = ""
        for a in string:
            b = ord(a)
            if (b > 96 and b < 123) or (b > 64 and b < 91) or b == 45 or b == 32 or b == 91 or b == 93 or (b > 47 and b < 58):
                parsed = parsed + a
        return parsed
    
    if not values:
        print('No data found.')
    else:
        kits = []
        count = 0
        for row in values:
            # Print columns A and E, which correspond to indices 0 and 4.
            #print(row)
            if count > 0:
                kit = {}
                kit["ID"] = row[0]
                kit["YFull-ID"] = row[3]
                kit["Lat"] = row[4]
                kit["Long"] = row[5]
                kit["Y-Haplogroup"] = row[6]
                print(row)
                kit["Year"] = row[15]
                kit["YFull-Branch"] = row[7]
                kit["Negative-SNPs"] = row[10]
                
                #print(kit)
                kits.append(kit)
            count += 1
        writeAncientSamplesFile(kits, outfile)
        
        

def downloadSupplemental(haplogroup, spreadsheetId, outfile, lines = 10000):
    SAMPLE_RANGE_NAME = haplogroup + "!A1:G" + str(lines)

    result = service.spreadsheets().values().get(spreadsheetId=spreadsheetId, range=SAMPLE_RANGE_NAME).execute()
    values = result.get('values', [])
    
    def parseAlphabet(string):
        parsed = ""
        for a in string:
            b = ord(a)
            if (b > 96 and b < 123) or (b > 64 and b < 91) or b == 45 or b == 32 or b == 91 or b == 93 or (b > 47 and b < 58):
                parsed = parsed + a
        return parsed
    
    if not values:
        print('No data found.')
    else:
        print('Name, Major:')
        kits = []
        count = 0
        for row in values:
            # Print columns A and E, which correspond to indices 0 and 4.
            print(row)
            if count > 0:
                kit = {}
                kit["id"] = row[0]
                kit["geocode"] = parseAlphabet(row[1])
                kit["latitude"] = ""
                kit["longitude"] = ""
                kit["clade"] = ""
                kit["negatives"] = ""
                kit["ybp"] = ""
                if len(row) > 2:            
                    kit["latitude"] = row[2]
                    kit["longitude"] = row[3]
                if len(row) > 4:            
                    kit["clade"] = row[4]
                if len(row) > 5:            
                    kit["negatives"] = row[5]
                if len(row) > 6:
                    kit["ybp"] = row[6]
                print(kit)
                kits.append(kit)
            count += 1
        writeSupplementalKitsFile(kits, outfile)
        
            #print('%s, %s, %s, %s, %s' % (row[0], row[1], row[2], row[3], row[4]))
def getSupplementalKitsFilename(haplogroup):
    return "C:\\PhyloGeographer\\preprocessing\\data\\supplemental\\" + haplogroup + "-supplemental-kits.txt"
  
#downloadSupplemental("J-M241",SPREADSHEET_ID, getSupplementalKitsFilename("J-M241")) 
