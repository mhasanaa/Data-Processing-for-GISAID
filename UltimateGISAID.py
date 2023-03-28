# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:54:59 2023

@author: Muhammad HASAN & KAN Chiming
"""
import argparse
import pandas as pd
import time
from datetime import datetime
import re

start_time = time.process_time()
current_time = datetime.now().strftime("%H:%M:%S")
print("Start =", current_time)
e = 0

def checkpoint(e):
    e = time.process_time() - start_time - e
    w = 'TIME TAKEN: ' + str(e) + 's'
    return w
    
def strip_n(f1): return [f1[i].rstrip('\r\n') for i in range(len(f1))]

def GISAIDProcessing(args):

    input_file = args.input_GISAID
    f0 = strip_n( open(input_file, errors="ignore").readlines() )
    print("Intial number from Raw Data %d Seq ..... %" %(len(f0)/2,checkpoint(e)) ) # how many input we have
    ID,Seq,ColDate,Date,Month,Year,EPI_ID= [],[],[],[],[],[],[]
    ColDate2,EPI_ID2,Seq2 = [],[],[]
    b = 0
    for i in np.arange(0,(len(f0)),2):
        S=f0[i+1]
        s= f0[i].replace('|', ' ')
        if ('X' not in S) and (len(S)>=1256): #Taken Length above 1256, only bottom limit
            Seq.append(S)
            w = f0[i].split("|")
            if "-" in w[2]:
                ColDate.append(w[2]) #collection Date
                EPI_ID.append(w[3])
            elif "-" in w[3]:
                ColDate.append(w[3]) #collection Date
                EPI_ID.append(w[4])
            b +=1
    
    print("Number of Sequences without X amino acid = %d Seq..... %" %(b,checkpoint(e)) )
    print("Number of Sequences with X amino acid = %d Seq..... %" %(len(f0)/2 - b,checkpoint(e)) )
    ErrorCount=0
    for i in range(len(ColDate)):
        w = ColDate[i].split("-")
        if ( len(w) >= 3 ):
            Date = w[2];month = w[1];year = w[0] #
        else:
            print('%d %s'%(i,ColDate[i]))
            ErrorCount+=1
            continue
        if month != "00" and int(year) >= 2019:
            if Date == "00": Date = "23"
            ColDate2.append(year+"-"+month+"-"+Date)
            EPI_ID2.append(EPI_ID[i])
            Seq2.append(Seq[i])
    
    df0 = pd.DataFrame({'ID':EPI_ID2,'ColDate': ColDate2,'sequence': Seq2})
    print("Length of Selected Data (after removing non-specific collection date): {}".format(len(df0)) )
    df0["ColDate"] = pd.to_datetime(df0["ColDate"])
    df0["Year"] = pd.DatetimeIndex(df0["ColDate"]).year
    df0["Month"] = pd.DatetimeIndex(df0["ColDate"]).month
    df0["MonthIndex"] = ( (df0["Year"]-2020)*12 ) + df0["Month"]
    df0.loc[(df0.MonthIndex == 0),'MonthIndex'] = 1
    df0['sequence'] = df0['sequence'].map(lambda x: x.rstrip('*'))
    df0.sort_values(by=['MonthIndex'],inplace=True, ascending=[True])
    GISAID_serial_number = re.findall(r'spikeprot\d+', args.input_GISAID)
    v = re.findall(r'\d+', GISAID_serial_number[0])
    
    if args.reference == None:
        print('There is no reference database...')
        print('Processsing all the unique sequences...')
        df0Final = df0.drop_duplicates(subset=['sequence'], keep="first")
        df0Final.to_excel("{}_{}_{}Seq.xlsx".format(args.output_name,v[0],len(df0)))
        
    else:
        print('The database directory :%s' % args.reference)
        COVIDMappingkeys=['sequence','mutation info|insertion info']
        df1 = pd.read_excel (args.reference,sheet_name ='Sheet1',usecols=COVIDMappingkeys)
        df0 = pd.merge (df0,df1,on="sequence",how="left").drop_duplicates (subset=['ID','sequence'], keep="first")
        df0Final = df0[df0['mutation info|insertion info'].isna()].drop_duplicates (subset=['sequence'], keep="first")
        df0Final.to_excel("{}_{}_{}Seq.xlsx".format(args.output_name,v[0],len(df0)))
    
    print(df0Final);print(checkpoint(e))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Data Processing for deLemus')
    parser.add_argument('-i', '--input', type=str, required=True,
                        dest='input_GISAID', help='Input file name, describe the directory if the input file in different directory (fasta)')
    parser.add_argument('-o', '--output_name', type=str, default='GISAID-Unique',
                        dest='output_name', help='Output filename without .xlsx (excel)')
    parser.add_argument('-r', '--reference', type=str, default=None,
                        dest='reference', help='Excel database path of previously collected sequences (excel). \
                        Use this whenever you want to update the data, so the code will only work on the new sequences only.\
                        Default is None')
    args = parser.parse_args()
    print('The input directory :%s' % args.input_GISAID)
    print('Processing GISAID Sequences...')
    GISAIDProcessing(args)