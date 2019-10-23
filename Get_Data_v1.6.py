# -*- coding: utf-8 -*-
"""
Created on Friday Oct 4 12:20:12 2019

@author: DuchJ
"""

import pandas as pd
import numpy as np
import re
import os

#function that deletes the object in the index of the array and creates a new array with a copy of that object
#this is essentially append for numpy arrays
def poprow(my_array,pr):
    """ row popping in numpy arrays
    Input: my_array - NumPy array, pr: row index to pop out
    Output: [new_array,popped_row] """
    i = pr
    pop = my_array[i]
    pop = np.array(pop)

    pop = np.array(pop).astype(int)

    new_array = np.delete(my_array,i)
    return [new_array,pop]

#function which deletes the top values in the array to a certain number, 
#but saves the deleted ones in a separate array using the poprow function
#this is to delete the HAZ hardness indents that landed in the FZ, 
#but we're gonna append the FZ indents to the FZ array later ;)
def del_FZ(hardness_array, not_HAZ_number):
    """input: the hardness array and the number of points we're deleting
    output: [the new hardness array, the deleted points or FZ array]"""
    FZ_array = []
    FZ_array = np.asarray (FZ_array, dtype=int)
    
    #not_HAZ_number = int(not_HAZ_number)
    #FZ_array = np.empty(shape=not_HAZ_number, dtype=int)

    hardness_array_cut = hardness_array

    x=0
    for x in range(not_HAZ_number):
        [hardness_array_cut, FZ_pop] = poprow(hardness_array_cut, 0)
        FZ_array = np.append(FZ_array, FZ_pop)

    return [hardness_array_cut, FZ_array]

#function which cleans up the data from the sheet, renames the columns, and returns the hardness array 
def trim_fat (hard_sheet):
    """Input: the excel spreadsheet from pandas, and the number of points that landed in the FZ (int)
    Output: the hardness in a numpy array, and the FZ hardness points also in a numpy array"""
    hard_sheet_length = len(hard_sheet.index)
    hard_sheet_new = hard_sheet.drop(
            [0,1,2,(hard_sheet_length-1),(hard_sheet_length-2),(hard_sheet_length-3)]
            )
    hard_sheet_new = hard_sheet_new.rename(
        columns={
            "":"index",
            "Name":"Mean Length",
            "Number of Included Replicates":"Hardness",
        })

    hardness_array =  hard_sheet_new.Hardness.to_numpy()
    
    return hardness_array

#this function takes the excel file name, opens the file, trims the data, removes the FZ from the HAZ and appends it to the FZ array
def get_hard_740H (filename, not_HAZ_number):
    """
    xlsx = pd.ExcelFile(filename)
    hardness_sheets = []
    for sheet in xlsx.sheet_names:
        hardness_sheets.append(xlsx.parse(sheet))
    """
    
    HAZ= pd.read_excel(filename, sheet_name=0)
    FZ= pd.read_excel(filename, sheet_name=1)
    
    HAZ_trim = trim_fat(HAZ)
    [HAZ_array, FZ_in_HAZ] = del_FZ(HAZ_trim, not_HAZ_number)
    
    
    FZ_trim = trim_fat(FZ)
    FZ_array = np.append(FZ_trim, FZ_in_HAZ)
    
    return [HAZ_array, FZ_array]

def get_hard_347H (filename, meta_data):
    #for 347H order is [material, vector_ID, time, time_type, temperature, constant, sheet_order, sample_info]
    #sheet_order is ['haz', 'fz', 'bm'] sample_info is a list with 
    
    '''
    funtion takes the filename of the excel document and the meta_data as input
    output is a list of lists (still trying to figure out best way to do this)
    of the order [material, temperature, time, strain, trace ID, hardness vector]
    '''
    
    sheet_number = how_many_sheets_347 (meta_data)
    sheet = 0
    haz_count = 0
    strain_count = 0
    
    #this is a preliminary storage device, i'm not conviced that a list of lists
    #is the best way to store this data, ideally we would run this once and just get
    #one document to save and call for the data analysis
    hardness_array = []
    
    #It appears that a dataframe would be the best method to store this data
    #we're putting the dataframe into another function because i love functions so much
    
    for sheet in range(sheet_number):
        vector_type = meta_data[6][sheet%3] #this will select haz, fz, or bm
        
        #want to increase by 1 every three counts
        if sheet !=0 and sheet%3 == 0:
            strain_count +=3
        
        if meta_data[1] == 'z' and sheet == 10:
            vector_type = 'fz'
            
            hardness_vector= pd.read_excel(filename, sheet_name=sheet)
            trimmed_vector = trim_fat(hardness_vector)
            
            #then we add FZ to datasheet and break, it's the end of the loop
            hardness_array.append(
                    [meta_data[0], meta_data[4], meta_data[2], meta_data[3], 
                     meta_data[7][strain_count], 'FZ', trimmed_vector]
                    )
            
            #return [material, temperature, time, time type, strain, trace ID, hardness vector]
            
            break
        
        hardness_vector= pd.read_excel(filename, sheet_name=sheet)
        trimmed_vector = trim_fat(hardness_vector)
        
        if vector_type == 'haz': #0 means a HAZ trace, but 
            #if there is a problem here, lets look at meta_347 and using += instead of append
            not_HAZ_number = int(meta_data[7][1 + haz_count])
            [HAZ_array, FZ_in_HAZ] = del_FZ(trimmed_vector, not_HAZ_number)
            haz_count +=3 #if sample info is constant, FZ number, sheet number, 
            #we need to progress 3 spaces to the next FZ number next loop
            
            hardness_array.append(
                    [meta_data[0], meta_data[4], meta_data[2], meta_data[3],
                     meta_data[7][strain_count], 'HAZ', HAZ_array]
                    )
            
        elif vector_type == 'fz':
            #add the trimmed vector to the FZ in HAZ from the previous loop
            FZ_array = np.append(trimmed_vector, FZ_in_HAZ)
            #then we add the FZ to the datasheet
            
            hardness_array.append(
                    [meta_data[0], meta_data[4], meta_data[2], meta_data[3], 
                     meta_data[7][strain_count], 'FZ', FZ_array]
                    )
            
        else: #it is bm otherwise
            BM_array = trimmed_vector
            hardness_array.append(
                    [meta_data[0], meta_data[4], meta_data[2], meta_data[3], 
                     meta_data[7][strain_count], 'BM', BM_array]
                    )
           
    return hardness_array
        
        #return [material, temperature, time, strain, trace ID, hardness vector]
        # if vector_ID == 'a' sample info is just FZ number
        # if vector_ID == 'c' order is [strain, FZ_number, first_vector]

#assuming the material is 740H, the filename should have all the meta data
def meta_740H(filename):
    mlist = re.split("( )",filename)
    material = mlist[0]
    temperature = int(mlist[2])
    strain = int(mlist[6])
    FZ_number = 2 * int(mlist[18])
    #time= int(mlist[26])
    
    if mlist[26] == '2k':
        time = 2000.0
    elif mlist[26] == '4k':
        time = 4000.0
    elif mlist[26] == '8k':
        time = 8000.0
    elif mlist[26] == '1k':
        time = 1000.0
    else:
        time = float(mlist[26])
    
    time_type = (re.split("\\.",mlist[-1]))[0]
    
    return [material, time, time_type, temperature, strain, FZ_number]

def meta_347H(filename):
    mlist = re.split("( )", filename)
    
    vector_ID = mlist [0]
    if vector_ID == 'a': #this was the standard before c, one sample, one condition etc
        material = mlist[2] #its 347H, shouldn't be different
        composition = mlist[4] #either 0.08c or 0.04c
        temperature = int(mlist[6]) #either 600 or 675
        strain = float(mlist[10]) #either 0, 5, 10, 15, or 20 % strain
        width = int(mlist[34]) #width of hardness trace
        FZ_number = width * int(mlist[22])
        time = int(mlist [38])
        
        if mlist[38] == '2k':
            time = 2000.0
        elif mlist[38] == '4k':
            time = 4000.0
        elif mlist[38] == '8k':
            time = 8000.0
        elif mlist[38] == '1k':
            time = 1000.0
        
        time_type = re.split("()", mlist[40])[2] + re.split("()", mlist[40])[4]
        sheet_order = ['haz', 'fz', 'bm']
        sample_info = [strain, FZ_number]
        
        constant = ['composition and strain', composition, strain]
        
        return [material, vector_ID, time, time_type, temperature, constant, \
                sheet_order, sample_info,width]

    elif vector_ID == 'b': #meta data 'b' means we have two 0 pct samples in the same sample
        material = '347H'
        
        time = int(mlist [2])
        
        if mlist[2] == '2k':
            time = 2000.0
        elif mlist[2] == '4k':
            time = 4000.0
        elif mlist[2] == '8k':
            time = 8000.0
        elif mlist[2] == '1k':
            time = 1000.0
        
        time_type = mlist[4]
        strain = int(mlist[6])
        temp = re.split("()", mlist[10])
        temperature = int(temp[2]+temp[4]+temp[6])
        sheet_1 = mlist[14]
        sheet_2 = mlist[16]
        sheet_3 = mlist[18]
        sheet_order = [sheet_1, sheet_2, sheet_3]
        
        composition_A = mlist[20]
        FZ_number_A = 3 * int(re.split("()", mlist[24])[2])
        A_first = int(re.split("()", mlist[22])[4])
        A_info = [composition_A, FZ_number_A, A_first]
        
        composition_B = mlist[26]
        FZ_number_B = 3 * int(re.split("()", mlist[30])[2])
        B_first = int(re.split("()", mlist[28])[4])
        B_info = [composition_B, FZ_number_B, B_first]
        
        sample_info = A_info + B_info #sample info will be a list with info in order
        
        constant = ['strain', strain]
        
        return [material, vector_ID, time, time_type, temperature, constant, \
                sheet_order, sample_info,2]
        
    elif vector_ID == 'c' or 'z': #meta data 'c' means we have all four samples in the same sample, z is a mistake at the end (experimental)
        material = '347H'
        #time = int(mlist[2])
        
        if mlist[2] == '2k':
            time = 2000.0
        elif mlist[2] == '4k':
            time = 4000.0
        elif mlist[2] == '8k':
            time = 8000.0
        else:
            time = float(mlist[2])
        
        time_type = mlist[4]
        composition = int(mlist[6])
        temp = re.split("()", mlist[10])
        temperature = int(temp[2]+temp[4]+temp[6])
        constant = ['composition', composition]
        
        sheet_1 = mlist[14]
        sheet_2 = mlist[16]
        sheet_3 = mlist[18]
        sheet_order = [sheet_1, sheet_2, sheet_3] 
        #the sheet_order is labels for the types of hardness traces
        #for example, the most common one will be: haz, fz, bm
        #so far there are no other orders, but this is to make it a bit easier if it changes in the future
        
        counting = 20
        
        sample_info = []

        #why four? there are four sets of hardness traces (order haz, fz, bm)
        #for each set of 3 hardness traces, we will get the strain value, strain
        #then we will get the number for the first vector, which sheet is it on?, first_vector
        #then we will get FZ number
        #if we have a vector ID z, that means that the fourth set is actually just a FZ trace
        for x in range (4):
            
            strain = re.split("()", mlist[counting])[2]+re.split("()", \
                              mlist[counting])[4]
            
            #if not strain.isdigit(): #if it is 5pct, strain above will return strain '5p', but 10pct will return '10'
                #strain = re.split("()", mlist[counting])[2]
            strain = int(strain)
            
            first_vector = re.split("()", mlist[counting+2])[4]+re.split("()", \
                                    mlist[counting+2])[6]
            if not first_vector.isdigit(): #if it is not starting with vector 10, this will run
                first_vector = re.split("()", mlist[counting+2])[4]
            first_vector = int(first_vector)
            
            if vector_ID=='z' and x==3:
                sample_info += [strain, 0, first_vector]
                counting += 6
            else:
                FZ_number = 3 * int(re.split("()", mlist[counting+4])[2])
                sample_info += [strain, FZ_number, first_vector]
                counting += 6
                
        return [material, vector_ID, time, time_type, temperature, constant, \
            sheet_order, sample_info,2]
        
    else:
        print ('you should debug') #this shouldn't happen, so there is a spelling mistake or code error somewhere

def get_meta(test_file):
    if re.split("( )", test_file)[0] == '740H':
        meta_data = meta_740H(test_file)
    else:
        meta_data = meta_347H(test_file)
        
        '''
        if meta_data [1] == 'a':
            print (meta_data)
        #'''
        
    return meta_data

#this function takes the vector ID and returns the number of sheets in the excel file
def how_many_sheets_347 (meta_data):
    vector_ID = meta_data[1]
    sheet_number = 2
    if vector_ID == 'a':
        sheet_number = 3
    elif vector_ID == 'b':
        sheet_number = 6
    elif vector_ID == 'z':
        sheet_number = 10
    else:
        sheet_number = 12
    return sheet_number

#this function returns the composition of the 347H (low or high). I messed up, so this is the easiest way
#NOTE: this only works for vector ID of a
def get_composition (meta_data):
    vector_ID = meta_data [1]
    
    if vector_ID == 'b':
        print ('You cannot use get_composition for a b type vector')
    else:
        composition = meta_data[5][1]
        if composition=='0.04c' or composition == 4:
            composition_347 = float('0.04')
        elif composition == '0.08c' or composition == 8:
            composition_347 = float('0.08')
            
        return composition_347

#this function will take the get_hard data and meta data and make it into a dataframe!
#'''
def make_data_frame_347 (hardness_array, meta_data, hardness_data):
    
    hardness_data = pd.DataFrame(
            columns=['Material','Composition','Temperature','Time','Strain',
                     'Type','Hardness','Width']
                )
    
    #lets start by getting our metadata again 
    vector_ID = meta_data[1]
    width=meta_data[-1]
    
    '''
    if vector_ID == 'a':
        print (hardness_array)
    
    #'''
    b_type_count=0
    
    for i in range(len(hardness_array)):
        material = hardness_array[i][0]
        
        #get the composition, which i kinda screwed up on earlier but this is fine
        if vector_ID != 'b':
            composition = get_composition(meta_data) 
            #for a, c, and z types composition is right in the metadata!
            strain = int(hardness_array[i][4])
            
        else: 
            #otherwise we have to go to the sample_info part of the metadata
            
            if i !=0 and i%3 == 0:
                b_type_count +=3
                
            composition = meta_data[7][b_type_count]
            strain = meta_data[5][1]
                    
            #[material, vector_ID, time, time_type, temperature, constant, sheet_order, sample_info]
            
            #sample info is labeled composition_A, FZ_number_A, A_first, composition_b... 
            #thus we want composition_A for the HAZ, FZ, and BM then after 3 
            #go to composition_b. every third file we increment this by 3 to 
            #get from composition a to b (three spaces)
            
        
        #print(strain)
                
        if composition in ('0.04C', '0.04c', '0.4C'):
            composition = 0.04
        elif composition in ('0.08C', '0.08c', '0.8C'):
            composition = 0.08
        else:
            composition = float(composition)
                
        temperature = int(hardness_array[i][1])
        
        if temperature == 670:
            temperature = 675
        
        if hardness_array[i][3]=='hr':
            time = float(hardness_array[i][2])
        elif hardness_array[i][3] =='min':
            time = float(hardness_array[i][2])/60.0
            
        hardness_type = hardness_array[i][5]        
        
        hardness = hardness_array[i][6]
        
        #pd.DataFrame({"a":[1, 2, 3, 4], "b":[5, 6, 7, 8]})
        
        hardness_data_2 = pd.DataFrame({
            "Material":[material],"Composition":[composition],
            "Temperature":[temperature],"Time":[time],"Strain":[strain],
            "Type":[hardness_type],"Hardness":[hardness],"Width":[width]
                })
        
        hardness_data = hardness_data.append(hardness_data_2, ignore_index = True)
    
    return hardness_data

def make_data_frame_740H (hardness_array, meta_data, hardness_data):
   
    #(material, temperature, strain, FZ_number, time, time_type)
    
    HAZ_hardness = hardness_array[0]
    FZ_hardness = hardness_array[1]
    #meta data order is [material, time, time_type, temperature, strain, FZ_number]
    
    material = '740H'
    temperature = meta_data[3]
    
    #reminder that meta data is [material, time, time_type, temperature, strain, FZ_number]
    
    if meta_data[2]=='hr':
        time = float(meta_data[1])
    elif meta_data[2] =='min':
        time = float(meta_data[1])/60.0
        
    if time==0:
        width=1
    else:
        width=2
    
    strain = meta_data[4]
    
    for i in range (2):
        if i ==0:
            hardness_type = 'HAZ'
            hardness = HAZ_hardness
            
        else:
            hardness_type = 'FZ'
            hardness = FZ_hardness
        
        hardness_data_2 = pd.DataFrame({
                "Material":[material],"Temperature":[temperature],
                "Time":[time],"Strain":[strain], "Type":[hardness_type],
                "Hardness":[hardness], "Width":[width]
                })
        
        hardness_data = hardness_data.append(hardness_data_2, ignore_index = True)
    
    return hardness_data

#''' Great! Now we have a data frame with all the information from one excel 
    #file, now lets do it with all the excel files in a directory!
def fill_data_frame(path):
    
    hardness_data = pd.DataFrame(
            columns=['Material','Composition','Temperature','Time','Strain',
                     'Type','Hardness','Width']
                )

    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if '.xlsx' in file:
                files.append(os.path.join(r, file))

    for f in files:
        
        hardness_data_i = pd.DataFrame(
            columns=['Material','Composition','Temperature','Time','Strain',
                     'Type','Hardness','Width']
                )
        
        f_title= f.split('\\')
        f_title = f_title[-1]
        
        print (f_title)
        
        #'''
        meta_data = get_meta(f_title)
        
        material_ID = meta_data[0]
        
        if material_ID == '347H' or material_ID == '347h':
            hardness_array = get_hard_347H (f, meta_data)
            hardness_data_i = make_data_frame_347 (hardness_array, meta_data, hardness_data_i)
        else:
            hardness_array = get_hard_740H (f, int(meta_data[5]))
            hardness_data_i = make_data_frame_740H (hardness_array, meta_data, hardness_data_i)
            
        hardness_data = hardness_data.append(hardness_data_i, ignore_index = True)
    
    hardness_data.to_csv(r'C:\Users\DuchJ\Desktop\Machine_Learning_Class\Final_Project\Final_Project_Code\all_data.csv')
    #return hardness_data
    #'''
    

'''    
#test_file = 'a 347H 0.08c 675 c 20 pct 0 v1 is fz 3 down into fz v1 is 3 wide 100 hr.xlsx'    
#test_file= 'b 500 hr 0 pctC 600degc order haz fz bm 0.4C v1to3 2fz 0.8C v4to6 3fz.xlsx'
test_file = 'a 347H 0.04c 600 c 5 pct 0 v1 is fz 4 down into fz v1 is 4 wide 100 hr.xlsx'
meta_data = get_meta(test_file)

hardness_array = get_hard_347H (test_file, meta_data)

#hardness_frame = make_data_frame_347(hardness_array, meta_data)
'''

path = r'C:\Users\DuchJ\Desktop\Machine_Learning_Class\Final_Project\Final_Project_Code\All_Data'
#path = r'G:\My Drive\347 and 740H Project\Hardness\Standard MetaData All Material'

hardness_data = fill_data_frame (path)


#eventually we want to pickle this bb
#df.to_pickle(file_name)