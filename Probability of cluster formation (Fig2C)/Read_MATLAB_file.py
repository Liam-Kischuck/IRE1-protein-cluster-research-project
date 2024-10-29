# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 16:24:29 2022

@author: Michael
"""

def read_matlab_data(filename):
    """filename is a .txt file. The data is in rows and can end with stuff like ,,,,\n. 
    This function deletes the trailing stuff."""
    count = 0
    record_timelines = []
    while True: 
        # Get next line from file
        a_timeline = filename.readline()
        # if line is empty
        # end of file is reached
        if not a_timeline:
            break
        #Otherwise continue
        record_timelines.append((a_timeline.split(',')))
        for i in range(len(record_timelines[count])):
            if len(record_timelines[count][i]) > 0:
                if record_timelines[count][i][0].isdigit() == True:
                    record_timelines[count][i] = float(record_timelines[count][i])
                else:
                    record_timelines[count] = record_timelines[count][:i]
                    break
            else:
                record_timelines[count] = record_timelines[count][:i]
                break
        count += 1
    return record_timelines
