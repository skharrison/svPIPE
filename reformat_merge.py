import pandas as pd
import io
import sys 

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def reformat_breakdancer(file):
    unformat_break = pd.read_csv(file, sep='\t', skiprows=[0,1,2,3])
    info_list = []
    for index,row in unformat_break.iterrows():
        info_line = "SVLEN=" + str(row["Size"]) + ":END=" + str(row["Pos2"]) + ":num_reads=" + str(row["num_Reads"])
        info_list.append(info_line)
    new_breakdancer = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO'])
    chrom_list = unformat_break['#Chr1'].to_list()
    pos_list = unformat_break['Pos1'].to_list()
    alt_list = unformat_break['Type'].to_list()
    qual_list = unformat_break['Score'].to_list()
    new_breakdancer['#CHROM'] = chrom_list
    new_breakdancer['POS'] = pos_list
    new_breakdancer['ID'] = "breakdancer"
    new_breakdancer['REF'] = "N"
    new_breakdancer['ALT'] = alt_list 
    new_breakdancer['QUAL'] = qual_list 
    new_breakdancer['FILTER'] = "."
    new_breakdancer['INFO'] = info_list
    return new_breakdancer 


def reformat_ss(manta,caller):  
    end_list = [] 
    if caller == 'breakdancer':
        for index,row in manta.iterrows():
            if row['ALT'] == 'CTX':
                manta.drop(index, inplace=True)
            else:
                end = row['INFO'].split('END=')[1].split(":")[0]
                end_list.append(end)
        chroms = manta['CHROM'].to_list()
        start_pos = manta['POS'].to_list()
        svtype_list = manta['ALT'].to_list()
        qual = manta['QUAL'].to_list()
        new_fr = pd.DataFrame(columns=['CHROM', 'START_POS', 'STOP_POS', 'SV_TYPE', 'QUAL', 'CALLER'])
        new_fr['CHROM'] = chroms
        new_fr['START_POS'] = start_pos
        new_fr['STOP_POS'] = end_list
        new_fr['SV_TYPE'] = svtype_list
        new_fr['QUAL'] = qual
        new_fr['CALLER'] = caller
        
    elif caller == 'pindel':
        for index,row in manta.iterrows():
            SV_type = row['INFO'].split('SVTYPE=')[1].split(';')[0]
            end = row['INFO'].split('END=')[1].split(";")[0]
            end_list.append(end)
            new = SV_type 
            manta.at[index,'ALT'] = new
        chroms = manta['CHROM'].to_list()
        start_pos = manta['POS'].to_list()
        svtype_list = manta['ALT'].to_list()
        new_fr = pd.DataFrame(columns=['CHROM', 'START_POS', 'STOP_POS', 'SV_TYPE', 'QUAL', 'CALLER'])
        new_fr['CHROM'] = chroms
        new_fr['START_POS'] = start_pos
        new_fr['STOP_POS'] = end_list
        new_fr['SV_TYPE'] = svtype_list
        new_fr['QUAL'] = "."
        new_fr['CALLER'] = caller
        
    elif caller == 'delly':
        ind_l = []
        for index,row in manta.iterrows():
            gt = str(manta.iloc[index,9])
            gt_l = gt.split(":")[0]
            if gt_l == '0/0':
                ind_l.append(index)
        for index,row in manta.iterrows():
            SV_type = row['INFO'].split('SVTYPE=')[1].split(';')[0]
            if SV_type == "BND":
                end_list.append("BND")
            else:
                end = row['INFO'].split('END=')[1].split(";")[0]
                end_list.append(end)
            new = SV_type
            manta.at[index,'ALT'] = new
        chroms = manta['CHROM'].to_list()
        start_pos = manta['POS'].to_list()
        svtype_list = manta['ALT'].to_list()
        qual = manta['QUAL'].to_list()
        new_fr = pd.DataFrame(columns=['CHROM', 'START_POS', 'STOP_POS', 'SV_TYPE', 'QUAL', 'CALLER'])
        new_fr['CHROM'] = chroms
        new_fr['START_POS'] = start_pos
        new_fr['STOP_POS'] = end_list
        new_fr['SV_TYPE'] = svtype_list
        new_fr['QUAL'] = qual
        new_fr['CALLER'] = caller
        for num in ind_l:
            new_fr.drop(num, inplace=True)
    else:

        for index,row in manta.iterrows():
            SV_type = row['INFO'].split('SVTYPE=')[1].split(';')[0]
            if SV_type == "BND":
                end_list.append("BND")
            else:
                end = row['INFO'].split('END=')[1].split(";")[0]
                end_list.append(end)
            new = SV_type
            manta.at[index,'ALT'] = new
        chroms = manta['CHROM'].to_list()
        start_pos = manta['POS'].to_list()
        svtype_list = manta['ALT'].to_list()
        qual = manta['QUAL'].to_list()
        new_fr = pd.DataFrame(columns=['CHROM', 'START_POS', 'STOP_POS', 'SV_TYPE', 'QUAL', 'CALLER'])
        new_fr['CHROM'] = chroms
        new_fr['START_POS'] = start_pos
        new_fr['STOP_POS'] = end_list
        new_fr['SV_TYPE'] = svtype_list
        new_fr['QUAL'] = qual
        new_fr['CALLER'] = caller
               
    return new_fr

def combine_SV(breakdancer_file,pindel_file,delly_file,manta_file,lumpy_file,gridss_file):
    breakdancer = pd.read_csv(breakdancer_file,sep='\t')
    pindel = pd.read_csv(pindel_file,sep='\t')
    lumpy = pd.read_csv(lumpy_file,sep='\t')
    delly = pd.read_csv(delly_file,sep='\t')
    manta = pd.read_csv(manta_file,sep='\t')
    gridss = pd.read_csv(gridss_file,sep='\t')
    all_types = pd.concat([breakdancer,pindel,lumpy,delly,manta,gridss], axis=0)
    return all_types

def make_merged(file, max_size, min_size, percent_threshold, caller_threshold):
    sample_SV = pd.read_csv(file,sep='\t')
    for index, row in sample_SV.iterrows():
        if row['STOP_POS'] == "BND":
            sample_SV.drop(index, inplace=True)
    for index, row in sample_SV.iterrows(): 
        sv_len = int(row['STOP_POS']) - int(row['START_POS'])
        if int(sv_len) > int(max_size) or int(sv_len) < int(min_size):
            sample_SV.drop(index, inplace=True)
    all_overlaps = []
    min_max = []
    for r1 in sample_SV.itertuples():
        pos_overlap = ""
        smallest = float(r1.START_POS)
        largest = float(r1.STOP_POS)
        for r2 in sample_SV.itertuples():  
            if r1.CHROM == r2.CHROM:
                if float(r2.START_POS) <= float(r1.START_POS) and float(r2.STOP_POS) < float(r1.STOP_POS) and float(r2.STOP_POS) > float(r1.START_POS):      
                    caller = r2.CALLER
                    start2 = float(r2.START_POS)
                    stop2 = float(r2.STOP_POS)
                    start1 = float(r1.START_POS)
                    stop1 = float(r1.STOP_POS)
                    sv_type = r2.SV_TYPE
                    overlap = (stop2 - start1) / (stop1 - start1) * 100
                    overlap_2 = (stop2 - start1) / (stop2 - start2) * 100
                    if overlap >= float(percent_threshold) and overlap_2 >= float(percent_threshold):
                        info = "ID=" + caller + ':' + "START=" + str(start2) + ':' + "END=" + str(stop2) + ':' + sv_type + ';' 
                        pos_overlap += info
                        if start2 < smallest:
                            smallest = start2
                        elif stop2 > largest:
                            largest = stop2
                elif float(r2.START_POS) > float(r1.START_POS) and float(r2.STOP_POS) >= float(r1.STOP_POS) and float(r2.START_POS) < float(r1.STOP_POS):
                    caller = r2.CALLER
                    start2 = float(r2.START_POS)
                    stop2 = float(r2.STOP_POS)
                    start1 = float(r1.START_POS)
                    stop1 = float(r1.STOP_POS)
                    sv_type = r2.SV_TYPE
                    overlap = (stop1 - start2) / (stop1 - start1) * 100
                    overlap_2 = (stop1 - start2) / (stop2 - start2) * 100
                    if overlap >= float(percent_threshold) and overlap_2 >= float(percent_threshold):
                        info = "ID=" + caller + ':' + "START=" + str(start2) + ':' + "END=" + str(stop2) + ':' + sv_type + ';' 
                        pos_overlap += info 
                        if start2 < smallest:
                            smallest = start2
                        elif stop2 > largest:
                            largest = stop2
                elif float(r2.START_POS) <= float(r1.START_POS) and float(r2.STOP_POS) >= float(r1.STOP_POS):
                    caller = r2.CALLER
                    start2 = float(r2.START_POS)
                    stop2 = float(r2.STOP_POS)
                    start1 = float(r1.START_POS)
                    stop1 = float(r1.STOP_POS)
                    sv_type = r2.SV_TYPE
                    overlap = (stop1 - start1) / (stop2 - start2) * 100
                    if overlap >= float(percent_threshold):
                        info = "ID=" + caller + ':' + "START=" + str(start2) + ':' + "END=" + str(stop2) + ':' + sv_type + ';' 
                        pos_overlap += info
                        if start2 < smallest:
                            smallest = start2
                        elif stop2 > largest:
                            largest = stop2

                elif float(r2.START_POS) >= float(r1.START_POS) and float(r2.STOP_POS) <= float(r1.STOP_POS):
                    caller = r2.CALLER
                    start2 = float(r2.START_POS)
                    stop2 = float(r2.STOP_POS)
                    start1 = float(r1.START_POS)
                    stop1 = float(r1.STOP_POS)
                    sv_type = r2.SV_TYPE
                    overlap = (stop2 - start2) / (stop1 - start1) * 100 

                    if overlap >= float(percent_threshold):
                        info = "ID=" + caller + ':' + "START=" + str(start2) + ':' + "END=" + str(stop2) + ':' + sv_type + ';' 
                        pos_overlap += info 
                        if start2 < smallest:
                            smallest = start2
                        elif stop2 > largest:
                            largest = stop2
                
        all_overlaps.append(pos_overlap)
        range_info = str(smallest) + '-' + str(largest)
        min_max.append(range_info)
    how_many = []
    for item in all_overlaps:
        each_type = item.split(';')
        del each_type[-1]
        all_types = []
        for thing in each_type:
            callerID = thing.split(':')[0]
            all_types.append(callerID)
        all_types = list(dict.fromkeys(all_types))
        total = len(all_types)  
        how_many.append(total)
    sample_SV['OVERLAPS'] = all_overlaps
    sample_SV['NUM_CALLERS'] = how_many
    sample_SV['RANGE'] = min_max
    dedup_SV = sample_SV.drop_duplicates('RANGE', keep='first')
    dedup_SV_2 = dedup_SV.drop_duplicates('RANGE', keep='first')
    high_num = dedup_SV_2.loc[dedup_SV_2['NUM_CALLERS'] >= caller_threshold]
    ranges = high_num['RANGE'].to_list()
    chroms = high_num['CHROM'].to_list()
    callers = high_num['OVERLAPS'].to_list()
    num_callers = high_num['NUM_CALLERS'].to_list()
    all_starts = []
    all_stops = []
    lengths = []
    for nums in ranges:
        start = nums.split('-')[0]
        stop = nums.split('-')[1]
        length = float(stop) - float(start)
        all_starts.append(start)
        all_stops.append(stop)
        lengths.append(length)    
    everything = {'CHROM': chroms, 'START': all_starts, 'STOP': all_stops, 'SV_CALLED': callers, 'NUM_CALLERS': num_callers, 'SV_LEN': lengths}
    new_frame = pd.DataFrame(everything)   
    return new_frame
  
def most_frequent(List): 
    return max(set(List), key = List.count)


def mergeTableStuff(r1,r2,sample_table):
    info1 = r1.SV_CALLED.split(';')
    info2 = r2.SV_CALLED.split(';')
    del info1[-1]
    del info2[-1]
    new_info = info1 + list(set(info2) - set(info1))
    info_str = ';'
    info_str = info_str.join(new_info) 
    info_str += ';'
    new_length = float(r1.STOP) - float(r2.START)
    sample_table.at[r1.Index, 'SV_CALLED'] = info_str
    sample_table.at[r1.Index, 'START'] = float(r2.START)
    sample_table.at[r1.Index, 'STOP'] = float(r1.STOP)
    sample_table.at[r1.Index, 'CHROM'] = r1.CHROM
    sample_table.at[r1.Index, 'SV_LEN'] = new_length
    sample_table.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
    sample_table.drop(r2.Index, inplace=True)
    return sample_table
def dedup_reformat(sample_table, over_t):

    over_t = float(over_t)
    for r1 in sample_table.itertuples():
        for r2 in sample_table.itertuples():
            if r1.Index != r2.Index:
                if r1.CHROM == r2.CHROM:
                    if float(r2.START) <= float(r1.START) and float(r2.STOP) < float(r1.STOP) and float(r2.STOP) > float(r1.START):
                        overlap = (float(r2.STOP) - float(r1.START)) / (float(r1.STOP) - float(r1.START)) * 100
                        overlap_2 = (float(r2.STOP) - float(r1.START)) / (float(r2.STOP) - float(r2.START)) * 100
                        if overlap >= over_t and overlap_2 >= over_t:
                            info1 = r1.SV_CALLED.split(';')
                            info2 = r2.SV_CALLED.split(';')
                            del info1[-1]
                            del info2[-1]
                            new_info = info1 + list(set(info2) - set(info1))
                            info_str = ';'
                            info_str = info_str.join(new_info) 
                            info_str += ';'
                            new_length = float(r1.STOP) - float(r2.START)
                            sample_table.at[r1.Index, 'SV_CALLED'] = info_str
                            sample_table.at[r1.Index, 'START'] = float(r2.START)
                            sample_table.at[r1.Index, 'STOP'] = float(r1.STOP)
                            sample_table.at[r1.Index, 'CHROM'] = r1.CHROM
                            sample_table.at[r1.Index, 'SV_LEN'] = new_length
                            sample_table.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                            sample_table.drop(r2.Index, inplace=True)
                            #sample_table = mergeTableStuff(r1,r2, sample_table)
                    elif float(r2.START) > float(r1.START) and float(r2.STOP) >= float(r1.STOP) and float(r2.START) < float(r1.STOP):
                        overlap = (float(r1.STOP) - float(r2.START)) / (float(r1.STOP) - float(r1.START)) * 100
                        overlap_2 = (float(r1.STOP) - float(r2.START)) / (float(r2.STOP) - float(r2.START)) * 100
                        if overlap >= over_t and overlap_2 >= over_t:
                            info1 = r1.SV_CALLED.split(';')
                            info2 = r2.SV_CALLED.split(';')
                            del info1[-1]
                            del info2[-1]
                            new_info = info1 + list(set(info2) - set(info1))
                            info_str = ';'
                            info_str = info_str.join(new_info) 
                            info_str += ';'
                            new_length = float(r2.STOP) - float(r1.START)
                            sample_table.at[r1.Index, 'SV_CALLED'] = info_str
                            sample_table.at[r1.Index, 'CHROM'] = r1.CHROM
                            sample_table.at[r1.Index, 'SV_LEN'] = new_length
                            sample_table.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                            sample_table.at[r1.Index, 'STOP'] = float(r2.STOP)
                            sample_table.at[r1.Index, 'START'] = float(r1.START)
                            sample_table.drop(r2.Index, inplace=True)
                    elif float(r2.START) <= float(r1.START) and float(r2.STOP) >= float(r1.STOP):
                        overlap = (float(r1.STOP) - float(r1.START)) / (float(r2.STOP) - float(r2.START)) * 100
                        if overlap >= over_t:
                            info1 = r1.SV_CALLED.split(';')
                            info2 = r2.SV_CALLED.split(';')
                            del info1[-1]
                            del info2[-1]
                            new_info = info1 + list(set(info2) - set(info1))
                            info_str = ';'
                            info_str = info_str.join(new_info)
                            info_str += ';'
                            new_length = float(r2.STOP) - float(r2.START)
                            sample_table.at[r1.Index, 'SV_CALLED'] = info_str
                            sample_table.at[r1.Index, 'START'] = float(r2.START)
                            sample_table.at[r1.Index, 'STOP'] = float(r2.STOP)
                            sample_table.at[r1.Index, 'CHROM'] = r1.CHROM
                            sample_table.at[r1.Index, 'SV_LEN'] = new_length
                            sample_table.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                            sample_table.drop(r2.Index, inplace=True)
                    elif float(r2.START) >= float(r1.START) and float(r2.STOP) <= float(r1.STOP):
                        overlap = (float(r2.STOP) - float(r2.START)) / (float(r1.STOP) - float(r1.START)) * 100
                        if overlap >= over_t:
                            info1 = r1.SV_CALLED.split(';')
                            info2 = r2.SV_CALLED.split(';')
                            del info1[-1]
                            del info2[-1]
                            new_info = info1 + list(set(info2) - set(info1))
                            info_str = ';'
                            info_str = info_str.join(new_info) 
                            info_str += ';'
                            sample_table.at[r1.Index, 'SV_CALLED'] = info_str
                            sample_table.at[r1.Index, 'START'] = float(r1.START)
                            sample_table.at[r1.Index, 'STOP'] = float(r1.STOP)
                            sample_table.at[r1.Index, 'CHROM'] = r1.CHROM
                            sample_table.at[r1.Index, 'SV_LEN'] = r1.SV_LEN
                            sample_table.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                            sample_table.drop(r2.Index, inplace=True)

    breakdancer = []
    pindel = []
    lumpy = []
    manta = []
    gridss = []
    delly = []

    for row1 in sample_table.itertuples():
        br_ss = []
        br_t = []
        pi_ss = []
        pi_t = []
        lu_ss = []
        lu_t = []
        ma_ss = []
        ma_t = []
        gr_ss = []
        gr_t = []
        de_ss = []
        de_t = []
        ind = row1.SV_CALLED.split(';')
        del ind[-1]
        for thing in ind:
            each = thing.split(':')
            ID = each[0].split('=')[1]
            START = each[1].split('=')[1]
            STOP = each[2].split('=')[1]
            TYPE = each[3]
            if ID == 'breakdancer':
                br_ss.append(float(START))
                br_ss.append(float(STOP))
                br_t.append(TYPE)
            elif ID == 'pindel':
                pi_ss.append(float(START))
                pi_ss.append(float(STOP))
                pi_t.append(TYPE)
            elif ID == 'lumpy':
                lu_ss.append(float(START))
                lu_ss.append(float(STOP))
                lu_t.append(TYPE)
            elif ID == 'delly':
                de_ss.append(float(START))
                de_ss.append(float(STOP))
                de_t.append(TYPE) 
            elif ID == 'manta':
                ma_ss.append(float(START))
                ma_ss.append(float(STOP))
                ma_t.append(TYPE)
            elif ID == 'gridss':
                gr_ss.append(float(START))
                gr_ss.append(float(STOP))
                gr_t.append(TYPE)      
        if len(br_ss) != 0:
            s = float(min(br_ss))
            e = float(max(br_ss))
            percent_o = (e - s) / (float(row1.STOP) - float(row1.START)) * 100
            ninfo = str(s) + '-' + str(e) + ':' + most_frequent(br_t) + ':' + str(round(percent_o, 2)) + '%'
            breakdancer.append(ninfo)
        if len(br_t) == 0:
            x = 'X'
            breakdancer.append(x)
        if len(pi_ss) != 0:
            s = float(min(pi_ss))
            e = float(max(pi_ss))
            percent_o = (e - s) / (float(row1.STOP) - float(row1.START)) * 100
            ninfo = str(s) + '-' + str(e) + ':' + most_frequent(pi_t) + ':' + str(round(percent_o, 2)) + '%'
            pindel.append(ninfo)
        if len(pi_t) == 0:
            x = 'X'
            pindel.append(x)
        if len(lu_ss) != 0:
            s = float(min(lu_ss))
            e = float(max(lu_ss))
            percent_o = (e - s) / (float(row1.STOP) - float(row1.START)) * 100
            ninfo = str(s) + '-' + str(e) + ':' + most_frequent(lu_t) + ':' + str(round(percent_o, 2)) + '%'
            lumpy.append(ninfo)
        if len(lu_t) == 0:
            x = 'X'
            lumpy.append(x)
        if len(de_ss) != 0:
            s = float(min(de_ss))
            e = float(max(de_ss))
            percent_o = (e - s) / (float(row1.STOP) - float(row1.START)) * 100
            ninfo = str(s) + '-' + str(e) + ':' + most_frequent(de_t) + ':' + str(round(percent_o, 2)) + '%'
            delly.append(ninfo)
        if len(de_t) == 0:
            x = 'X'
            delly.append(x)
        if len(ma_ss) != 0:
            s = float(min(ma_ss))
            e = float(max(ma_ss))
            percent_o = (e - s) / (float(row1.STOP) - float(row1.START)) * 100
            ninfo = str(s) + '-' + str(e) + ':' + most_frequent(ma_t) + ':' + str(round(percent_o, 2)) + '%'
            manta.append(ninfo)
        if len(ma_t) == 0:
            x = 'X'
            manta.append(x)
        if len(gr_ss) != 0:
            s = float(min(gr_ss))
            e = float(max(gr_ss))
            percent_o = (e - s) / (float(row1.STOP) - float(row1.START)) * 100
            ninfo = str(s) + '-' + str(e) + ':' + most_frequent(gr_t) + ':' + str(round(percent_o, 2)) + '%'
            gridss.append(ninfo)
        if len(gr_t) == 0:
            x = 'X'
            gridss.append(x)
    starts = sample_table['START'].to_list()
    stops = sample_table['STOP'].to_list()
    chroms = sample_table['CHROM'].to_list()
    num_callers = sample_table['NUM_CALLERS'].to_list()
    all_lengths = sample_table['SV_LEN'].to_list()
    all_reformat = {'CHROM': chroms, 'START': starts, 'STOP': stops,'NUM_CALLERS': num_callers, 'SV_LEN': all_lengths, 'BREAKDANCER': breakdancer, 'PINDEL': pindel, 'LUMPY': lumpy, 'MANTA': manta, 'DELLY': delly, 'GRIDSS': gridss}
    reformat_frame = pd.DataFrame(all_reformat)
    return reformat_frame



def main():
    data = sys.argv[1]
    out = sys.argv[2]



