import pandas as pd
import io
import os 
import sys 


def overlap_samples(all_samples, over_t):
    all_samples['SAMPLE_COUNTS'] = ''
    over_t = int(over_t)
    for r1 in all_samples.itertuples():
        all_callers = []
        for r2 in all_samples.itertuples():
            if r1.Index != r2.Index:
                if r1.SAMPLE != r2.SAMPLE:
                    if r1.CHROM == r2.CHROM:
                        if float(r2.START) <= float(r1.START) and float(r2.STOP) < float(r1.STOP) and float(r2.STOP) > float(r1.START):
                            overlap = (float(r2.STOP) - float(r1.START)) / (float(r1.STOP) - float(r1.START)) * 100
                            overlap_2 = (float(r2.STOP) - float(r1.START)) / (float(r2.STOP) - float(r2.START)) * 100
                            if overlap >= over_t and overlap_2 >= over_t:
                                if r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE) + ':' + str(r2.SAMPLE)
                                elif r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS != '':
                                    new_sample_info = str(r1.SAMPLE) + ':' + str(r2.SAMPLE_COUNTS)
                                elif r1.SAMPLE_COUNTS != '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE_COUNTS) + ':' + str(r2.SAMPLE)
                                else:
                                    info1 = r1.SAMPLE_COUNTS.split(':')
                                    info2 = r2.SAMPLE_COUNTS.split(':')
                                    new_info = info1 + list(set(info2) - set(info1))
                                    new_sample_info = ':'
                                    new_sample_info = new_sample_info.join(new_info)

                                new_b = r1.BREAKDANCER + ';' + r2.BREAKDANCER
                                new_p = r1.PINDEL + ';' + r2.PINDEL
                                new_l = r1.LUMPY + ';' + r2.LUMPY
                                new_m = r1.MANTA + ';' + r2.MANTA
                                new_d = r1.DELLY + ';' + r2.DELLY
                                new_g = r1.GRIDSS + ';' + r2.GRIDSS
                                new_length = float(r1.STOP) - float(r2.START)
                                all_samples.at[r1.Index, 'START'] = r2.START
                                all_samples.at[r1.Index, 'SAMPLE'] = r1.SAMPLE
                                all_samples.at[r1.Index, 'SAMPLE_COUNTS'] = new_sample_info
                                all_samples.at[r1.Index, 'STOP'] = r1.STOP
                                all_samples.at[r1.Index, 'CHROM'] = r1.CHROM
                                all_samples.at[r1.Index, 'SV_LEN'] = new_length
                                all_samples.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                                all_samples.at[r1.Index, 'BREAKDANCER'] = new_b
                                all_samples.at[r1.Index, 'PINDEL'] = new_p
                                all_samples.at[r1.Index, 'LUMPY'] = new_l
                                all_samples.at[r1.Index, 'MANTA'] = new_m
                                all_samples.at[r1.Index, 'DELLY'] = new_d
                                all_samples.at[r1.Index, 'GRIDSS'] = new_d
                                all_samples.drop(r2.Index, inplace=True)

                        elif float(r2.START) > float(r1.START) and float(r2.STOP) >= float(r1.STOP) and float(r2.START) < float(r1.STOP):
                            overlap = (float(r1.STOP) - float(r2.START)) / (float(r1.STOP) - float(r1.START)) * 100
                            overlap_2 = (float(r1.STOP) - float(r2.START)) / (float(r2.STOP) - float(r2.START)) * 100
                            if overlap >= over_t and overlap_2 >= over_t:
                                if r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE) + ':' + str(r2.SAMPLE)
                                elif r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS != '':
                                    new_sample_info = str(r1.SAMPLE) + ':' + str(r2.SAMPLE_COUNTS)
                                elif r1.SAMPLE_COUNTS != '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE_COUNTS) + ':' + str(r2.SAMPLE)
                                else:
                                    info1 = r1.SAMPLE_COUNTS.split(':')
                                    info2 = r2.SAMPLE_COUNTS.split(':')
                                    new_info = info1 + list(set(info2) - set(info1))
                                    new_sample_info = ':'
                                    new_sample_info = new_sample_info.join(new_info)
                                new_length = float(r2.STOP) - float(r1.START)
                                new_b = r1.BREAKDANCER + ';' + r2.BREAKDANCER
                                new_p = r1.PINDEL + ';' + r2.PINDEL
                                new_l = r1.LUMPY + ';' + r2.LUMPY
                                new_m = r1.MANTA + ';' + r2.MANTA
                                new_d = r1.DELLY + ';' + r2.DELLY
                                new_g = r1.GRIDSS + ';' + r2.GRIDSS
                                all_samples.at[r1.Index, 'SAMPLE_COUNTS'] = new_sample_info
                                all_samples.at[r1.Index, 'CHROM'] = r1.CHROM
                                all_samples.at[r1.Index, 'SV_LEN'] = new_length
                                all_samples.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                                all_samples.at[r1.Index, 'STOP'] = r2.STOP
                                all_samples.at[r1.Index, 'START'] = r1.START
                                all_samples.at[r1.Index, 'SAMPLE'] = r1.SAMPLE
                                all_samples.at[r1.Index, 'BREAKDANCER'] = new_b
                                all_samples.at[r1.Index, 'PINDEL'] = new_p
                                all_samples.at[r1.Index, 'LUMPY'] = new_l
                                all_samples.at[r1.Index, 'MANTA'] = new_m
                                all_samples.at[r1.Index, 'DELLY'] = new_d
                                all_samples.at[r1.Index, 'GRIDSS'] = new_d
                                all_samples.drop(r2.Index, inplace=True)

                        elif float(r2.START) <= float(r1.START) and float(r2.STOP) >= float(r1.STOP):

                            overlap = (float(r1.STOP) - float(r1.START)) / (float(r2.STOP) - float(r2.START)) * 100
                            if overlap >= over_t:
                                if r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE) + ':' +  str(r2.SAMPLE)
                                elif r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS != '':
                                    new_sample_info = str(r1.SAMPLE) + ':' + str(r2.SAMPLE_COUNTS)
                                elif r1.SAMPLE_COUNTS != '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE_COUNTS) + ':' +  str(r2.SAMPLE)
                                else:
                                    info1 = r1.SAMPLE_COUNTS.split(':')
                                    info2 = r2.SAMPLE_COUNTS.split(':')
                                    new_info = info1 + list(set(info2) - set(info1))
                                    new_sample_info = ':'
                                    new_sample_info = new_sample_info.join(new_info)
                                new_length = float(r2.STOP) - float(r2.START)
                                new_b = r1.BREAKDANCER + ';' + r2.BREAKDANCER
                                new_p = r1.PINDEL + ';' + r2.PINDEL
                                new_l = r1.LUMPY + ';' + r2.LUMPY
                                new_m = r1.MANTA + ';' + r2.MANTA
                                new_d = r1.DELLY + ';' + r2.DELLY
                                new_g = r1.GRIDSS + ';' + r2.GRIDSS
                                all_samples.at[r1.Index, 'SAMPLE_COUNTS'] = new_sample_info
                                all_samples.at[r1.Index, 'START'] = r2.START
                                all_samples.at[r1.Index, 'STOP'] = r2.STOP
                                all_samples.at[r1.Index, 'CHROM'] = r1.CHROM
                                all_samples.at[r1.Index, 'SV_LEN'] = new_length
                                all_samples.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                                all_samples.at[r1.Index, 'SAMPLE'] = r1.SAMPLE
                                all_samples.at[r1.Index, 'BREAKDANCER'] = new_b
                                all_samples.at[r1.Index, 'PINDEL'] = new_p
                                all_samples.at[r1.Index, 'LUMPY'] = new_l
                                all_samples.at[r1.Index, 'MANTA'] = new_m
                                all_samples.at[r1.Index, 'DELLY'] = new_d
                                all_samples.at[r1.Index, 'GRIDSS'] = new_d
                                all_samples.drop(r2.Index, inplace=True)


                        elif float(r2.START) >= float(r1.START) and float(r2.STOP) <= float(r1.STOP):

                            overlap = (float(r2.STOP) - float(r2.START)) / (float(r1.STOP) - float(r1.START)) * 100
                            if overlap >= over_t:
                                if r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE) + ':' +  str(r2.SAMPLE)
                                elif r1.SAMPLE_COUNTS == '' and r2.SAMPLE_COUNTS != '':
                                    new_sample_info = str(r1.SAMPLE) + ':' + str(r2.SAMPLE_COUNTS)
                                elif r1.SAMPLE_COUNTS != '' and r2.SAMPLE_COUNTS == '':
                                    new_sample_info = str(r1.SAMPLE_COUNTS) + ':' + str(r2.SAMPLE)
                                else:
                                    info1 = r1.SAMPLE_COUNTS.split(':')
                                    info2 = r2.SAMPLE_COUNTS.split(':')
                                    new_info = info1 + list(set(info2) - set(info1))
                                    new_sample_info = ':'
                                    new_sample_info = new_sample_info.join(new_info)
                                new_b = r1.BREAKDANCER + ';' + r2.BREAKDANCER
                                new_p = r1.PINDEL + ';' + r2.PINDEL
                                new_l = r1.LUMPY + ';' + r2.LUMPY
                                new_m = r1.MANTA + ';' + r2.MANTA
                                new_d = r1.DELLY + ';' + r2.DELLY
                                new_g = r1.GRIDSS + ';' + r2.GRIDSS
                                all_samples.at[r1.Index, 'SAMPLE_COUNTS'] = new_sample_info
                                all_samples.at[r1.Index, 'START'] = r1.START
                                all_samples.at[r1.Index, 'STOP'] = r1.STOP
                                all_samples.at[r1.Index, 'CHROM'] = r1.CHROM
                                all_samples.at[r1.Index, 'SV_LEN'] = r1.SV_LEN
                                all_samples.at[r1.Index, 'NUM_CALLERS'] = r1.NUM_CALLERS
                                all_samples.at[r1.Index, 'SAMPLE'] = r1.SAMPLE
                                all_samples.at[r1.Index, 'BREAKDANCER'] = new_b
                                all_samples.at[r1.Index, 'PINDEL'] = new_p
                                all_samples.at[r1.Index, 'LUMPY'] = new_l
                                all_samples.at[r1.Index, 'MANTA'] = new_m
                                all_samples.at[r1.Index, 'DELLY'] = new_d
                                all_samples.at[r1.Index, 'GRIDSS'] = new_d
                                all_samples.drop(r2.Index, inplace=True)

    all_samples
    how_many = []
    dedup_samples = []
    for row in all_samples.itertuples():
        if row.SAMPLE_COUNTS == '':
            how_many.append(1)
            dedup_samples.append(str(row.SAMPLE))
        else:
            called_sams = row.SAMPLE_COUNTS.split(':')
            called_sams = list(dict.fromkeys(called_sams))
            total = len(called_sams)
            how_many.append(total)
            the_samples = ':'
            the_samples = the_samples.join(called_sams)
            dedup_samples.append(the_samples)


    all_samples['SAMPLE_COUNTS'] = dedup_samples
    all_samples['NUM_SAMPLES'] = how_many
    all_samples.sort_values(by='CHROM', inplace=True)
    return all_samples


def main():
    mergeDir = sys.argv[1]
    print(mergeDir)
    allFrames = []
    for file in os.listdir(mergeDir):
        i = 0
        if file.endswith("_merged.txt"):
            frame = pd.read_csv(file, sep='\t')
            allFrames.append(frame)

    all_samples = pd.concat(allFrames, axis=0, ignore_index=True)
    all_samples = overlap_samples(all_samples, 75)
    all_samples.sort_values(by=['CHROM','START'], inplace=True)

    remove = all_samples.loc[all_samples['NUM_SAMPLES'] < int(sys.argv[2])]
    remove.to_csv(mergeDir + "/noShareSV.txt", sep="\t", index=False)
    bed_format = remove.filter(['CHROM','START','STOP'], axis=1)
    bed_format.to_csv(mergeDir + "/noShareSV.bed", sep='\t', index=False, header=False)         




main()