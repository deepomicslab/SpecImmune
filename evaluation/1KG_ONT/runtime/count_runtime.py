import os,re
import pandas as pd

def extract_time(time_file): #log file obtained by /usr/bin/time -v
        #if no time available
    if os.path.isfile(time_file):
        for line in open(time_file):
            time_re = re.search('User time \(seconds\):(.*?)$', line)
            if time_re:
                user_time =  time_re.group(1).strip()

            time_sys = re.search('System time \(seconds\):(.*?)$', line)
            if time_sys:
                sys_time = time_sys.group(1).strip()
        # print (user_time, sys_time)
        all_time = float(user_time) + float(sys_time)
        final_time = round(all_time/3600, 2)
        return final_time
    else:
        return None
    

def extract_wall_clock_time(time_file): #log file obtained by /usr/bin/time -v
        #if no time available
    for line in open(time_file):
        time_re = re.search('Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):(.*?)$', line)
        if time_re:
            wall_block_time = time_re.group(1).strip()
    array = wall_block_time.split(":")
    hours = int(array[0].strip()) + float(array[1].strip())/60
    return hours

def extract_mem(time_file):
    used_mem = 0 #if no time available
    for line in open(time_file):
        time_re = re.search('Maximum resident set size \(kbytes\):(.*?)$', line)
        if time_re:
            used_mem =  time_re.group(1).strip()
    final_mem = round(float(used_mem)/1000000, 2)
    return final_mem

def for_each_file(time_file):
    time = extract_time(time_file)
    wall_clock_time = extract_wall_clock_time(time_file)
    mem = extract_mem(time_file)
    return time, wall_clock_time, mem

def batch_count(folder, output):
    data = []
    # folder = "/home/wangshuai/00.hla/long/experiments/vdj/vdj_results/hprc_hifi/"
    for subfolder in os.listdir(folder):
        if os.path.isdir(folder+subfolder):
            time_file = folder+f"/{subfolder}.time"
            time, wall_clock_time, mem = for_each_file(time_file)
            print(subfolder, time, wall_clock_time, mem)
            data.append([subfolder, time, wall_clock_time, mem])
    df = pd.DataFrame(data, columns = ['sample', 'time', 'wall_clock_time', 'mem'])
    df.to_csv(output, index=False)


if __name__ == '__main__':
    # time_file = '/mnt/d/HLAPro_backup/Nanopore_optimize/test.time'
    # time, wall_clock_time, mem = for_each_file(time_file)
    # print(time, wall_clock_time, mem)

    folder = "/home/wangshuai/00.hla/long/experiments/vdj/vdj_results/hprc_hifi/"
    output = "vdj_time.csv"
    batch_count(folder, output)

    folder = "/home/wangshuai/00.hla/long/experiments/kir/kir_results/hprc_hifi_time/"
    output = "kir_time.csv"
    batch_count(folder, output)

    folder = "/disk2/workspace/wangshuai/03.hla/cyp_results/spec_hprc_time/"
    output = "cyp_time.csv"
    batch_count(folder, output)


