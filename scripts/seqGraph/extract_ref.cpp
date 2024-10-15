#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <sys/stat.h>
#include <pthread.h>
#include <thread>
#include <sstream>
#include <map>
#include <mutex>


using namespace std;

const short c = 300;
const char coder_num = 3;
const unsigned char least_depth = 3;
const int k = 32;
long array_size = pow(2, k);
char *kmer_count_table = new char[array_size];

int MIN_BASE_NUM = 6; // 6
int REF_NEAR = 500; // 500
int DIFF = 2; // 2
int PEAK_W = 5; // 5
int NEAR = 0; // 0
int SKIP_N = 2*k; // 2*k
int SKIP_S = k; //k
int SKIP_A = 1;
int MIN_READS = 1; // 1
int MAX_PEAK_NUM = 300000000;
int thread_num;
std::mutex mtx;

class Split_reads{
public:
    map<int, int> chr_kmer_count;
    map<int, int> filter_kmer_count;
    map<int, int> chr_peak_index;
    string read_name;
    int min_base_num = MIN_BASE_NUM;
    int base_hits = 0;
    int base_kmer[3];
    int base_chr[3];
    int tem_index = 0;
    bool valid_base = false;
    void count_peak_kmer(int peak_chr, int peak_pos, int peak_index, int coder_index);
    void check_split(unsigned char* peak_filter);
    void init_array(void);
    void judge_base(void);
};

void Split_reads::init_array(void){
    memset(base_kmer, 0, sizeof(int)*3);
}

void Split_reads::count_peak_kmer(int peak_chr, int peak_pos, int peak_index, int coder_index){
    base_kmer[coder_index] = peak_index;
    base_chr[coder_index] = peak_chr;
}

void Split_reads::judge_base(void){
    map<int, int>::iterator chr_iter;
    int select_chr = 0;
    int select_index = 0;
    int select_num = 0;
    bool flag = false;

    for (int i = 0; i < 3; i++){
        // //cout << base_kmer[i] << "\t" << base_chr[i] << endl;
        if (base_kmer[i] != 0){
            flag = true;
            chr_iter = chr_kmer_count.find(base_chr[i]);
            if (chr_iter != chr_kmer_count.end()){
                if (chr_iter->second >= select_num){
                    select_index = base_kmer[i];
                    select_chr = base_chr[i];
                    select_num = chr_iter->second;
                    // //cout << i << "##\t" << base_kmer[i]<<"\t" << base_chr[i] << "\t" << chr_iter->second << endl;
                }
            }
            else{
                // //cout << i << "##\t" << base_kmer[i]<<"\t" << base_chr[i] << "\t" << endl;
                if (select_index == 0){
                    select_index = base_kmer[i];
                    select_chr = base_chr[i];
                    select_num = 0;
                }
            }
        }
    }

    if (flag){
        if (chr_kmer_count.find(select_chr) == chr_kmer_count.end()){
            chr_kmer_count[select_chr] = 1;
            chr_peak_index[select_chr] = select_index; // the first peak of the genome
        }
        else{
            chr_kmer_count[select_chr] += 1;
        }
        base_hits += 1;
    }
}

void Split_reads::check_split(unsigned char* peak_filter){
    map<int, int>::iterator iter;
    int chr, peak_index;

    int largest = 0; // top2
    int second = 0;

    iter = chr_kmer_count.begin();
    while(iter != chr_kmer_count.end()){
        // //cout << iter->first << ":"<< iter->second <<endl;
        if (iter->second >= min_base_num){
            filter_kmer_count[iter->first] = iter->second;
            if (iter->second >= largest){
                second = largest;
                largest = iter->second;
            }
            else{
                if (iter->second >= second){
                    second = iter->second;
                }
            }
        }
        iter ++ ;
    }

    if (filter_kmer_count.size() > 1){
        // map<int, int>::iterator iter;
        iter = filter_kmer_count.begin();
        while(iter != filter_kmer_count.end()){
            // //cout << iter->first << ":"<< iter->second <<endl;
            if (iter->second == largest || iter->second == second){
                chr = iter->first;
                peak_index = chr_peak_index[chr];
                if ((int) peak_filter[peak_index] < 254){
                    peak_filter[peak_index] += 1;
                }
            }
            iter ++ ;
        }
        // //cout << read_name << "------------------"<<endl;
    }
}

class Peaks{
public:
    // int my_peak_index = 0;
    int *peak_index_array;
    int filter_peak_num = 0;
    int near = NEAR;
    int merge_close_peak = 50;
    int ref_gap = 500;
    int ref_near = REF_NEAR;
    int max_peak_num = MAX_PEAK_NUM;
    int each_thread_peak_num;
    int *peak_loci = new int[max_peak_num*2];
    unsigned char *peak_filter = new unsigned char[max_peak_num];
    unsigned int *peak_kmer = new unsigned int[array_size];
    int thread_num;
    void init(int thread_num);
    void add_peak(int ref_index, int pos, unsigned int* record_ref_index,int ref_len,
                  long & total_peak_num, unsigned char* record_ref_hit, int thread_index);
    void delete_array(void);
    void slide_reads(string fastq_file, string fastq_file_2, char* coder, int* base,
                     char* comple, short *choose_coder, int down_sam_ratio, long start, long end);
    void count_filtered_peak(string interval_name);
    bool merge_peak(int ref_index, int pos, int my_peak_index);
};

void Peaks::init(int obtain_thread_num){
    thread_num = obtain_thread_num;
    peak_index_array = new int[thread_num];
    each_thread_peak_num = max_peak_num/thread_num;
    for (int j = 0; j < thread_num; j++){
        peak_index_array[j] = 0 + each_thread_peak_num*j;
    }

}

void Peaks::add_peak(int ref_index, int pos, unsigned int* record_ref_index, int ref_len,
                     long & total_peak_num, unsigned char* record_ref_hit, int thread_index){
    // //cout << Peaks::merge_peak(ref_index, pos)<<endl;
    int my_peak_index = peak_index_array[thread_index];
    if (Peaks::merge_peak(ref_index, pos,my_peak_index)){
        int index;
        for (int near_pos = pos - near; near_pos < pos + 1; near_pos++){
            if (near_pos>=0 & near_pos<=ref_len){
                for (int p = 0; p < 3; p++){
                    index = coder_num*near_pos+p;
                    if (record_ref_hit[index] > 0){
                        peak_kmer[record_ref_index[index]] = my_peak_index - 1;
                    }
                }
            }
        }
    }
    else{
        peak_loci[2*my_peak_index] = ref_index;
        peak_loci[2*my_peak_index+1] = pos;
        int index;
        for (int near_pos = pos - near; near_pos < pos + 1; near_pos++){
            if (near_pos>=0 & near_pos<=ref_len){
                for (int p = 0; p < 3; p++){
                    index = coder_num*near_pos+p;
                    if (record_ref_hit[index] > 0){
                        peak_kmer[record_ref_index[index]] = my_peak_index;
                    }
                }
            }
        }
        // if (peak_index_array[thread_index] - each_thread_peak_num * thread_index >= each_thread_peak_num){
        if (my_peak_index >= max_peak_num){
            //cout <<"Too many peaks! We recommand reduce the sampling size, or you can appoint larger max_peak_num."<<endl;
        }
        my_peak_index += 1;
        total_peak_num += 1;
    }
    peak_index_array[thread_index] = my_peak_index;

}

bool Peaks::merge_peak(int ref_index, int pos, int my_peak_index){
    bool exist_peak = false;
    if (my_peak_index > 0){
        // if (ref_index == peak_loci[2*my_peak_index-2] & pos - peak_loci[2*my_peak_index-1] < merge_close_peak){
        //     exist_peak = true;
        // }
        // with density peaks, the position of all the loci is the first one, the error can be heavy.
        // so use window~50bp, the peaks locate in the same window are merged as one, which the pos is that of the first adding peak's.
        if (ref_index == peak_loci[2*my_peak_index-2] & pos/merge_close_peak == peak_loci[2*my_peak_index-1]/merge_close_peak){
            exist_peak = true;
        }
    }
    return exist_peak;
}

string get_read_ID(string reads_seq){
    string delimiter = "/";
    string read_name_forward = reads_seq.substr(0, reads_seq.find(delimiter));
    delimiter = " ";
    read_name_forward = read_name_forward.substr(0, read_name_forward.find(delimiter));
    delimiter = "\t";
    read_name_forward = read_name_forward.substr(0, read_name_forward.find(delimiter));
    return read_name_forward;
}

void Peaks::slide_reads(string fastq_file, string fastq_file_2, char* coder, int* base,
                        char* comple, short *choose_coder, int down_sam_ratio, long start, long end){

    time_t t0 = time(0);
    ifstream fq_file;
    ifstream fq_file_2;
    fq_file.open(fastq_file);
    fq_file_2.open(fastq_file_2);
    string reads_seq, reads_seq_2;
    int reads_int [300];
    int reads_comple_int [300];
    unsigned int lines = 0;
    int converted_reads [900];
    int complemented_reads [900];
    int m;
    int n;
    unsigned int kmer_index, comple_kmer_index, real_index, b;
    int r ;
    // short read_len = 0;
    short read_len_1, read_len_2;
    int chr_index, peak_locus,peak_index;
    long pos = 0;
    long second_pos;


    string read_name;

    for (long i = start; i>0; i--){
        fq_file.seekg(i, ios::beg);
        char j;
        fq_file.get(j);
        if (j == '@'){ //not only read name has this symbol.
            if (i == 1){
                pos = i;
                break;
            }
            else{
                fq_file.seekg(i-1, ios::beg);
                char j;
                fq_file.get(j);
                if (j == '\n'){ //the @ should be the first char in a line
                    fq_file.seekg(i, ios::beg);
                    getline(fq_file,reads_seq);
                    string read_name_forward = get_read_ID(reads_seq);
                    if (read_name_forward.length() < 40){
                        pos = i;
                        break;
                    }
                }
            }
        }
    }
    fq_file.seekg(pos, ios::beg);
    fq_file_2.seekg(pos, ios::beg);
    long add_size = start;
    // while (fq_file >> reads_seq)
    while (getline(fq_file,reads_seq))
    {
        // fq_file_2>>reads_seq_2;
        getline(fq_file_2, reads_seq_2);
        if (start == 0 & lines % 1000000 == 1000000-1){
            //cout <<start<< " recheck reads\t"<<lines<<endl;
        }
        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length();
        if (lines == 0){
            string read_name_forward = get_read_ID(reads_seq);
            string read_name_reverse = get_read_ID(reads_seq_2);
            // //cout << read_name_forward << "\t" << read_name_reverse << endl;

            if (read_name_forward != read_name_reverse){
                // throw "Paired-end reads not consistent!";
                second_pos = pos - 1000000000;
                if (second_pos < 1){
                    second_pos = 1;
                }
                fq_file_2.seekg(second_pos, ios::beg);
                // //cout << read_name_forward<<"\tshould be the same with\t"<< read_name_reverse<< endl;
            }
            // else{
            //     //cout << read_name_forward<<"\tis same with\t"<< read_name_reverse<< endl;
            // }
            int while_times = 0;
            while (read_name_forward != read_name_reverse){

                getline(fq_file_2, reads_seq_2);
                read_name_reverse = get_read_ID(reads_seq_2);
                while_times += 1;
                if (while_times > 1000000000){
                    //cout << "Too many iterations to make paired-end reads be consistent!\
                     Please use single thread for this sample to avoid the problem." << endl;
                    //cout << read_name_forward<<"\t<=>\t"<< read_name_reverse<< endl;
                    break;
                }
            }

            //cout << read_name_forward<<"\t<=>\t"<< read_name_reverse<< endl;


        }
        if (lines % 4 == 0){
            read_name = get_read_ID(reads_seq);
        }
        if (lines % 4 == 1){
            time_t t1 = time(0);
            // if (lines == 1){
            // if (reads_seq.length() <= reads_seq_2.length()){
            //     read_len = reads_seq.length();//cal read length
            // }
            // else{
            //     read_len = reads_seq_2.length();//cal read length
            // }
            read_len_1 = reads_seq.length();
            read_len_2 = reads_seq_2.length();

            // }
            r = rand() % 100;
            if (r < down_sam_ratio){
                Split_reads each_read;
                each_read.read_name = read_name;

                //forward reads
                for (int j = 0; j < read_len_1; j++){
                    reads_int[j] = (int)reads_seq[j];
                    reads_comple_int[j] = comple[reads_int[j]];
                }
                for (int j = 0; j < read_len_1-k+1; j++){
                    each_read.init_array();
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                            n = coder[c*choose_coder[(k-1-z)*3+i]+reads_comple_int[j+z]];
                            if (m == 5){
                                all_valid = false;
                                break;
                            }
                            kmer_index += m*base[z];
                            comple_kmer_index += n*base[(k-1-z)];

                        }
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }
                        else{
                            real_index = kmer_index;
                        }
                        // //cout << j << "\t" << i << "\t" << peak_kmer[real_index] << "\t" << peak_loci[2*peak_kmer[real_index]]<<endl;
                        if (all_valid & peak_kmer[real_index] != 0){
                            each_read.count_peak_kmer(peak_loci[2*peak_kmer[real_index]],
                                                      peak_loci[2*peak_kmer[real_index]+1], peak_kmer[real_index], i);
                        }
                    }
                    each_read.judge_base();
                }
                //reverse reads
                for (int j = 0; j < read_len_2; j++){
                    reads_int[j] = (int)reads_seq_2[j];
                    reads_comple_int[j] = comple[reads_int[j]];
                }
                for (int j = 0; j < read_len_2-k+1; j++){
                    each_read.init_array();
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                            n = coder[c*choose_coder[(k-1-z)*3+i]+reads_comple_int[j+z]];
                            if (m == 5){
                                all_valid = false;
                                break;
                            }
                            kmer_index += m*base[z];
                            comple_kmer_index += n*base[(k-1-z)];
                        }
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }
                        else{
                            real_index = kmer_index;
                        }
                        // //cout << j << "\t" << i << "\t" << peak_kmer[real_index] << "\t" << peak_loci[2*peak_kmer[real_index]]<<endl;
                        if (all_valid & peak_kmer[real_index] != 0){
                            each_read.count_peak_kmer(peak_loci[2*peak_kmer[real_index]],
                                                      peak_loci[2*peak_kmer[real_index]+1], peak_kmer[real_index], i);
                        }
                    }
                    each_read.judge_base();
                }
                if (each_read.base_hits >= each_read.min_base_num){
                    each_read.check_split(peak_filter);
                }

            }
        }
        lines ++;
    }
    fq_file.close();
    fq_file_2.close();
}

void Peaks::delete_array(void){
    delete [] peak_loci;
    delete [] peak_kmer;
    delete [] peak_filter;
    delete [] peak_index_array;
}

void Peaks::count_filtered_peak(string interval_name){
    ofstream interval_file;
    interval_file.open(interval_name, ios::out | ios::trunc);

    long final_ref_len = 0;
    for (int j = 0; j < thread_num; j++){
        int my_peak_index = peak_index_array[j];
        int start = 1;
        int end = 1;
        int chr = 1;
        for (int i = each_thread_peak_num*j; i < my_peak_index; i++){
            if ((int)peak_filter[i] >= MIN_READS ){
                filter_peak_num += 1;
                if (chr == peak_loci[2*i] & peak_loci[2*i+1]-ref_near - end < ref_gap){
                    end = peak_loci[2*i+1]+ref_near;
                }
                else{
                    interval_file << chr << "\t" << start << "\t"<< end << endl;
                    final_ref_len += (end - start);
                    chr = peak_loci[2*i];
                    start = peak_loci[2*i+1]-ref_near;
                    end = peak_loci[2*i+1]+ref_near;
                }

            }

        }
        interval_file << chr << "\t" << start << "\t"<< end << endl;
        final_ref_len += (end - start);
    }

    interval_file.close();
    //cout <<"final_ref_len is:\t"<<final_ref_len<<endl;
}

void slide_window(unsigned char* record_ref_hit, int ref_len, int ref_index, long & extract_ref_len,
                  ofstream & interval_file, float hit_ratio, float perfect_hit_ratio,
                  long & total_peak_num,Peaks & MyPeak,unsigned int* record_ref_index, int thread_index){
    //find density hits regions

    int frag_index = 0;
    int window = 500;
    int one_coder_bases = 0; // the number of bases supported by at least one coder in a window
    int three_coder_bases = 0; // the number of bases supported by three coder in a window
    int one_coder_min = window * hit_ratio;
    int three_coder_min = window * perfect_hit_ratio;
    bool good_window;
    bool conti_flag = false;
    int start = 0 ;
    int end = 0;
    int* save_good_intervals = new int[2*ref_len/window];
    short *single_hit_num = new short[ref_len];
    short *trio_hit_num = new short[ref_len];
    bool *peak_hit = new bool[ref_len];
    int* save_peak_intervals = new int[2*ref_len];
    int peak_index = 0;
    int left, right, diff;

    for (int j = 0; j < ref_len; j++){
        short hit_coder_num = 0;
        for (int p = 0; p < 3; p++){
            // if (record_ref_hit[ref_len*p+j] == true){
            if (record_ref_hit[coder_num*j+p] == least_depth){
                hit_coder_num += 1;
            }
        }
        if (hit_coder_num == 3){
            trio_hit_num[j] = 1;
        }
        else{
            trio_hit_num[j] = 0;
        }
        if (hit_coder_num > 0){
            single_hit_num[j] = 1;
        }
        else{
            single_hit_num[j] = 0;
        }

        if (j < window){
            if (hit_coder_num >= 1){
                one_coder_bases += 1;
            }
            if (hit_coder_num == 3){
                three_coder_bases += 1;
            }
        }
        else{
            one_coder_bases = one_coder_bases - single_hit_num[j-window] + single_hit_num[j];
            three_coder_bases = three_coder_bases - trio_hit_num[j-window] + trio_hit_num[j];
        }

        if (one_coder_bases >= one_coder_min & three_coder_bases >= three_coder_min){
            good_window = true;
        }
        else{
            good_window = false;
        }

        if (conti_flag == false & good_window == true){
            start = j - 2 * window;
            if (start < 1){
                start = 1;
            }
            conti_flag = true;
        }
        if (conti_flag == true & good_window == false){
            end = j + 2 * window;
            if (end > ref_len){
                end = ref_len;
            }
            if (frag_index > 0 & start - save_good_intervals[2*frag_index-1] < window ){
                save_good_intervals[2*frag_index-1] = end;
            }
            else{
                save_good_intervals[2 * frag_index] = start;
                save_good_intervals[2 * frag_index+1] = end;
                frag_index += 1;
            }
            conti_flag = false;
        }
    }
    // if (ref_index == 5){
    //     for (int j = 0; j < ref_len; j++){
    //         if(j > 133808 - 100 & j < 133808 + 100){
    //             //cout << j<<"\t"<<single_hit_num[j] << "\t"<<(int)single_hit_num[j]<<"\t"<<peak_hit[j]<<endl;
    //         }
    //     }
    // }

    if (conti_flag == true & good_window == true){
        end = ref_len;
        if (frag_index > 0 & start - save_good_intervals[2*frag_index-1] < window ){
            save_good_intervals[2*frag_index-1] = end;
        }
        else{
            save_good_intervals[2 * frag_index] = start;
            save_good_intervals[2 * frag_index+1] = end;
            frag_index += 1;
        }
    }

    int el = 0;
    for (int i = 0; i < frag_index; i++){
        el += (save_good_intervals[2*i+1] - save_good_intervals[2*i]);
    }
    auto el_ratio = float(el)/float(ref_len);
    if(el>0 && el_ratio>0.7)
        cout<<"ref_index\t"<<ref_index<<"\t"<<frag_index<<"\t"<<el<<"\t"<<ref_len<<"\t"<<el_ratio<<"\n";
    delete [] save_good_intervals;
    delete [] single_hit_num;
    delete [] trio_hit_num;
    delete [] save_peak_intervals;
    delete [] peak_hit;
    // return extract_ref_len;
}

float cal_tab_empty_rate(){
    double empty_num = 0;
    float empty_rate;
    float weak_rate = 0;
    double weak_num = 0;
    for (unsigned int j = 0; j < array_size; j++){
        if ((int)kmer_count_table[j] == 0){
            empty_num += 1;
        }
        // else{
        //     //cout << (int)kmer_count_table[j] <<empty_num<< endl;
        // }
        if ((int)kmer_count_table[j] != least_depth){
            weak_num += 1;
        }
        // //cout << j << (int)kmer_count_table[j]<< empty_num <<endl;


    }
    empty_rate = empty_num/array_size;
    weak_rate = weak_num/array_size;
    //cout << array_size << "\t" << weak_num<< "\t" << empty_num<<endl;
    //cout << array_size << "\t" << weak_rate<< "\t" << empty_rate<<endl;
    return empty_rate;
}

void read_ref(string fasta_file, char* coder, int* base, int k, char* comple,
              string index_name, short *choose_coder)
{
    ifstream fa_file;
    ofstream index_file;
    ofstream len_file;
    index_file.open(index_name, ios::out | ios::binary);
    len_file.open(fasta_file+".genome.len.txt", ios::out);
    fa_file.open(fasta_file, ios::in);
    string ref_seq, line_seq;
    ref_seq = "\0";
    int ref_len;
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index, real_index;
    int ref_index = 0;
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    string chr_name ;
    string pre_name = "start";
    time_t t0 = time(0);
    int covert_num, comple_num;
    short convert_ref[300];
    short complemented_ref[300];
    char support_coder;

    // save random coder
    for (int j = 0; j < 100; j++){
        index_file.write((char *)(&choose_coder[j]), sizeof(unsigned int));
    }
    //cout <<"Start index ref..."<<endl;

    // while (fa_file >> line_seq){ //remember the last chr
    while (getline(fa_file,line_seq)){
        if (line_seq[0] == '>'){
            chr_name = pre_name;
            pre_name = get_read_ID(line_seq).substr(1);
            // pre_name = line_seq.substr(1);
            // if (ref_index > 10){
            //     break;
            // }
            ref_len= ref_seq.length();
            slide_ref_len += ref_len;
            // //cout << chr_name <<"\t"<< ref_index <<"\t" << ref_len << endl;
            if (ref_len > k){   // just to skip the first empty seq
                len_file << chr_name <<"\t"<< ref_index <<"\t" << ref_len <<"\t" << slide_ref_len << endl;
                // /*
                for (int i = 0; i < 3; i++){
                    kmer_index = 0; // start kmer
                    comple_kmer_index = 0;
                }
                int *ref_int = new int[ref_len];
                int *ref_comple_int = new int[ref_len];
                for (int j = 0; j < ref_len; j++){
                    ref_int[j] = (int)ref_seq[j];
                    ref_comple_int[j] = comple[ref_int[j]];
                }
                index_file.write((char *)(&ref_len), sizeof(unsigned int));
                for (int j = 0; j < ref_len-k+1; j++){
                    support_coder = 0;
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+ref_int[j+z]];
                            if (m == 5){
                                all_valid = false;
                                break;
                            }
                            kmer_index += m*base[z];
                            comple_kmer_index += coder[c*choose_coder[(k-1-z)*3+i]+ref_comple_int[j+z]]*base[(k-1-z)];
                        }
                        // //cout <<kmer_index<<"\t"<< j<<"\t"<<comple_kmer_index<<"\t"<<reads_seq[j] <<endl;
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }
                        else{
                            real_index = kmer_index;
                        }
                        if (all_valid == false){
                            real_index = 0;
                        }
                        index_file.write((char *)(&real_index), sizeof(real_index));
                    }
                }
                delete [] ref_int;
                delete [] ref_comple_int;
                // */
            }
            // else{
            //     //cout << chr_name<<"\t is an incorrect genome, too short."<<endl;
            // }
            if (ref_index % 1000 == 0){
                time_t t1 = time(0);
                //cout << chr_name<<"\t" << ref_index << "\t" <<ref_len <<" bp\t"<<slide_ref_len<<" bp\t" <<t1-t0<< endl;
            }
            ref_index += 1;
            ref_seq = "\0";
        }
        else{
            ref_seq += line_seq;
        }
    }
    // for the last chr
    chr_name = pre_name;
    ref_len= ref_seq.length();
    slide_ref_len += ref_len;
    if (ref_len > k){
        for (int i = 0; i < 3; i++){
            kmer_index = 0; // start kmer
            comple_kmer_index = 0;
        }
        int *ref_int = new int[ref_len];
        int *ref_comple_int = new int[ref_len];
        for (int j = 0; j < ref_len; j++){
            ref_int[j] = (int)ref_seq[j];
            ref_comple_int[j] = comple[ref_int[j]];
        }
        index_file.write((char *)(&ref_len), sizeof(unsigned int));
        for (int j = 0; j < ref_len-k+1; j++){
            support_coder = 0;
            for (int i = 0; i < 3; i++){
                kmer_index = 0;
                comple_kmer_index = 0;
                bool all_valid = true;
                for (int z = 0; z<k; z++){
                    m = coder[c*choose_coder[z*3+i]+ref_int[j+z]];
                    if (m == 5){
                        all_valid = false;
                        break;
                    }
                    kmer_index += m*base[z];
                    comple_kmer_index += coder[c*choose_coder[(k-1-z)*3+i]+ref_comple_int[j+z]]*base[(k-1-z)];
                }
                if (kmer_index > comple_kmer_index){ //use a smaller index
                    real_index = comple_kmer_index;
                }
                else{
                    real_index = kmer_index;
                }
                if (all_valid == false){
                    real_index = 0;
                }
                index_file.write((char *)(&real_index), sizeof(real_index));
            }
        }
        delete [] ref_int;
        delete [] ref_comple_int;
        time_t t1 = time(0);
        len_file << chr_name <<"\t"<< ref_index <<"\t" << ref_len <<"\t" << slide_ref_len << endl;
        //cout << "Last chr:\t" << chr_name<<"\t" << ref_index << "\t" <<ref_len <<" bp\t"<<slide_ref_len<<" bp\t" <<t1-t0<< endl;
    }

    //cout << "Index is done."<< "\t" << extract_ref_len << endl;
    fa_file.close();
    index_file.close();
    len_file.close();
}

void read_index(char* coder, int* base, int k, char* comple, string index_name, string interval_name,
                short *choose_coder, float hit_ratio, float perfect_hit_ratio, Peaks & MyPeak, long start,
                long end, int start_ref_index, long&  extract_ref_len, long & slide_ref_len, long & total_peak_num, int thread_index){

    // long extract_ref_len = 0;
    // long slide_ref_len = 0;
    // long total_peak_num = 0;

    unsigned int *record_ref_index;
    ifstream index_file;
    ofstream interval_file;
    index_file.open(index_name, ios::in | ios::binary);
    interval_file.open(interval_name, ios::out | ios::trunc);
    unsigned int real_index;
    unsigned int kmer_index;
    long i = 0;
    int ref_len, ref_index;
    ref_index = start_ref_index;

    unsigned char *record_ref_hit;
    time_t t0 = time(0);

    long start_point = start;
    int buffer_size;
    //cout <<"Start slide ref..."<<endl;

    index_file.seekg (0, index_file.end);
    long length = index_file.tellg();
    //cout <<start_ref_index <<"\t"<< start_point<<endl;


    // /*
    index_file.seekg(start, ios::beg);
    while (!index_file.eof()){
        index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
        ref_len = real_index;
        // //cout <<start<< "####### each ref length" << ref_len << endl;
        buffer_size = (ref_len-k+1)*coder_num*4;
        char *each_ref_buffer = new char[buffer_size];
        index_file.read(each_ref_buffer, buffer_size);

        int n = 0;
        int j = 0;
        record_ref_hit = new unsigned char[ref_len*coder_num];
        record_ref_index = new unsigned int [ref_len*coder_num];
        while (n < buffer_size){
            memcpy(&kmer_index, each_ref_buffer, sizeof(unsigned int));
            record_ref_index[j] = kmer_index;
            if (kmer_index != 0){
                record_ref_hit[j] = kmer_count_table[kmer_index];
            }
            else{
                record_ref_hit[j] = 0;
            }
            j += 1;
            each_ref_buffer = each_ref_buffer + sizeof(unsigned int);
            n = n + sizeof(unsigned int);
        }
        // mtx.lock();
        slide_window(record_ref_hit, ref_len, ref_index, extract_ref_len, interval_file, hit_ratio,
                     perfect_hit_ratio,total_peak_num,MyPeak,record_ref_index,thread_index);
        // mtx.unlock();
        slide_ref_len += ref_len;
        start_point += 4;
        start_point += buffer_size;
        each_ref_buffer = each_ref_buffer - buffer_size; //move the pointer to the start
        delete [] each_ref_buffer;
        delete [] record_ref_hit;
        delete [] record_ref_index;
        if (ref_index % 10000 == 0){
            time_t t1 = time(0);
            //cout <<"#\t" <<start<<"\t"<<ref_index << "\t" <<slide_ref_len << " bp\t" <<
//                 extract_ref_len <<" bp\t" << total_peak_num << "peaks\t"<<t1-t0<< "\t single ref"<<ref_len << endl;
        }
        ref_index += 1;
        if (start_point >= end){
            // //cout << length << "\t" << start_point <<endl;
            // //cout << "read index ending." <<endl;
            break;
        }

        i += 1;
        real_index = 0;
    }
    // */
    index_file.close();
    interval_file.close();
    time_t t1 = time(0);
    //cout << slide_ref_len << " bp\t" << extract_ref_len <<" bp\t" <<total_peak_num<<"peaks\t" <<t1-t0<< endl;
    // //cout << "##########end" << end << "\t"<< start_point<<endl;
}

void read_fastq(string fastq_file, int k, char* coder, int* base, char* comple,
                short *choose_coder, int down_sam_ratio, long start, long end)
{
    time_t t0 = time(0);
    ifstream fq_file;
    fq_file.open(fastq_file);

    long pos = 0;
    for (long i = start; i>0; i--){
        fq_file.seekg(i, ios::beg);
        char j;
        fq_file.get(j);
        if (j == '@'){ //only read name has this symbol.
            pos = i;
            break;
        }
    }
    fq_file.seekg(pos, ios::beg);
    long add_size = start;


    string reads_seq;
    int reads_int [300];
    int reads_comple_int [300];

    unsigned int lines = 0;
    int converted_reads [900];
    int complemented_reads [900];
    int m;
    int n;
    unsigned int kmer_index, comple_kmer_index, real_index, b;
    int r ;
    short read_len = 0;

    // while (fq_file >> reads_seq)
    while (getline(fq_file,reads_seq))
    {
        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length();

        if (lines % 4 == 1){
            time_t t1 = time(0);
            if (lines % 10000000 == 10000000-1){
                //cout <<lines<<"reads\t" << t1-t0 <<endl;
            }
            // if (lines == 1){
            read_len = reads_seq.length();//cal read length
            // }
            r = rand() % 100 ;
            // if (lines < 100){
            //     //cout <<r << "r"<<endl;
            // }

            if (r < down_sam_ratio){
                for (int j = 0; j < read_len; j++){
                    reads_int[j] = (int)reads_seq[j];
                    reads_comple_int[j] = comple[reads_int[j]];
                }
                for (int j = 0; j < read_len-k+1; j++){
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        // //cout <<j<<"\t"<<i<<"\t";
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                            n = coder[c*choose_coder[(k-1-z)*3+i]+reads_comple_int[j+z]];
                            // //cout <<m<<" ";
                            // if (reads_int[j+z]==78 || reads_int[j+z]==110 ){
                            //     //cout << "N\t="<< (int)m <<endl;
                            // }
                            if (m == 5){
                                all_valid = false;
                                // //cout << "N !" << endl;
                                break;
                            }
                            kmer_index += m*base[z];
                            comple_kmer_index += n*base[(k-1-z)];

                        }
                        // //cout<<endl;
                        // //cout <<kmer_index<<"\t"<< j<<"\t"<<comple_kmer_index<<"\t"<<reads_seq[j] <<endl;
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }
                        else{
                            real_index = kmer_index;
                        }
                        if ((int)kmer_count_table[real_index] < least_depth & all_valid == true ){
                            kmer_count_table[real_index] += 1;
                            // //cout << (int)kmer_count_table[real_index] << "\t" << endl;
                        }
                    }
                }
            }
        }
        lines++;
    }
    fq_file.close();
    // return kmer_count_table;

}

char * generate_coder(bool coder_num)
{
    // A:65 97 T:116 84 C:99 67 G: 103 71
    static char coder [1000];
    for (int j = 0; j < 1000; j++){
        coder[j] = 5;
    }
    coder[65] = 1;
    coder[97] = 1;

    coder[84] = 1;
    coder[116] = 1;

    coder[67] = 0;
    coder[99] = 0;

    coder[71] = 0;
    coder[103] = 0;

    coder[65+c] = 1;
    coder[97+c] = 1;

    coder[84+c] = 0;
    coder[116+c] = 0;

    coder[67+c] = 1;
    coder[99+c] = 1;

    coder[71+c] = 0;
    coder[103+c] = 0;

    coder[65+2*c] = 1;
    coder[97+2*c] = 1;

    coder[84+2*c] = 0;
    coder[116+2*c] = 0;

    coder[67+2*c] = 0;
    coder[99+2*c] = 0;

    coder[71+2*c] = 1;
    coder[103+2*c] = 1;

    return coder;
}

int * generate_base(int k)
{
    static int base [32];
    for (int i = 0; i<k; i++){
        base[i] = pow(2, k-i-1);
    }
    return base;
}

char * generate_complement(void)
{
    static char comple [256];
    for (int j = 0; j < 256; j++){
        comple[j] = 0;
    }
    comple[65] = 84;
    comple[97] = 84;
    comple[116] = 65;
    comple[84] = 65;
    comple[99] = 71;
    comple[67] = 71;
    comple[103] = 67;
    comple[71] = 67;
    return comple;
}

short * random_coder(int k)
{
    static short permu[18] = {0,1,2,0,2,1,1,2,0,1,0,2,2,0,1,2,1,0};
    static short choose_coder[100];
    int r;
    unsigned seed;
    seed = time(0);
    srand(seed);
    for (int j = 0; j < k; j++){
        r = rand() % 6;
        // //cout << r << endl;
        for (int i = 0; i < coder_num; i++){
            choose_coder[j*3+i] = permu[r*3+i];
            // //cout << choose_coder[j*3+i]<<"#"<<permu[r*3+i] << endl;
        }
    }
    // for (int j = 0; j < 100; j++){
    //         //cout << j <<"\t" << choose_coder[j] << endl;
    // }
    return choose_coder;
}

short * saved_random_coder(string index_name){
    ifstream index_file;
    static short read_choose_coder[100];
    long i = 0;
    unsigned int real_index;
    index_file.open(index_name, ios::in | ios::binary);
    while (!index_file.eof()){
        index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
        if (i < 100){
            read_choose_coder[i] = real_index;
        }
        else{
            break;
        }
        i += 1;
    }
    index_file.close();
    return read_choose_coder;
}

int cal_sam_ratio(string fq1, long down_sampling_size){
    int down_sam_ratio = 0;
    int read_len = 0;
    long i = 0;
    long sample_size = 0;
    string reads_seq;

    ifstream fq_file;
    fq_file.open(fq1);
    while (getline(fq_file, reads_seq)){
        if (i == 1){
            read_len = reads_seq.length();
            //cout <<"The first read length is "<<read_len<<endl;
        }
        if (i % 4 == 1){
            sample_size += reads_seq.length();
        }
        i += 1;
    }
    fq_file.close();
    sample_size = sample_size * 2; //pair
    down_sam_ratio = 100*down_sampling_size/sample_size;
    //cout<<i<<"\t"<<sample_size<<" Down-sampling ratio is "<<down_sam_ratio<< "%." <<endl;
    return down_sam_ratio;
}

long file_size(string filename)
{
    struct stat statbuf;
    stat(filename.c_str(),&statbuf);
    long size=statbuf.st_size;
    return size;
}

long * split_ref(string index_name, string fasta_file, int thread_num){
    ifstream len_file;
    len_file.open(fasta_file+".genome.len.txt", ios::in);

    static long split_ref_cutoffs[300];
    int cut_index = 0;
    long index_size = file_size(index_name);
    long each_index_size = index_size/thread_num + 1;

    int ref_len, ref_index;
    long slide_ref_len;
    string chr_name;

    long add;
    long pos = 100 * 4;
    long start_byte, end_byte;
    start_byte = pos;
    long start_ref_index = 1;
    long count_ref_index = 0;

    while(len_file >> chr_name >> ref_index >> ref_len >> slide_ref_len ){
        count_ref_index += 1;

        add = 4*((ref_len-k+1)*coder_num+1); //the size of the genome.
        if (pos -start_byte > each_index_size){
            end_byte = pos + add;
            //cout << start_byte <<"--"<<end_byte<<"\t"<<ref_len<<"index"<< start_ref_index<<endl;
            split_ref_cutoffs[3*cut_index] = start_byte;
            split_ref_cutoffs[3*cut_index+1] = end_byte;
            split_ref_cutoffs[3*cut_index+2] = start_ref_index;
            cut_index += 1;
            start_byte = end_byte;
            start_ref_index = count_ref_index + 1;
        }
        pos += add;
        if (pos >= index_size){
            break;
        }

    }
    if (start_byte != index_size){
        end_byte = index_size;
        split_ref_cutoffs[3*cut_index] = start_byte;
        split_ref_cutoffs[3*cut_index+1] = end_byte;
        split_ref_cutoffs[3*cut_index+2] = start_ref_index;
        split_ref_cutoffs[3*cut_index+3] = 0; // indicate the end
        //cout << start_byte <<"\t--\t"<<end_byte<<"\t"<<ref_len<<"index"<< start_ref_index<<endl;
    }
    len_file.close();
    return split_ref_cutoffs;
}

int main( int argc, char *argv[])
{
    char *coder;
    int *base;
    char *comple;
    short *choose_coder;
    coder = generate_coder(3);
    base = generate_base(k);
    comple = generate_complement();
    time_t now1 = time(0);

    string fasta_file = argv[3];
    string fq1 = argv[1];
    string fq2 = argv[2];
    string interval_name = argv[4];
    string accept_hit_ratio = argv[5];
    string accept_perfect_hit_ratio = argv[6];
    string accept_thread_num = argv[7];
    float hit_ratio = stod(accept_hit_ratio);
    float perfect_hit_ratio = stod(accept_perfect_hit_ratio);
    long down_sampling_size = 2000000000; //2G bases
    thread_num = stod(accept_thread_num);


    long start = 0;
    long end = 0;

    unsigned seed;
    // seed = time(0);
    seed = 1;
    srand(seed);

    int down_sam_ratio = cal_sam_ratio(fq1, down_sampling_size); //percent of downsampling ratio (1-100).
    // int down_sam_ratio = 0;
    //index
    string index_name = fasta_file + ".k" + to_string(k) + ".index.dat";
    ifstream findex(index_name);
    if (! findex.good()){
        //cout << "Reference index not detected, start index..." << endl;
        choose_coder = random_coder(k);
        read_ref(fasta_file, coder, base, k, comple, index_name, choose_coder);
    }
    else{
        //cout << "Reference index is detected." << endl;
    }

    //cout << "Start extract..."<<endl;
    memset(kmer_count_table, 0, sizeof(char)*array_size);
    choose_coder = saved_random_coder(index_name);

    // fq1 = "/mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.05_depth50_reads150_sample_0_high.1.fq";
    // fq2 = "/mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.05_depth50_reads150_sample_0_high.2.fq";

    long size = file_size(fq1);
    long each_size = size/thread_num;
    // //cout <<size<<endl;

    std::vector<std::thread>threads;

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(read_fastq, fq1, k, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
    for (auto&th : threads)
        th.join();
    threads.clear();

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(read_fastq, fq2, k, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
    for (auto&th : threads)
        th.join();
    threads.clear();
    time_t now2 = time(0);
    //cout << "reads finish.\t" << now2 - now1 << endl;


    Peaks MyPeak;
    MyPeak.init(thread_num);
    memset(MyPeak.peak_filter, 0, sizeof(unsigned char)*MyPeak.max_peak_num);
    memset(MyPeak.peak_kmer, 0, sizeof(unsigned int)*array_size);

    int ref_thread_num;
    // if (thread_num > 5){
    //     ref_thread_num = 5;
    // }
    // else{
    //     ref_thread_num = thread_num;
    // }
    ref_thread_num = thread_num;
    long * split_ref_cutoffs = split_ref(index_name, fasta_file, ref_thread_num);
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    long total_peak_num = 0;
    // std::mutex mylock;
    for (int i=0; i<ref_thread_num; i++){
        start = split_ref_cutoffs[3*i];
        if (start == 0){
            break;
        }
        end = split_ref_cutoffs[3*i+1];
        int start_ref_index = split_ref_cutoffs[3*i+2];
        //cout << "Thread N.O."<<i << "\t"<< start << "\t" <<end<< "\t"<< start_ref_index << endl;
        // read_index(coder, base, k, comple, index_name, interval_name, choose_coder, hit_ratio, perfect_hit_ratio, std::ref(MyPeak), start, end);
        threads.push_back(thread(read_index, coder, base, k, comple, index_name, interval_name,
                                 choose_coder, hit_ratio, perfect_hit_ratio, std::ref(MyPeak), start, end,
                                 start_ref_index,std::ref(extract_ref_len),std::ref(slide_ref_len),std::ref(total_peak_num), i));
    }
    for (auto&th : threads)
        th.join();
    threads.clear();
    //cout << "raw peaks is done."<<endl;
    // down_sam_ratio = 50;

    // fq1 = "/mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.05_depth50_reads150_sample_0.1.fq";
    // fq2 = "/mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.05_depth50_reads150_sample_0.2.fq";

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        // //cout << "Reads thread N.O."<<i << "\t"<< start << "\t" <<end<< endl;
        threads.push_back(thread(&Peaks::slide_reads, MyPeak, fq1, fq2, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
    for (auto&th : threads)
        th.join();
    threads.clear();
    //cout << "filtering peaks is done."<<endl;

    MyPeak.count_filtered_peak(interval_name);
    //cout<<"filtered peak number:" <<MyPeak.filter_peak_num<<endl;
    MyPeak.delete_array();
    delete [] kmer_count_table;
    time_t now3 = time(0);
    //cout << "Finish with time:\t" << now3-now1<<endl;
    return 0;
}