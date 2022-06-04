#include <iostream>
#include <algorithm>
#include <string.h>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <set>
#include <vector>
#include <map>

using namespace std;
#pragma warning(disable : 4996)

//#define FILEOUT
#define SMALL_DATA //С���ݲ��Կ���
#define DEBUG
//#define OVERALL_DEBUG //�Ƚ�������
#define LOCAL_DEBUG //�Ƚϲ����˵���
#define PRINT_RESULT //�Ƿ�����������
//#define PRINT_LOOP_TIMES//������loop�Ĵ���
//#define PRINT_TERMINATE_DATA//���TERMINATE������
//#define PRINT_CAUGHT_FLOW_NUM//ÿ����100�������һ����Ϣ

// ���ѳ�
// Encoded in CRLF UTF-8

// ���� 8��byte
// ʱ��� 8��byte
#ifndef SMALL_DATA
const int ID_length = 8;
const int TimeStamp_length = 8;
#else
const int ID_length = 13;
const int TimeStamp_length = 0;
#endif // !SMALL_DATA

//---------------------   �ڴ˵���   --------------------------//

const int DATA_FILE_NUM = 10;     // Ҫ�����ļ�����
double POSSIBLE_THRESHOLD = 0.99; // hat_p����ֵ���������ṩ����0.99
const int STAR_THRESHOLD = 11;    // ���һ���������ʽ�г�������ô���*��������Ϊû�д���
int THRESHOLD = 10000;             // չʾ������ô��ļ�¼������

const double MY_ERROR_THRESHOLD_SKETCH = 2.0; // �����ֵ�߹���Сsketch����ô�౶������Ϊ�ܿ����Ǽ�����
const double MY_ERROR_THRESHOLD_V0 = 0.95;    // �����ֵ�߹���Сsketch����ô�౶������Ϊ�ܿ����Ǽ�����

// ��l,r,c������hash������ǰ,�Է���ʹ��
const int l = 8 * ID_length;            // ����bit��
const int r = 3;                        // sketch������
const int c = 9000;                     // sketch������
const int heavy_hash_length = 18000;    // heavy_flow�ĳ���
double eta = 50;                        // disagree / agree ����
double lowest_agree_ratio = 2.87915509; // ��һ����ͼ��������agree ratio

const double bit_ostracism_lamda = 0.01;
int bit_ostracism_0[l + 1], bit_ostracism_1[l + 1];
bool Is_kth_freed[l + 1] = { false };
int ostracism_est[l + 1];
//---------------------   �ڴ˵���   --------------------------//

char my_buggy_flow[ID_length] = {-114, -95, 23, 76, -58, 88, -37, 49};

int get_bit(unsigned char* a, int pos);

bool my_debug_cmp(char* s)
{
    for(int i = 0; i < ID_length; i++)
    {
        if(s[i] != my_buggy_flow[i])
        {
            return false;
        }
    }
    int buggy_bit = 40;
    printf("THE %d BIT IS: %d",buggy_bit, get_bit((unsigned char *)my_buggy_flow, buggy_bit - 1));
    return true;
}

bool my_cmp(char*, char*);
void Flow_out(char* s);

int num_of_star;

class ID_input
{
public:
    char x[ID_length + TimeStamp_length + 1];
    operator char* () const
    {
        return (char*)x;
    }
    operator unsigned char* () const
    {
        return (unsigned char*)x;
    }
    bool operator<(const ID_input& _f) const
    {
        for (int i = 0; i < ID_length; i++)
        {
            if (x[i] < _f.x[i])
            {
                return true;
            }
            else if (x[i] > _f.x[i])
            {
                return false;
            }
        }
        return false;
    }
};

// h1,h2...,hr �±��1��ʼ
uint32_t(*hash_function[r + 1])(char*);

uint64_t AwareHash(unsigned char* data, uint64_t n,
    uint64_t hash, uint64_t scale, uint64_t hardener)
{

    while (n)
    {
        hash *= scale;
        hash += *data++;
        n--;
    }
    return hash ^ hardener;
}

// ���Թ�ϣ��������������Ϊ���ѡȡ
// ��ϣ����Ϊ1 ~ c
uint32_t test_hash_0(char* f)
{
    return AwareHash((unsigned char*)f, ID_length, 354289553, 354289627, 1054289603) % c + 1;
}
uint32_t test_hash_1(char* f)
{
    return AwareHash((unsigned char*)f, ID_length, 554289569, 554289613, 2054289649) % c + 1;
}
uint32_t test_hash_2(char* f)
{
    return AwareHash((unsigned char*)f, ID_length, 654289577, 654289631, 2354289697) % c + 1;
}

uint32_t heavy_hash(char* f)
{
    return AwareHash((unsigned char*)f, ID_length, 354289577, 354289631, 1754289697) % heavy_hash_length + 1;
}

struct flow_tin
{
    ID_input f;
    uint32_t agree;
    uint32_t disagree;
    bool ex;
} heavy_flow[heavy_hash_length + 1];
bool heavy_flag[heavy_hash_length + 1];

vector<ID_input> all_id_flow;
vector<ID_input> smaller_id_flow;
// ��V[k][i][j]���� k = 0 Ϊ�ܵĲ� i j ����1��ʼ
unsigned int V[l + 1][r + 1][c + 1];
unsigned int V_initial[l + 1][r + 1][c + 1];
// ��ֵ
double p[l + 1];
double p_initial[l + 1];
// ��׼��,ע���Ƿ����
double sigma[l + 1];
double sigma_initial[l + 1];

void Insert_Flow(ID_input f)
{
    int h;

    h = heavy_hash(f);

    if (false == heavy_flag[h])
    {
        heavy_flow[h].f = f;
        heavy_flow[h].agree = 1;
        heavy_flow[h].ex = false;
        heavy_flow[h].disagree = 0;
        heavy_flag[h] = true;
    }
    else
    {
        if (my_cmp(f.x, heavy_flow[h].f.x))
        {
            heavy_flow[h].agree += 1;
        }
        else
        {
            heavy_flow[h].disagree += 1;
            if (((double)heavy_flow[h].disagree / (double)heavy_flow[h].agree) < eta)
            {
                smaller_id_flow.push_back(f);
            }
            else
            {
                if (heavy_flow[h].agree > 10000)
                {
                    printf("flow > 10000 but it is expel\n");
                }
                for (size_t i = 0; i < heavy_flow[h].agree; i++)
                {
                    smaller_id_flow.push_back(heavy_flow[h].f);
                }
                heavy_flow[h].f = f;
                heavy_flow[h].agree = 1;
                heavy_flow[h].ex = true;
                heavy_flow[h].disagree = 1;
            }
        }
    }
}

int k_count = 0; // ͳ��������
//.dat ת vector<ID>
int Read_Flowdata()
{
    char datafileName[100];
    // ע���ļ�·��
    int read_file_no;
#ifndef SMALL_DATA
    read_file_no = 1;
    sprintf(datafileName, "./formatted00.dat");
#else
    read_file_no = 11;
    sprintf(datafileName, "./data/0.dat");
#endif

    ID_input tmp_five_tuple;
    for(int i = 0; i < read_file_no; i++)
    {
        printf("READ %d FILES!\n", i);
        #ifdef SMALL_DATA
        sprintf(datafileName, "./data/%d.dat",i);
        #endif// SMALL_DATA
        FILE* fin = fopen(datafileName, "rb");
        if (NULL != fin)
        {
            int buggy_num = 0;
            fread(&tmp_five_tuple, TimeStamp_length, 1, fin);
            while (fread(&tmp_five_tuple, ID_length, 1, fin)) // ��8byte
            {
                // ����ʱ���
                Insert_Flow(tmp_five_tuple);
                all_id_flow.push_back(tmp_five_tuple);
                if(my_debug_cmp(tmp_five_tuple.x))
                {
                    printf("GOT BUGGY FLOW! SIZE: %d\n",++buggy_num);
                }
                k_count++;
                fread(&tmp_five_tuple, TimeStamp_length, 1, fin);
            }

            fclose(fin);

    #ifdef DEBUG
            int test_a = smaller_id_flow.size();
    #endif // DEBUG

        }
        else
        {
            return -1;
        }
    }
    return 0;
/*
#else
    for (int kk = 0; kk <= DATA_FILE_NUM; kk++)
    {
        sprintf(datafileName, "./data/%d.dat", kk);

        ID_input tmp_five_tuple;

        FILE* fin = fopen(datafileName, "rb");
        if (NULL != fin)
        {
            int k_count = 0;
            while (fread(&tmp_five_tuple, ID_length, 1, fin)) // ��13byte
            {
                smaller_id_flow.push_back(tmp_five_tuple);
                k_count++;
            }

            fclose(fin);

#ifdef DEBUG
            // Լ 2000 0000
            int test_a = smaller_id_flow.size();
            printf("READ %d FLOWS, %d TIMES\n", test_a, kk);
#endif // DEBUG
        }
        else
        {
            return -1;
        }
    }
    return 0;
    // sprintf(datafileName, "./data/all1.dat");
#endif // !SMALL_DATA
*/
}

//ʹ���˴�˷�
int get_bit(unsigned char* a, int pos)
{
    int byte = pos / 8;
    int bit = pos % 8;
    if ((a[byte] & (1 << (7 - bit))) == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}
void set_bit(unsigned char* a, int pos, int v)
{
    int byte = pos / 8;
    int bit = pos % 8;
    if (v == 1)
    {
        a[byte] = a[byte] | (1 << (7 - bit));
    }
    else
    {
        a[byte] = a[byte] & ~(1 << (7 - bit));
    }
}

void Flow2Sketch()
{
    memset(V, 0, sizeof(V));
    uint32_t tmp_hash[r + 1];
    for (ID_input tmp_flow : smaller_id_flow)
        // for (ID_input tmp_flow : all_id_flow)
    {
        for (size_t i = 1; i <= r; i++)
        {
            tmp_hash[i] = hash_function[i](tmp_flow);
            ++V[0][i][tmp_hash[i]];
        }
        for (size_t k = 1; k <= ID_length * 8; k++)
        {
            // 0 ��� V[k] ��Ӱ��
            if (0 != get_bit((unsigned char*)tmp_flow, k - 1))
            {
                bit_ostracism_1[k]++;
                for (size_t i = 1; i <= r; i++)
                {
                    ++V[k][i][tmp_hash[i]];
                }
            }
            else bit_ostracism_0[k]++;
        }
    }
    for (size_t k = 1; k <= ID_length * 8; k++)
    {
        double k_lamda = (double)bit_ostracism_0[k] / (double)(bit_ostracism_1[k] + 1);
        printf("bit_ostracism_0[%d]=%d, bit_ostracism_1[%d]=%d, k_lamda=%lf \n", (int)k, bit_ostracism_0[k], (int)k, bit_ostracism_1[k], k_lamda);
        if (k_lamda < bit_ostracism_lamda) {
            Is_kth_freed[k] = true;
            ostracism_est[k] = 1;
            printf("Is_kth_freed[%d] is true, ostracism_est[%d] is %d\n", (int)k,(int)k, ostracism_est[k]);
            //printf("bit_ostracism_0[%d]=%d, bit_ostracism_1[%d]=%d\n", (int)k, bit_ostracism_0[k], (int)k, bit_ostracism_1[k]);
        }
        if (k_lamda > 1.0 / bit_ostracism_lamda) {
            Is_kth_freed[k] = true;
            ostracism_est[k] = 0;
            printf("Is_kth_freed[%d] is true,ostracism_est[%d] is %d\n", (int)k, (int)k, ostracism_est[k]);
            //printf("bit_ostracism_0[%d]=%d,bit_ostracism_1[%d]=%d\n", (int)k, bit_ostracism_0[k], (int)k, bit_ostracism_1[k]);
        }
    }
    //printf("\nBit_Ostracism completed\n");
}

void Sketch2N_p_sigma()
{
    double sum;
    double square_sum;
    double tmp_r;
    for (size_t k = 0; k <= ID_length * 8; k++)
    {
        sum = 0;
        square_sum = 0;
        if (k > 0 && Is_kth_freed[k]) {
            for (size_t i = 1; i <= r; i++)
            {
                for (size_t j = 1; j <= c; j++)
                {
                    tmp_r = (ostracism_est[k] == 1) ? 1.0 : 0.0;
                    sum += tmp_r;
                    square_sum += tmp_r * tmp_r;
                }
            }
        }
        else {
            for (size_t i = 1; i <= r; i++)
            {
                for (size_t j = 1; j <= c; j++)
                {
                    tmp_r = (double)(V[k][i][j]) / (double)(V[0][i][j]);
                    sum += tmp_r;
                    square_sum += tmp_r * tmp_r;
                }
            }
        }
        p[k] = (double)sum / (double)(r * c);
        sigma[k] = sqrt(square_sum / (double)(r * c) - p[k] * p[k]);
    }
}

// ���˲���������bit_flow��prob_vector��ÿbit�����Ǵ�1��ʼ�ģ��������㷨�еı���
//         flow�Ǵ�0��ʼ�ģ������ϸ��������ݺ͹�ϣ����

double normalCFD(double value) // �����׼��̬�ֲ�
{
    /*
    printf("------------------------------VAL1: %lf\n", value);
    printf("------------------------------VAL2: %lf\n", -value / sqrt(2));
    printf("------------------------------RSLT: %lf\n", erfc(-value / sqrt(2)));
    */
    return 0.5 * erfc(-value / sqrt(2));
}

struct ans_t // ��extract large flow����������Ԫ�ص�type
{
    char bit_flow[l + 2];       // �������1��byte��'0'��'1'��ʾ1��bit
    char flow[(l + 7) / 8 + 1]; // ����Ǹ�Դ������ͬ�Ľṹ
    unsigned int size;          // ���Ĵ�С
    double prob_vector[l + 2];  // ����������
    ans_t(char* bbit_flow, char* fflow, unsigned int ssize = 0, double* pprob_vector = NULL)
    {
        strcpy(bit_flow, bbit_flow);
        for (int i = 0; i < ID_length; i++)
        {
            flow[i] = fflow[i];
        }
        if (pprob_vector != NULL)
        {
            for (int i = 1; i <= l; i++)
            {
                prob_vector[i] = pprob_vector[i];
            }
        }
        size = ssize;
    }
};

double cal_hat_p(double theta, int i, int j,
    unsigned int V[][r + 1][c + 1], double* p, double* sigama, int k) // ���������k��bitΪ1�ĸ���
{
    if (Is_kth_freed[k]) {
        if (ostracism_est[k] == 0)return 0;
        return 1;
    }
    double r = (double)V[k][i][j] / V[0][i][j];
    if (r < theta)
    {
        return 0;
    }
    if (1 - r < theta)
    {
        return 1;
    }
    double ans = 0;
    double prob_1 = (V[k][i][j] - theta * V[0][i][j]) /
        (V[0][i][j] - theta * V[0][i][j]);
    double prob_0 = (V[k][i][j]) /
        (V[0][i][j] - theta * V[0][i][j]);
    double normal_val1 = normalCFD((prob_1 - p[k]) / sigama[k]);
    double normal_val0 = normalCFD((prob_0 - p[k]) / sigama[k]);
    /*
    if (normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) < 0.99 && normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) > 0.01 ||
        normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) < 0 || normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) > 1)
    {
        printf("Sigama: %lf, Pk : %lf, i: %d, j: %d, theta: %lf, Vk: %d, V0: %d\n", sigama[k], p[k], i, j, theta, V[k][i][j], V[0][i][j]);
        printf("PROB1 : %lf, PROB2 : %lf, N0 : %lf, N1: %lf, HAT_P : %lf\n", prob_0, prob_1, normal_val0, normal_val1, normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]));
    }
    */
    return normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]);
}

double cal_hat_p_print(double theta, int i, int j,
    unsigned int V[][r + 1][c + 1], double* p, double* sigama, int k) // ���������k��bitΪ1�ĸ���
{
    if (Is_kth_freed[k]) {
        if (ostracism_est[k] == 0)return 0;
        return 1;
    }
    double r = (double)V[k][i][j] / V[0][i][j];
    if (r < theta)
    {
        return 0;
    }
    if (1 - r < theta)
    {
        return 1;
    }
    double ans = 0;
    double prob_1 = (V[k][i][j] - theta * V[0][i][j]) /
        (V[0][i][j] - theta * V[0][i][j]);
    double prob_0 = (V[k][i][j]) /
        (V[0][i][j] - theta * V[0][i][j]);
    double normal_val1 = normalCFD((prob_1 - p[k]) / sigama[k]);
    double normal_val0 = normalCFD((prob_0 - p[k]) / sigama[k]);
    if (normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) < 0.99 && normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) > 0.01 ||
        normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) < 0 || normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]) > 1)
    {
        printf("ERFC1: %lf, ERFC0: %lf\n", erfc(-((prob_1 - p[k]) / sigama[k]) / sqrt(2)), erfc(-((prob_0 - p[k]) / sigama[k]) / sqrt(2)));
        printf("NORMAL ARG1: %lf, ARG0: %lf\n", (prob_1 - p[k]) / sigama[k], (prob_0 - p[k]) / sigama[k]);
        printf("Sigama: %lf, Pk : %lf, i: %d, j: %d, theta: %lf, Vk: %d, V0: %d\n", sigama[k], p[k], i, j, theta, V[k][i][j], V[0][i][j]);
        printf("PROB1 : %lf, PROB2 : %lf, N0 : %lf, N1: %lf, HAT_P : %lf\n", prob_0, prob_1, normal_val0, normal_val1, normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]));
    }
    return normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]);
}

struct two_types_of_flow // ��������ʾ
{
    char bit_flow[l + 2];       // һ��'1'����'0'��ʾһ��bit
    char flow[(l + 7) / 8 + 1]; // Ӧ���Ǻ͸���������һ���Ľṹ(��֪���᲻����bug...)
    two_types_of_flow(char* bbit_flow, char* fflow)
    {
        strcpy(bit_flow, bbit_flow);
        for (int i = 0; i < (l + 7) / 8; i++)
        {
            flow[i] = fflow[i];
        }
        flow[(l + 7) / 8] = '\0';
    }
};

vector<two_types_of_flow> possible_flows; // һ������������ȫ�ֱ���
char current_T[l + 2];                    // һ������������ȫ�ֱ���

int flag_ = 0;

int loop_num = 0;

void find_possible_flows(int i, int j, int k, char* T) // �ҵ��������ʽ�����п��ܵ���
{
    if (k == l + 1)
    {
        if ((++loop_num) % 15000 == 0 && loop_num > 15000)
        {
            printf("%d stars: LOOP %d TIMES!\n", num_of_star, loop_num);
        }
        char ans[(l + 7) / 8 + 1];

        for (int kk = 1; kk <= l; kk++)
        {
            set_bit((unsigned char*)ans, kk - 1, current_T[kk] == '1' ? 1 : 0);
        }

        ans[(l + 7) / 8] = '\0';
        if (hash_function[i](ans) == j)
        {
            int flag = 0;
            /*
            for (auto iter : smaller_id_flow)
            {
                if (my_cmp(iter.x, ans))
                {
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                printf("----------------NOT IN!-------------\n");
                printf("T: %s\n", T);
            }
            else if (flag == 1)
            {
                printf("----------------IN!-------------\n");
                Flow_out(ans);
                printf("\n");
            }
            */
            possible_flows.push_back(two_types_of_flow(current_T, ans));
            // Flow_out(possible_flows.back().flow);
        }
        return;
    }
    else
    {
        if (T[k] != '*')
        {
            current_T[k] = T[k];
            find_possible_flows(i, j, k + 1, T);
        }
        else
        {
            current_T[k] = '0';
            find_possible_flows(i, j, k + 1, T);
            current_T[k] = '1';
            find_possible_flows(i, j, k + 1, T);
        }
    }
    return;
}

/* ����ans_t��һ��������ÿ��Ԫ�ش洢�ˣ�
 *
 *      ��byte��ʾbit������
 *      ����������ͬ�ṹ������
 *      ���Ĵ�С
 *      ����������
 *
 */
vector<ans_t> ExtractLargeFlows(double theta, int i, int j,
    unsigned int V[][r + 1][c + 1], double* p, double* sigama)
{
    // ��һ��������ÿ��bit�ĸ��ʹ�ֵ
    double hat_p[l + 1];
    for (int k = 1; k <= l; k++)
    {
        hat_p[k] = cal_hat_p(theta, i, j, V, p, sigama, k);
    }
    // printf("V[%d][%d] step 1 completed\n", i, j);
    //  �ڶ������ҵ����к�ѡ�Ĵ���������possible_flows����
    char T[l + 2];
    for (int k = 1; k <= l; k++)
    {
        if (hat_p[k] > POSSIBLE_THRESHOLD)
        {
            T[k] = '1';
        }
        else if (1 - hat_p[k] > POSSIBLE_THRESHOLD)
        {
            T[k] = '0';
        }
        else
        {
            T[k] = '*';
        }
    }
    T[l + 1] = '\0';
    T[0] = '#';
    num_of_star = 0;
    loop_num = 0;
    for (int kk = 1; kk <= l; kk++)
    {
        if (T[kk] == '*')
        {
            num_of_star++;
        }
    }
    vector<ans_t> result;
    if (num_of_star > STAR_THRESHOLD)
    {
        return result;
        for (int k = 1; k <= l; k++)
        {
            hat_p[k] = cal_hat_p_print(theta, i, j, V, p, sigama, k);
        }
    }
    current_T[l + 1] = '\0';
    current_T[0] = '#';
    possible_flows.clear();
    find_possible_flows(i, j, 1, T);
    // printf("V[%d][%d] step 2 completed\n", i, j);
    //  �����������ƴ�����Ƶ�ʺͿ���������
    vector<double> estimated_frequency;
    double estimated_p[l + 1];
    for (vector<two_types_of_flow>::iterator item = possible_flows.begin();
        item != possible_flows.end(); item++)
    {
        estimated_frequency.clear();
        int min_sketch = 0xfffffff;
        for (int k = 1; k <= l; k++)
        {
            if (Is_kth_freed[k]) {
                if (ostracism_est[k] == 1) {
                    estimated_p[k] = hat_p[k];
                }
                else if (ostracism_est[k] == 0) {
                    estimated_p[k] = 1 - hat_p[k];
                }
            }
            else if (!Is_kth_freed[k] && item->bit_flow[k] == '1')
            {
                min_sketch = (V[k][i][j] < min_sketch) ? V[k][i][j] : min_sketch;
                double r = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency.push_back(((r - p[k]) / (1 - p[k])) * V[0][i][j]);
                estimated_p[k] = hat_p[k];
            }
            else
            {
                min_sketch = (V[0][i][j] - V[k][i][j] < min_sketch) ? V[0][i][j] - V[k][i][j] : min_sketch;
                double r = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency.push_back((1 - r / p[k]) * V[0][i][j]);
                estimated_p[k] = 1 - hat_p[k];
            }
        }
        sort(estimated_frequency.begin(), estimated_frequency.end());
        double ans_estimated_frequency = estimated_frequency[(estimated_frequency.size() - 1) / 2];
        if (ans_estimated_frequency > min_sketch)
        {
            if (ans_estimated_frequency > MY_ERROR_THRESHOLD_SKETCH * min_sketch &&
                ans_estimated_frequency > MY_ERROR_THRESHOLD_V0 * V[0][i][j])
            {
                break;
            }
            ans_estimated_frequency = min_sketch;
            /*
            printf("ESTIMATED LARGER THAN MIN SKETCH! ESTIMATED: %lf, MIN: %d, V0: %d\n", ans_estimated_frequency, min_sketch, V[0][i][j]);
            ans_estimated_frequency = min_sketch;
            */

            /*
            int flag = 0;
            for (auto iter : smaller_id_flow)
            {
                if (my_cmp(iter.x, item->flow))
                {
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                printf("----------------NOT IN!-------------\n");
                printf("T: %s\n", T);
            }
            */
        }
        // Flow_out(item->flow);
        result.push_back(ans_t(item->bit_flow, item->flow, ans_estimated_frequency, estimated_p));
        // Flow_out(result.back().flow);
    }
    // printf("V[%d][%d] step 3 completed\n", i, j);
    //  ���Ĳ���ȥsketch����ѡ�������ݣ�ɾ����С��
    if (r == 1)
    {
        return result;
    }
    for (vector<ans_t>::iterator item = result.begin(); item != result.end();)
    {
        for (int ii = 1; ii <= r; ii++)
        {
            if (ii == i)
            {
                continue;
            }
            int jj = hash_function[ii](item->flow);
            for (int k = 1; k <= l; k++)
            {
                if (Is_kth_freed[k]) {
                    continue;
                }
                else if (item->bit_flow[k] == '0' && V[0][ii][jj] - V[k][ii][jj] < item->size)
                {
                    item->size = V[0][ii][jj] - V[k][ii][jj];
                }
                else if (item->bit_flow[k] == '1' && V[k][ii][jj] < item->size)
                {
                    item->size = V[k][ii][jj];
                }
            }
        }
        if (item->size < theta * V[0][i][j])
        {
            item = result.erase(item);
            if (item == result.end())
                break;
        }
        else
        {
            item++;
        }
    }
    return result;
}

//��sketch��ɾȥ����
void RemoveFlows(vector<ans_t> FF)
{
    int FF_size = FF.size();
    uint32_t tmp_hash[r + 1];
    for (int it = 0; it < FF.size(); it++)
    {
        char* ans = FF[it].flow;

        for (size_t i = 1; i <= r; i++)
        {
            tmp_hash[i] = hash_function[i](ans);
            V[0][i][tmp_hash[i]] -= FF[it].size;
        }
        for (size_t k = 1; k <= ID_length * 8; k++)
        {
            if (0 != get_bit((unsigned char*)ans, k - 1) && !Is_kth_freed[k])
            {
                for (size_t i = 1; i <= r; i++)
                {
                    V[k][i][tmp_hash[i]] -= FF[it].size;
                }
            }
        }
    }
}

//��鵱ǰsketch�Ƿ���ϸ�˹�ֲ�
//�������������L��sketchÿһ��������1sigma,2sigma,3sigma������Ҫ��
bool Terminate(double theta)
{
    double STEP = 0.005;
    double RATE1 = 0.6826 + STEP * log2(theta);
    double RATE2 = 0.9544 + STEP * log2(theta);
    double RATE3 = 0.9973 + STEP * log2(theta);

    for (size_t k = 1; k <= ID_length * 8; k++)
    {
        size_t sigma_num1 = 0, sigma_num2 = 0, sigma_num3 = 0;
        for (int i = 1; i <= r; i++)
        {
            for (int j = 1; j <= c; j++)
            {
                double r;
                if (Is_kth_freed[k]) {
                    sigma_num3++;
                    sigma_num2++;
                    sigma_num1++;
                }
                else if (!Is_kth_freed[k]) {
                    r = (double)V[k][i][j] / V[0][i][j];
                    if (r <= p[k] + 3.0 * sigma[k] && r >= p[k] - 3.0 * sigma[k])
                        sigma_num3++;
                    if (r <= p[k] + 2.0 * sigma[k] && r >= p[k] - 2.0 * sigma[k])
                        sigma_num2++;
                    if (r <= p[k] + 1.0 * sigma[k] && r >= p[k] - 1.0 * sigma[k])
                        sigma_num1++;
                }
            }
        }
        if (sigma_num1 == 0 && sigma_num2 == 0 && sigma_num3 == 0)
        {
            printf("I WONDER WHY PROGRAME RUNNING REACH HERE\n");
        }
        double rate1 = (double)sigma_num1 / (double)(r * c);
        double rate2 = (double)sigma_num2 / (double)(r * c);
        double rate3 = (double)sigma_num3 / (double)(r * c);
#ifdef PRINT_TERMINATE_DATA
        printf("V[%d] sigma1_num=%d,sigma2_num=%d,sigma3_num=%d\n", (int)k, (int)sigma_num1, (int)sigma_num2, (int)sigma_num3);
        printf("V[%d] rate1=%lf,rate2=%lf,rate3=%lf\n", (int)k, rate1, rate2, rate3);
#endif // PRINT_TERMINATE_DATA
        if (rate1 < RATE1)
            return false;
        if (rate2 < RATE2)
            return false;
        if (rate3 < RATE3)
            return false;
    }
    return true;
}

bool my_cmp(char* s1, char* s2)
{
    for (size_t i = 0; i < ID_length; i++)
    {
        if (s1[i] != s2[i])
        {
            return false;
        }
    }
    return true;
}

void Flow_out(char* s)
{
    printf("ID: ");
    for (size_t i = 0; i < ID_length; i++)
    {
        printf("%d ", s[i]);
    }
    printf("\n");
}

struct two_int
{
    int i1;
    int i2;
    double ratio;
    two_int() : i1(0), i2(0) {};
};

int main()
{
#ifdef FILEOUT
    freopen("out.txt", "w", stdout);
#endif // FILEOUT

    // �����ϣ����
    hash_function[1] = test_hash_0;
    hash_function[2] = test_hash_1;
    hash_function[3] = test_hash_2;
    // FΪ���д������ϣ�FFΪÿ��ѭ���ҳ��Ĵ�������
    vector<ans_t> F;

    int lowest_agree = (int)(((double)k_count / (double)heavy_hash_length) * lowest_agree_ratio); // ��һ����ͼ��������agree number

    if (0 == Read_Flowdata()) //�����ݶ���
    {
        {
            char tmp_bit_flow[l + 2];
            double tmp_prob_vector[l + 2];
            uint32_t tmp_sk_n;
            for (size_t k = 1; k <= l; k++)
            {
                tmp_prob_vector[k] = 1;
            }

            for (size_t j = 1; j <= heavy_hash_length; j++)
            {
                if (my_debug_cmp(heavy_flow[j].f.x))
                {
                    printf("HEAVY PART CAUGHT THE BUGGY FLOW, SIZE: %d\n", heavy_flow[j].agree);
                }
                if (lowest_agree > heavy_flow[j].agree)
                {
                    for (size_t i = 0; i < heavy_flow[j].agree; i++)
                    {
                        smaller_id_flow.push_back(heavy_flow[j].f);
                    }
                    continue;
                }
                for (size_t tmp_l = 1; tmp_l <= l; tmp_l++)
                {
                    tmp_bit_flow[tmp_l] = (get_bit(heavy_flow[j].f, tmp_l - 1) == 1) ? '1' : '0';
                }
                F.push_back(ans_t(tmp_bit_flow, heavy_flow[j].f, heavy_flow[j].agree, tmp_prob_vector));
            }
        }
        // sketch ����N(p,sigma)
        Flow2Sketch();
        // ����heavy_flow
        Sketch2N_p_sigma();

        printf("Read in flow data success\n");
        memcpy(V_initial, V, sizeof(V));
        memcpy(p_initial, p, sizeof(p));
        memcpy(sigma_initial, sigma, sizeof(sigma));
    }
    else
    {
        printf("Read in flow data error\n");
    }

    double theta = 0.5;
    int nnnn = 1;
    while (1)
    {
        printf("LOOP %d START!\n", nnnn);
        for (auto iter : F)
        {
            if (my_debug_cmp(iter.flow))
            {
                printf("BUGGY FLOW IS IN! SIZE: %d\n", iter.size);
            }
        }
        int my_flow_num = 0;
        vector<ans_t> FF;
        for (int i = 1; i <= r; i++)
        {
            for (int j = 1; j <= c; j++)
            {
                if (0 == V[0][i][j])
                {
                    continue;
                }
                vector<ans_t> temp_F = ExtractLargeFlows(theta, i, j,
                    V, p, sigma);
                if (!temp_F.empty())
                {
                    my_flow_num++;
                    for (vector<ans_t>::iterator it = temp_F.begin(); it < temp_F.end(); it++)
                    {
                        bool temp_Fin = false;
                        if (!FF.empty())
                        {
                            for (vector<ans_t>::iterator iter = FF.begin(); iter < FF.end(); iter++)
                            {
                                if (strcmp(iter->bit_flow, it->bit_flow) == 0)
                                {
                                    temp_Fin = true;
                                    break;
                                }
                            }
                        }
                        if (!temp_Fin)
                            FF.push_back(*it);
                    }
#ifdef PRINT_CAUGHT_FLOW_NUM
                    if (my_flow_num % 100 == 0)
                    {
                        printf("CAUGHT %d FLOWS!\n", my_flow_num);
                    }
#endif
                }
            }
        }
        //����ѭ���ҳ�����ʱ���޳����������¼�������������
        if (!FF.empty())
        {
            printf("SIZE : %d\n", FF.size());
            for (vector<ans_t>::iterator it = FF.begin(); it < FF.end(); it++)
            {
                bool FF_in = false;
                vector<ans_t>::iterator temp_pos = FF.begin();
                if (!F.empty())
                {
                    for (vector<ans_t>::iterator iter = F.begin(); iter < F.end(); iter++)
                    {
                        if (my_cmp(iter->flow, it->flow))
                            // if (strcmp(iter->bit_flow, it->bit_flow) == 0)
                        {
                            FF_in = true;
                            temp_pos = iter;
                            break;
                        }
                    }
                }
                if (!FF_in)
                    F.push_back(*it);
                else if (FF_in)
                {
                    temp_pos->size += it->size;
                    if (my_debug_cmp(temp_pos->bit_flow))
                    {
                        printf("ADD BUGGY FLOW FINISH, SIZE: %d\n", temp_pos->size);
                    }
                }
            }
            RemoveFlows(FF);
            printf("RemoveFlowscompleted\n");
            Sketch2N_p_sigma();
        }
        printf("%d loop is completed______________, theta = %lf\n\n", nnnn, theta);
        nnnn++;
        /*
        if (nnnn > 10)
            break;
        */

        if (Terminate(theta))
            break;
        //û���ҳ�������theta����
        if (FF.empty())
            theta /= 2;
    }

#ifdef DEBUG

    printf("\n--------------------------------------------------------\n");
#endif // DEBUG
    printf("\n########################################################\n");

#ifdef DEBUG
    {
        class flow_debug
        {
        public:
            ID_input f;
            uint32_t s;
            bool operator<(const flow_debug& _flow)
            {
                return (s < _flow.s);
            }
        };
        map<ID_input, two_int> flow_queue;
        vector<ID_input> Flow_sort = all_id_flow;
        printf("SORT BEGIN!\n");
        sort(Flow_sort.begin(), Flow_sort.end());
        printf("SORT END!\n");
        printf("NEED %d LOOPS\n", Flow_sort.size());
        ID_input tmp_flow;
        uint32_t flow_size;
        auto iter = Flow_sort.begin();
        tmp_flow = *Flow_sort.begin();
        flow_size = 1;
        ++iter;

        int loop_time = 0;

        flow_debug tmp;
        while (iter < Flow_sort.end())
        {
            if (my_cmp(tmp.f.x, iter->x))
            {
                ++flow_size;
            }
            else
            {
                tmp.f = tmp_flow;
                tmp.s = flow_size;
                flow_queue[tmp.f].i1 = tmp.s;
                tmp_flow = *iter;
                flow_size = 1;
            }
            ++iter;
            loop_time++;
#ifdef PRINT_LOOP_TIMES
            if (loop_time % 100000 == 0)
            {
                printf("LOOP %d TIMES!\n", loop_time);
                printf("FLOW QUEUE LENGTH: %d\n", flow_queue.size());
            }
#endif // PRINT_LOOP_TIMES
        }
        flow_queue[tmp.f].i1 = flow_size;
        for (auto item : F)
        {
            ID_input x;
            for (int i = 0; i <= ID_length; i++)
            {
                x.x[i] = item.flow[i];
            }
            flow_queue[x].i2 = item.size;
            if (flow_queue[x].i2 == 797 && flow_queue[x].i1 == 15119)
            {
                Flow_out(x.x);
            }
            flow_queue[x].ratio = (double)flow_queue[x].i1 / flow_queue[x].i2;
        }
        double error_rate = 0;
        for (auto i : flow_queue)
        {
            if (i.second.i1 > THRESHOLD || i.second.i2 > THRESHOLD)
            {
                if(i.second.i1 == 18031)
                {
                    Flow_out(i.first);
                }
                // Flow_out(i.first);
#ifdef PRINT_RESULT
                printf("appear  %d  times, SKETCH CATCH %d TIMES, RATIO: %lf\n", i.second.i1, i.second.i2, i.second.ratio);
#endif // PRINT_RESULT
                error_rate += i.second.i1 * (i.second.ratio - 1.0) * (i.second.ratio - 1.0);
            }
        }

        printf("\nDEBUG END!!! \nTHRESHOLD: %d\nERROR RATIO: %lf\nPOSSIBLE_THRESHOLD: %lf\nSTAR_THRESHOLD: %d\nMY_ERROR_THRESHOLD_SKETCH: %lf\nMY_ERROR_THRESHOLD_V0: %lf\nSKETCH LENGTH: %d, HEIGHT: %d, WIDTH: %d\n",
            THRESHOLD, error_rate, POSSIBLE_THRESHOLD, STAR_THRESHOLD, MY_ERROR_THRESHOLD_SKETCH, MY_ERROR_THRESHOLD_V0, l, r, c);
        return 0;
    }
#endif // DEBUG
    return 0;
}