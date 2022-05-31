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
#define SMALL_DATA //小数据测试开关
#define DEBUG
//#define OVERALL_DEBUG //比较所有流
#define LOCAL_DEBUG //比较捕获了的流

// 兰佳晨
// Encoded in CRLF UTF-8

// 流键 13个byte
// 时间戳 13个byte
#ifndef SMALL_DATA
const int ID_length = 8;
const int TimeStamp_length = 8;
#else
const int ID_length = 13;
const int TimeStamp_length = 0;
#endif // !SMALL_DATA

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

// 将l,r,c参数及hash函数提前,以方便使用
const int l = 8 * ID_length;
const int r = 3;   // 2是我瞎写的
const int c = 5000; // 100是我瞎写的

const double MY_ERROR_THRESHOLD_SKETCH = 2.0; // 如果估值高过最小sketch的这么多倍，则认为是假阳性
const double MY_ERROR_THRESHOLD_V0 = 0.95; // 如果估值高过最小sketch的这么多倍，则认为是假阳性

// h1,h2...,hr 下标从1开始
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

// 测试哈希函数，六个质数为随机选取
// 哈希返回为1 ~ c
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

vector<ID_input> all_id_flow;
// 按V[k][i][j]排列 k = 0 为总的层 i j 都从1开始
unsigned int V[l + 1][r + 1][c + 1];
unsigned int V_initial[l + 1][r + 1][c + 1];
// 均值
double p[l + 1];
double p_initial[l + 1];
// 标准差,注意是方差开方
double sigma[l + 1];
double sigma_initial[l + 1];

//.dat 转 vector<ID>
int Read_Flowdata()
{
    char datafileName[100];
    // 注意文件路径
#ifndef SMALL_DATA
    sprintf(datafileName, "./formatted00.dat"); 

    ID_input tmp_five_tuple;

    FILE* fin = fopen(datafileName, "rb");
    if (NULL != fin)
    {
        int k_count = 0;
        while (fread(&tmp_five_tuple, ID_length, 1, fin)) // 读13byte
        {
            all_id_flow.push_back(tmp_five_tuple);
            // 跳过时间戳
            fread(&tmp_five_tuple, TimeStamp_length, 1, fin);
            k_count++;
        }

        fclose(fin);

#ifdef DEBUG
        // 约 2000 0000
        int test_a = all_id_flow.size();
#endif // DEBUG

        return 0;
    }
    else
    {
        return -1;
    }
#else
    for(int kk = 0; kk <= 10; kk++)
    {
        sprintf(datafileName, "./data/%d.dat",kk);

        ID_input tmp_five_tuple;

        FILE* fin = fopen(datafileName, "rb");
        if (NULL != fin)
        {
            int k_count = 0;
            while (fread(&tmp_five_tuple, ID_length, 1, fin)) // 读13byte
            {
                all_id_flow.push_back(tmp_five_tuple);
                // 跳过时间戳
                fread(&tmp_five_tuple, TimeStamp_length, 1, fin);
                k_count++;
            }

            fclose(fin);

    #ifdef DEBUG
            // 约 2000 0000
            int test_a = all_id_flow.size();
            printf("READ %d FLOWS, %d TIMES\n",test_a, kk);
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
}

//使用了大端法
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
    uint32_t tmp_hash[r + 1];
    for (ID_input tmp_flow : all_id_flow)
    {
        for (size_t i = 1; i <= r; i++)
        {
            tmp_hash[i] = hash_function[i](tmp_flow);
            ++V[0][i][tmp_hash[i]];
        }
        for (size_t k = 1; k <= ID_length * 8; k++)
        {
            // 0 则对 V[k] 无影响
            if (0 != get_bit((unsigned char*)tmp_flow, k - 1))
            {
                for (size_t i = 1; i <= r; i++)
                {
                    ++V[k][i][tmp_hash[i]];
                }
            }
        }
    }
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
        for (size_t i = 1; i <= r; i++)
        {
            for (size_t j = 1; j <= c; j++)
            {
                tmp_r = (double)(V[k][i][j]) / (double)(V[0][i][j]);
                sum += tmp_r;
                square_sum += tmp_r * tmp_r;
            }
        }
        p[k] = (double)sum / (double)(r * c);
        sigma[k] = sqrt(square_sum / (double)(r * c) - p[k] * p[k]);
    }
}

// 杨仕博：代码里bit_flow、prob_vector、每bit概率是从1开始的，以契合算法中的表述
//         flow是从0开始的，以契合给定的数据和哈希函数

double normalCFD(double value) // 计算标准正态分布
{
    /*
    printf("------------------------------VAL1: %lf\n", value);
    printf("------------------------------VAL2: %lf\n", -value / sqrt(2));
    printf("------------------------------RSLT: %lf\n", erfc(-value / sqrt(2)));
    */
    return 0.5 * erfc(-value / sqrt(2));
}

struct ans_t // 是extract large flow返回向量中元素的type
{
    char bit_flow[l + 2];       // 这个是用1个byte的'0'或'1'表示1个bit
    char flow[(l + 7) / 8 + 1]; // 这个是跟源数据相同的结构
    unsigned int size;          // 流的大小
    double prob_vector[l + 2];  // 可能性向量
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
    unsigned int V[][r + 1][c + 1], double* p, double* sigama, int k) // 计算大流第k个bit为1的概率
{
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
    unsigned int V[][r + 1][c + 1], double* p, double* sigama, int k) // 计算大流第k个bit为1的概率
{
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


struct two_types_of_flow // 两种流表示
{
    char bit_flow[l + 2];       // 一个'1'或者'0'表示一个bit
    char flow[(l + 7) / 8 + 1]; // 应该是和给定的数据一样的结构(不知道会不会有bug...)
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

vector<two_types_of_flow> possible_flows; // 一个辅助函数的全局变量
char current_T[l + 2];                    // 一个辅助函数的全局变量

int flag_ = 0;

int loop_num = 0;

void find_possible_flows(int i, int j, int k, char* T) // 找到正则表达式中所有可能的流
{
    if (k == l + 1)
    {
        if ((loop_num++) % 10000 == 0 && loop_num > 100000)
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
            for (auto iter : all_id_flow)
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
            //Flow_out(possible_flows.back().flow);
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

/* 返回ans_t的一个向量，每个元素存储了：
 *
 *      用byte表示bit的流键
 *      和输入数据同结构的流键
 *      流的大小
 *      可能性向量
 *
 */
vector<ans_t> ExtractLargeFlows(double theta, int i, int j,
    unsigned int V[][r + 1][c + 1], double* p, double* sigama)
{
    // 第一步，计算每个bit的概率估值
    double hat_p[l + 1];
    for (int k = 1; k <= l; k++)
    {
        hat_p[k] = cal_hat_p(theta, i, j, V, p, sigama, k);
    }
    // printf("V[%d][%d] step 1 completed\n", i, j);
    //  第二步，找到所有候选的大流，存在possible_flows里面
    char T[l + 2];
    for (int k = 1; k <= l; k++)
    {
        if (hat_p[k] > 0.99)
        {
            T[k] = '1';
        }
        else if (1 - hat_p[k] > 0.99)
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
    //printf("There are %d stars\n", num_of_star);
    vector<ans_t> result;
    if (num_of_star > 15)
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
    //  第三步，估计大流的频率和可能性向量
    double estimated_frequency[l + 1];
    double estimated_p[l + 1];
    for (vector<two_types_of_flow>::iterator item = possible_flows.begin();
        item != possible_flows.end(); item++)
    {
        int min_sketch = 0xfffffff;
        for (int k = 1; k <= l; k++)
        {
            if (item->bit_flow[k] == '1')
            {
                min_sketch = (V[k][i][j] < min_sketch)? V[k][i][j] : min_sketch;
                double r = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency[k] = ((r - p[k]) / (1 - p[k])) * V[0][i][j];
                estimated_p[k] = hat_p[k];
            }
            else
            {
                min_sketch = (V[0][i][j] - V[k][i][j] < min_sketch)? V[0][i][j] - V[k][i][j] : min_sketch;
                double r = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency[k] = (1 - r / p[k]) * V[0][i][j];

                estimated_p[k] = 1 - hat_p[k];
            }
        }
        sort(estimated_frequency + 1, estimated_frequency + 1 + l);
        double ans_estimated_frequency = estimated_frequency[l / 2];
        if(ans_estimated_frequency > min_sketch)
        {
            if(ans_estimated_frequency > MY_ERROR_THRESHOLD_SKETCH * min_sketch && 
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
            for (auto iter : all_id_flow)
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
        //Flow_out(item->flow);
        result.push_back(ans_t(item->bit_flow, item->flow, ans_estimated_frequency, estimated_p));
        //Flow_out(result.back().flow);
    }
    // printf("V[%d][%d] step 3 completed\n", i, j);
    //  第四步，去sketch里查候选流的数据，删掉过小的
    if (r == 1)
    {
        return result;
    }
    for (vector<ans_t>::iterator item = result.begin(); item != result.end(); )
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
                if (item->bit_flow[k] == '0' && V[0][ii][jj] - V[k][ii][jj] < item->size)
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
            if (item == result.end())break;
        }
        else {
            item++;
        }
    }
    return result;
}

//从sketch中删去大流
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
            if (0 != get_bit((unsigned char*)ans, k - 1))
            {
                for (size_t i = 1; i <= r; i++)
                {
                    V[k][i][tmp_hash[i]] -= FF[it].size;
                }
            }
        }
    }
}

double MY_THETA_THRESHOLD = 1.0/(double)(1<<18);

//检查当前sketch是否符合高斯分布
//我理解的是所有L个sketch每一个都符合1sigma,2sigma,3sigma的数量要求
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
                double r = (double)V[k][i][j] / V[0][i][j];
                if (r <= p[k] + 3.0 * sigma[k] && r >= p[k] - 3.0 * sigma[k])
                    sigma_num3++;
                if (r <= p[k] + 2.0 * sigma[k] && r >= p[k] - 2.0 * sigma[k])
                    sigma_num2++;
                if (r <= p[k] + 1.0 * sigma[k] && r >= p[k] - 1.0 * sigma[k])
                    sigma_num1++;
            }
        }
        if(sigma_num1 == 0 && sigma_num2 == 0 && sigma_num3 == 0)
        {
            printf("I WONDER WHY PROGRAME RUNNING REACH HERE\n");
        }
        printf("V[%d] sigma1_num=%d,sigma2_num=%d,sigma3_num=%d\n", (int)k, (int)sigma_num1, (int)sigma_num2, (int)sigma_num3);
        double rate1 = (double)sigma_num1 / (double)(r * c);
        double rate2 = (double)sigma_num2 / (double)(r * c);
        double rate3 = (double)sigma_num3 / (double)(r * c);
        printf("V[%d] rate1=%lf,rate2=%lf,rate3=%lf\n", (int)k, rate1, rate2, rate3);
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
#endif                        // FILEOUT
    if (0 == Read_Flowdata()) //流数据读入
    {
        // 赋予哈希函数
        hash_function[1] = test_hash_0;
        hash_function[2] = test_hash_1;
        hash_function[3] = test_hash_2;
        // 流转sketch
        Flow2Sketch();
        // sketch 生成N(p,sigma)
        Sketch2N_p_sigma();

#ifdef DEBUG
        printf("Read in flow data success\n");
        memcpy(V_initial, V, sizeof(V));
        memcpy(p_initial, p, sizeof(p));
        memcpy(sigma_initial, sigma, sizeof(sigma));
#endif // DEBUG
    }
    else
    {
        printf("Read in flow data error\n");
    }

    // F为所有大流集合，FF为每次循环找出的大流集合
    vector<ans_t> F;
    double theta = 0.5;
    int nnnn = 1;
    while (1)
    {
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
                    if(my_flow_num % 100 == 0)
                    {
                        printf("CAUGHT %d FLOWS!\n",my_flow_num);
                    }
                }
            }
        }
        //本次循环找出大流时，剔除大流，重新计算期望、方差
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
                        if (strcmp(iter->bit_flow, it->bit_flow) == 0)
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
        //没有找出大流，theta减半
        if (FF.empty())
            theta /= 2;
    }

#ifdef DEBUG

    int threshold = 1000;

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
        printf("NEED %d LOOPS\n",Flow_sort.size());
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
                if(flow_size > threshold)
                    flow_queue[tmp.f].i1 = tmp.s;
                tmp_flow = *iter;
                flow_size = 1;
            }
            ++iter;
            loop_time ++;
            if(loop_time % 100000 == 0)
            {
                printf("LOOP %d TIMES!\n",loop_time);
                printf("FLOW QUEUE LENGTH: %d\n",flow_queue.size());
            }
        }
        flow_queue[tmp.f].i1 = flow_size;
        for(auto item : F)
        {
            if(item.size>threshold)
            {
                ID_input x;
                for(int i=0;i<=ID_length;i++)
                {
                    x.x[i]=item.flow[i];
                }
                flow_queue[x].i2 = item.size;
                flow_queue[x].ratio = (double)flow_queue[x].i1/flow_queue[x].i2;
            }
        }
        for (auto i : flow_queue)
        {
            //Flow_out(i.first);
            printf("appear  %d  times, SKETCH CATCH %d TIMES, RATIO: %lf\n", i.second.i1, i.second.i2, i.second.ratio);
        }
        printf("\nDEBUG END!!!\n");
    }
#endif // DEBUG
    return 0;
}